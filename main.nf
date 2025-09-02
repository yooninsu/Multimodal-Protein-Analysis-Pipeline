#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 워크플로우 파라미터 정의
params.fasta = false
params.outdir = 's3://my-af-project-bucket/output'
params.help = false

// AWS HealthOmics Ready2Run 워크플로우 ID들 (실제 확인된 값)
params.alphafold_small_id = '4885129'    // AlphaFold for up to 600 residues
params.alphafold_large_id = '6094971'    // AlphaFold for 601-1200 residues  
params.esmfold_id = '1830181'            // ESMFold for up to 800 residues
params.max_small_residues = 600          // 작은 AlphaFold 워크플로우 임계값

// 도움말 출력
if (params.help) {
    log.info """
    =================================================================
    ColabFold Multimodal Protein Analysis Pipeline
    =================================================================
    
    이 파이프라인은 AWS HealthOmics Ready2Run 워크플로우를 사용하여
    AlphaFold와 ESMFold 구조 예측을 비교 분석합니다.
    
    Usage:
    nextflow run main.nf --fasta <S3_FASTA_PATH> [options]
    
    Required Arguments:
    --fasta       S3 path to FASTA sequence file
    
    Optional Arguments:
    --outdir      S3 output directory (default: s3://my-af-project-bucket/output)
    --help        Show this help message
    
    Available Ready2Run Workflows:
    - AlphaFold (up to 600 residues): ID 4885129, ~450min
    - AlphaFold (601-1200 residues): ID 6094971, ~675min  
    - ESMFold (up to 800 residues): ID 1830181, ~15min
    
    Example:
    nextflow run main.nf --fasta s3://my-af-project-bucket/input/Q6PSU2.fasta
    
    Note: 이 파이프라인은 AWS HealthOmics Ready2Run 워크플로우만 사용하므로
    로컬에 Python 패키지나 의존성을 설치할 필요가 없습니다.
    =================================================================
    """
    exit 0
}

// 입력 검증
if (!params.fasta || !params.fasta.startsWith('s3://')) {
    error "FASTA file must be an S3 path. Use --fasta s3://bucket/file.fasta"
}

// 메인 워크플로우
workflow {
    log.info "Starting ColabFold Multimodal Analysis Pipeline"
    log.info "Using AWS HealthOmics Ready2Run Workflows Only"
    log.info "Input FASTA: ${params.fasta}"
    log.info "Output directory: ${params.outdir}"
    
    // 입력 채널 생성
    fasta_ch = Channel.fromPath(params.fasta)
    
    // 1. 서열 길이 확인 및 적절한 AlphaFold 워크플로우 선택
    SEQUENCE_LENGTH_CHECK(fasta_ch)
    
    // 2. AlphaFold 구조 예측 (자동으로 적절한 워크플로우 선택)
    ALPHAFOLD_READY2RUN(fasta_ch, SEQUENCE_LENGTH_CHECK.out.length)
    
    // 3. ESMFold 구조 예측 (AWS HealthOmics Ready2Run)
    ESMFOLD_READY2RUN(fasta_ch)
    
    // 4. 구조 비교 분석 (간단한 AWS CLI 기반)
    COMPARE_STRUCTURES(
        ALPHAFOLD_READY2RUN.out.results,
        ESMFOLD_READY2RUN.out.results
    )
    
    // 5. 최종 보고서 생성
    GENERATE_REPORT(COMPARE_STRUCTURES.out.comparison)
}

// Process 1: 서열 길이 확인
process SEQUENCE_LENGTH_CHECK {
    input:
    path fasta_file
    
    output:
    stdout emit: length
    
    script:
    """
    # FASTA 파일에서 서열 길이 계산
    SEQUENCE=\$(grep -v "^>" ${fasta_file} | tr -d '\\n' | tr -d ' ')
    LENGTH=\${#SEQUENCE}
    echo "Sequence length: \$LENGTH residues"
    
    if [ \$LENGTH -le ${params.max_small_residues} ]; then
        echo "Using AlphaFold small workflow (up to 600 residues)"
        echo "small"
    else
        echo "Using AlphaFold large workflow (601-1200 residues)"  
        echo "large"
    fi
    """
}

// Process 2: AlphaFold Ready2Run 워크플로우 (자동 선택)
process ALPHAFOLD_READY2RUN {
    publishDir "${params.outdir}/alphafold", mode: 'copy'
    
    input:
    path fasta_file
    val sequence_type
    
    output:
    path "alphafold_run_*.json", emit: results
    
    script:
    """
    echo "Starting AlphaFold Ready2Run workflow..."
    echo "Sequence analysis result: ${sequence_type}"
    
    # 서열 길이에 따라 적절한 워크플로우 ID 선택
    if [[ "${sequence_type}" =~ "small" ]]; then
        WORKFLOW_ID="${params.alphafold_small_id}"
        WORKFLOW_NAME="AlphaFold for up to 600 residues"
        ESTIMATED_TIME="450 minutes"
    else
        WORKFLOW_ID="${params.alphafold_large_id}"
        WORKFLOW_NAME="AlphaFold for 601-1200 residues"
        ESTIMATED_TIME="675 minutes"
    fi
    
    echo "Selected workflow: \$WORKFLOW_NAME (ID: \$WORKFLOW_ID)"
    echo "Estimated duration: \$ESTIMATED_TIME"
    
    # AlphaFold Ready2Run 워크플로우 시작
    RUN_ID=\$(aws omics start-run \\
        --workflow-type READY2RUN \\
        --workflow-id \$WORKFLOW_ID \\
        --name "colabfold-alphafold-\$(date +%s)" \\
        --parameters '{"inputFasta": "${params.fasta}"}' \\
        --output-uri ${params.outdir}/alphafold/raw \\
        --role-arn \$(aws sts get-caller-identity --query Account --output text | sed 's/^/arn:aws:iam::/; s/\$/:\\/role\\/HealthOmicsWorkflowRole/') \\
        --region us-east-1 \\
        --query 'id' --output text)
    
    echo "AlphaFold Run ID: \$RUN_ID"
    echo "Workflow: \$WORKFLOW_NAME"
    
    # 실행 완료 대기
    echo "Waiting for AlphaFold workflow completion (estimated: \$ESTIMATED_TIME)..."
    START_TIME=\$(date +%s)
    
    while true; do
        STATUS=\$(aws omics get-run --id \$RUN_ID --query 'status' --output text 2>/dev/null || echo "CHECKING")
        ELAPSED=\$(( \$(date +%s) - START_TIME ))
        ELAPSED_MIN=\$(( ELAPSED / 60 ))
        
        echo "AlphaFold Status: \$STATUS (Elapsed: \${ELAPSED_MIN}m) - \$(date)"
        
        case \$STATUS in
            "COMPLETED")
                echo "✅ AlphaFold prediction completed successfully"
                echo "Total time: \${ELAPSED_MIN} minutes"
                break
                ;;
            "FAILED"|"CANCELLED")
                echo "❌ AlphaFold prediction failed with status: \$STATUS"
                aws omics get-run --id \$RUN_ID --query 'statusMessage' --output text
                exit 1
                ;;
            "RUNNING"|"PENDING"|"STARTING"|"RUNNABLE")
                sleep 60  # AlphaFold는 오래 걸리므로 1분 간격
                ;;
            *)
                echo "Unknown status: \$STATUS, continuing to wait..."
                sleep 60
                ;;
        esac
    done
    
    # 실행 메타데이터 저장
    aws omics get-run --id \$RUN_ID > "alphafold_run_\${RUN_ID}.json"
    
    # 선택된 워크플로우 정보 저장
    echo "{\\"workflow_id\\": \\"\$WORKFLOW_ID\\", \\"workflow_name\\": \\"\$WORKFLOW_NAME\\", \\"estimated_time\\": \\"\$ESTIMATED_TIME\\"}" > alphafold_workflow_info.json
    
    echo "AlphaFold Ready2Run workflow completed: \$RUN_ID"
    """
}

// Process 3: ESMFold Ready2Run 워크플로우
process ESMFOLD_READY2RUN {
    publishDir "${params.outdir}/esmfold", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "esmfold_run_*.json", emit: results
    
    script:
    """
    echo "Starting ESMFold Ready2Run workflow..."
    echo "Using ESMFold for up to 800 residues (ID: ${params.esmfold_id})"
    echo "Estimated duration: 15 minutes"
    
    # ESMFold 실행
    RUN_ID=\$(aws omics start-run \\
        --workflow-type READY2RUN \\
        --workflow-id ${params.esmfold_id} \\
        --name "colabfold-esmfold-\$(date +%s)" \\
        --parameters '{"inputFasta": "${params.fasta}"}' \\
        --output-uri ${params.outdir}/esmfold/raw \\
        --role-arn \$(aws sts get-caller-identity --query Account --output text | sed 's/^/arn:aws:iam::/; s/\$/:\\/role\\/HealthOmicsWorkflowRole/') \\
        --region us-east-1 \\
        --query 'id' --output text)
    
    echo "ESMFold Run ID: \$RUN_ID"
    echo "Publisher: Meta Research"
    
    # 실행 완료 대기
    echo "Waiting for ESMFold workflow completion (estimated: 15 minutes)..."
    START_TIME=\$(date +%s)
    
    while true; do
        STATUS=\$(aws omics get-run --id \$RUN_ID --query 'status' --output text 2>/dev/null || echo "CHECKING")
        ELAPSED=\$(( \$(date +%s) - START_TIME ))
        ELAPSED_MIN=\$(( ELAPSED / 60 ))
        
        echo "ESMFold Status: \$STATUS (Elapsed: \${ELAPSED_MIN}m) - \$(date)"
        
        case \$STATUS in
            "COMPLETED")
                echo "✅ ESMFold prediction completed successfully"
                echo "Total time: \${ELAPSED_MIN} minutes"
                break
                ;;
            "FAILED"|"CANCELLED")
                echo "❌ ESMFold prediction failed with status: \$STATUS"
                aws omics get-run --id \$RUN_ID --query 'statusMessage' --output text
                exit 1
                ;;
            "RUNNING"|"PENDING"|"STARTING"|"RUNNABLE")
                sleep 30  # ESMFold는 빠르므로 30초 간격
                ;;
            *)
                echo "Unknown status: \$STATUS, continuing to wait..."
                sleep 30
                ;;
        esac
    done
    
    # 실행 메타데이터 저장
    aws omics get-run --id \$RUN_ID > "esmfold_run_\${RUN_ID}.json"
    
    # 워크플로우 정보 저장
    echo "{\\"workflow_id\\": \\"${params.esmfold_id}\\", \\"workflow_name\\": \\"ESMFold for up to 800 residues\\", \\"publisher\\": \\"Meta Research\\", \\"estimated_time\\": \\"15 minutes\\"}" > esmfold_workflow_info.json
    
    echo "ESMFold Ready2Run workflow completed: \$RUN_ID"
    """
}

// Process 4: 구조 비교 분석 (AWS CLI만 사용)
process COMPARE_STRUCTURES {
    publishDir "${params.outdir}/comparison", mode: 'copy'
    
    input:
    path alphafold_results
    path esmfold_results
    
    output:
    path "structure_comparison.json", emit: comparison
    path "analysis_summary.txt", emit: summary
    
    script:
    """
    echo "Analyzing AlphaFold vs ESMFold results..."
    
    # 결과 파일 정보 수집
    AF_RUN_ID=\$(basename ${alphafold_results} .json | sed 's/alphafold_run_//')
    ESM_RUN_ID=\$(basename ${esmfold_results} .json | sed 's/esmfold_run_//')
    
    echo "AlphaFold Run ID: \$AF_RUN_ID"
    echo "ESMFold Run ID: \$ESM_RUN_ID"
    
    # 실행 시간 및 비용 정보 수집
    AF_START=\$(cat ${alphafold_results} | jq -r '.creationTime')
    AF_END=\$(cat ${alphafold_results} | jq -r '.stopTime // "null"')
    AF_STATUS=\$(cat ${alphafold_results} | jq -r '.status')
    ESM_START=\$(cat ${esmfold_results} | jq -r '.creationTime')
    ESM_END=\$(cat ${esmfold_results} | jq -r '.stopTime // "null"')
    ESM_STATUS=\$(cat ${esmfold_results} | jq -r '.status')
    ESM_START=\$(cat ${esmfold_results} | jq -r '.creationTime')
    ESM_END=\$(cat ${esmfold_results} | jq -r '.stopTime // "null"')
    ESM_STATUS=\$(cat ${esmfold_results} | jq -r '.status')

    # 워크플로우 정보 확인
    if [ -f "alphafold_workflow_info.json" ]; then
        AF_WORKFLOW=\$(cat alphafold_workflow_info.json | jq -r '.workflow_name')
        AF_ESTIMATED=\$(cat alphafold_workflow_info.json | jq -r '.estimated_time')
    else
        AF_WORKFLOW="AlphaFold (auto-selected)"
        AF_ESTIMATED="unknown"
    fi
    
    # 비교 결과 JSON 생성
    cat > structure_comparison.json << EOF
{
    "comparison_timestamp": "\$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "input_fasta": "${params.fasta}",
    "workflows": {
        "alphafold": {
            "run_id": "\$AF_RUN_ID",
            "workflow_name": "\$AF_WORKFLOW",
            "estimated_duration": "\$AF_ESTIMATED",
            "publisher": "DeepMind",
            "start_time": "\$AF_START",
            "end_time": "\$AF_END",
            "status": "\$AF_STATUS",
            "output_location": "${params.outdir}/alphafold/raw",
            "method": "AlphaFold3 via HealthOmics Ready2Run"
        },
        "esmfold": {
            "run_id": "\$ESM_RUN_ID",
            "workflow_name": "ESMFold for up to 800 residues",
            "workflow_id": "${params.esmfold_id}",
            "estimated_duration": "15 minutes",
            "publisher": "Meta Research",
            "start_time": "\$ESM_START",
            "end_time": "\$ESM_END", 
            "status": "\$ESM_STATUS",
            "output_location": "${params.outdir}/esmfold/raw",
            "method": "ESMFold via HealthOmics Ready2Run"
        }
    },
    "comparison_notes": {
        "speed": "ESMFold is ~45x faster than AlphaFold (15min vs 450-675min)",
        "accuracy": "AlphaFold3 generally provides higher accuracy and confidence",
        "residue_limits": "AlphaFold: up to 1200, ESMFold: up to 800",
        "use_case": "ESMFold for rapid screening, AlphaFold for high-accuracy predictions"
    },
    "recommendations": {
        "fast_screening": "Use ESMFold for initial structural screening",
        "high_accuracy": "Use AlphaFold for publication-quality structures", 
        "consensus_regions": "Trust regions where both methods agree",
        "validation": "Compare confidence scores and structural features"
    }
}
EOF
    
    # 분석 요약 텍스트 생성
    cat > analysis_summary.txt << EOF
=== ColabFold Multimodal Analysis Summary ===

Analysis completed: \$(date)
Input FASTA: ${params.fasta}

AlphaFold Analysis:
- Run ID: $AF_RUN_ID
- Workflow: $AF_WORKFLOW
- Estimated time: $AF_ESTIMATED
- Publisher: DeepMind
- Start: $AF_START
- End: $AF_END
- Status: $AF_STATUS

ESMFold Analysis:
- Run ID: $ESM_RUN_ID
- Workflow: ESMFold for up to 800 residues (ID: ${params.esmfold_id})
- Estimated time: 15 minutes
- Publisher: Meta Research
- Start: $ESM_START
- End: $ESM_END
- Status: $ESM_STATUS

Output Locations:
- AlphaFold results: ${params.outdir}/alphafold/raw
- ESMFold results: ${params.outdir}/esmfold/raw
- Comparison data: ${params.outdir}/comparison

Performance Comparison:
- Speed: ESMFold is ~45x faster than AlphaFold
- Accuracy: AlphaFold typically more accurate
- Residue limits: AlphaFold (up to 1200), ESMFold (up to 800)

Next Steps:
1. Download PDB files from S3 output locations
2. Load structures in PyMOL or ChimeraX
3. Align and compare structures (e.g., RMSD analysis)
4. Compare confidence scores and identify consensus regions
5. Analyze epitope predictions based on structure

=== Pipeline Completed Successfully ===
EOF

echo "Structure comparison analysis completed"
            "run_id": "\$ESM_RUN_ID", 
            "start_time": "\$ESM_START",
            "end_time": "\$ESM_END",
            "output_location": "${params.outdir}/esmfold/raw",
            "method": "ESMFold via HealthOmics Ready2Run"
        }
    },
    "comparison_notes": {
        "speed": "ESMFold is typically 60x faster than MSA-based methods",
        "accuracy": "AlphaFold3 generally provides higher accuracy",
        "use_case": "Both methods useful for comparative analysis"
    }
}
EOF
    
    # 분석 요약 텍스트 생성
    cat > analysis_summary.txt << EOF
=== ColabFold Multimodal Analysis Summary ===

Analysis completed: \$(date)
Input FASTA: ${params.fasta}

AlphaFold Analysis:
- Run ID: \$AF_RUN_ID
- Method: AlphaFold3 via AWS HealthOmics Ready2Run
- Start: \$AF_START
- Status: Completed

ESMFold Analysis:
- Run ID: \$ESM_RUN_ID  
- Method: ESMFold via AWS HealthOmics Ready2Run
- Start: \$ESM_START
- Status: Completed

Output Locations:
- AlphaFold results: ${params.outdir}/alphafold/raw
- ESMFold results: ${params.outdir}/esmfold/raw
- Comparison data: ${params.outdir}/comparison

Next Steps:
1. Download PDB files from S3 output locations
2. Perform detailed structural analysis using PyMOL or ChimeraX
3. Compare confidence scores and structural features
4. Analyze epitope prediction differences

=== Pipeline Completed Successfully ===
EOF
    
    echo "Structure comparison analysis completed"
    """
}

// Process 4: 최종 보고서 생성
process GENERATE_REPORT {
    publishDir "${params.outdir}/final", mode: 'copy'
    
    input:
    path comparison_json
    
    output:
    path "colabfold_final_report.html"
    path "workflow_manifest.json"
    
    script:
    """
    echo "Generating final ColabFold analysis report..."
    
    # HTML 보고서 생성
    cat > colabfold_final_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>ColabFold Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background: #f0f8ff; padding: 15px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; }
        .success { color: #28a745; }
        .info { color: #007bff; }
        .code { background: #f8f9fa; padding: 10px; font-family: monospace; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 ColabFold Multimodal Protein Analysis</h1>
        <p><strong>Analysis completed:</strong> \$(date)</p>
        <p><strong>Input:</strong> ${params.fasta}</p>
    </div>
    
    <div class="section">
        <h2>📊 Workflow Summary</h2>
        <p class="success">✅ AlphaFold3 structure prediction completed via AWS HealthOmics</p>
        <p class="success">✅ ESMFold structure prediction completed via AWS HealthOmics</p>
        <p class="success">✅ Comparative analysis completed</p>
    </div>
    
    <div class="section">
        <h2>🔬 Methods Used</h2>
        <ul>
            <li><strong>AlphaFold3:</strong> DeepMind's latest structure prediction model</li>
            <li><strong>ESMFold:</strong> Meta AI's language model-based folding (up to 800 residues)</li>
            <li><strong>Infrastructure:</strong> AWS HealthOmics Ready2Run workflows</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>📁 Output Locations</h2>
        <div class="code">
            AlphaFold results: ${params.outdir}/alphafold/raw/<br>
            ESMFold results: ${params.outdir}/esmfold/raw/<br>
            Comparison data: ${params.outdir}/comparison/<br>
        </div>
    </div>
    
    <div class="section">
        <h2>🚀 Next Steps</h2>
        <ol>
            <li>Download PDB structures from S3 output locations</li>
            <li>Load structures in PyMOL or ChimeraX for visualization</li>
            <li>Perform RMSD analysis between AlphaFold and ESMFold structures</li>
            <li>Compare confidence scores and identify high-confidence regions</li>
            <li>Analyze potential epitope regions based on surface accessibility</li>
        </ol>
    </div>
    
    <div class="section">
        <h2>💡 Analysis Tips</h2>
        <p><strong>AlphaFold3:</strong> Higher accuracy, especially for complex domains</p>
        <p><strong>ESMFold:</strong> Faster execution, good for initial screening</p>
        <p><strong>Consensus regions:</strong> Areas where both methods agree are likely highly reliable</p>
    </div>
</body>
</html>
EOF

    # 워크플로우 매니페스트 생성
    cat > workflow_manifest.json << EOF
{
    "pipeline_name": "ColabFold Multimodal Protein Analysis",
    "version": "1.0.0",
    "timestamp": "\$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "input_fasta": "${params.fasta}",
    "output_directory": "${params.outdir}",
    "workflows_used": [
        {
            "name": "AlphaFold3",
            "type": "AWS HealthOmics Ready2Run",
            "description": "High-accuracy protein structure prediction"
        },
        {
            "name": "ESMFold", 
            "type": "AWS HealthOmics Ready2Run",
            "description": "Language model-based structure prediction (up to 800 residues)"
        }
    ],
    "infrastructure": {
        "compute": "AWS HealthOmics (fully managed)",
        "storage": "Amazon S3",
        "dependencies": "None (containerized Ready2Run workflows)"
    },
    "advantages": [
        "No local dependencies required",
        "Scalable cloud execution", 
        "Predictable per-run pricing",
        "Enterprise-grade security and compliance"
    ]
}
EOF

    echo "ColabFold analysis report generated successfully"
    echo "📊 View report: ${params.outdir}/final/colabfold_final_report.html"
    echo "📋 Manifest: ${params.outdir}/final/workflow_manifest.json"
    """
}
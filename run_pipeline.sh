#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "사용법:"
    echo "  predict  - AlphaFold 구조 예측 실행"
    echo "  analyze  - 에피토프 분석 실행 (PDB 파일 필요)"
    echo "  full     - 전체 파이프라인 실행"
    echo "  bash     - 대화형 모드"
    exit 1
fi

case "$1" in
    "predict")
        echo "=== AlphaFold 구조 예측 시작 ==="
        if ls /app/input/*.fasta 1> /dev/null 2>&1; then
            colabfold_batch /app/input/*.fasta /app/output --num-models 1
        else
            echo "오류: input 디렉토리에 FASTA 파일이 없습니다."
            exit 1
        fi
        ;;
    "analyze")
        echo "=== 에피토프 분석 시작 ==="
        PDB_FILE=$(ls /app/output/*.pdb 2>/dev/null | head -1)
        CSV_FILE=$(ls /app/input/*.csv 2>/dev/null | head -1)
        
        if [ -z "$PDB_FILE" ]; then
            echo "오류: output 디렉토리에 PDB 파일이 없습니다."
            exit 1
        fi
        
        if [ -z "$CSV_FILE" ]; then
            echo "오류: input 디렉토리에 CSV 파일이 없습니다."
            exit 1
        fi
        
        python /app/analyze_structure.py \
            --pdb_file "$PDB_FILE" \
            --truth_file "$CSV_FILE" \
            --output_file /app/output/epitope_analysis_results.csv
        ;;
    "full")
        echo "=== 전체 파이프라인 시작 ==="
        /app/run_pipeline.sh predict
        if [ $? -eq 0 ]; then
            /app/run_pipeline.sh analyze
        else
            echo "구조 예측 실패로 인해 분석을 중단합니다."
            exit 1
        fi
        ;;
    "bash")
        exec /bin/bash
        ;;
    *)
        echo "알 수 없는 명령어: $1"
        exit 1
        ;;
esac

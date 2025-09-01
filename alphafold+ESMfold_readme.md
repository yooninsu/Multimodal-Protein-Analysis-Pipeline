# GraphBepi Multimodal Protein Analysis Pipeline

A comprehensive Nextflow pipeline for multimodal protein epitope prediction combining structure prediction, sequence embeddings, and graph neural networks on AWS cloud infrastructure.

## 🧬 Overview

This pipeline integrates multiple state-of-the-art approaches for protein analysis:
- **AlphaFold3** structure prediction via AWS HealthOmics Ready2Run workflows
- **ESMFold** comparative structure prediction
- **ESM-2** protein language model embeddings
- **Graph Neural Networks** for epitope prediction through GraphBepi
- **Multimodal fusion** of structural and sequence features

## 🏗️ Architecture

```
Input FASTA → [AlphaFold3] → Structure Features
             ↓
             [ESMFold] → Comparative Analysis
             ↓
             [ESM-2] → Sequence Embeddings
             ↓
             [Graph Construction] → Contact Graph
             ↓
             [GraphBepi GNN] → Epitope Predictions
             ↓
             [Result Integration] → Final Analysis
```

## 📋 Prerequisites

### AWS Infrastructure
- **AWS Account** with HealthOmics access
- **S3 Bucket** for input/output storage
- **ECR Repository** for custom containers
- **IAM Roles** with appropriate permissions:
  - `omics:StartRun`, `omics:GetRun`
  - `s3:GetObject`, `s3:PutObject`
  - `batch:SubmitJob` (if using AWS Batch)

### Software Requirements
- [Nextflow](https://www.nextflow.io/) >= 22.04.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)
- AWS CLI configured with appropriate credentials

### Container Dependencies
The pipeline uses a custom container with:
- **PyTorch** >= 1.12.0
- **ESM** (fair-esm) for protein language models
- **BioPython** for structural analysis
- **ESMFold** for comparative structure prediction
- **NumPy**, **Pandas** for data processing

## 🚀 Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/your-username/graphbepi-pipeline.git
cd graphbepi-pipeline
```

### 2. Configure AWS Credentials
```bash
aws configure
# or use IAM roles/instance profiles
```

### 3. Prepare Input Data
Upload your FASTA file to S3:
```bash
aws s3 cp your_protein.fasta s3://your-bucket/input/
```

### 4. Run Pipeline
```bash
nextflow run main.nf \
  --fasta s3://your-bucket/input/your_protein.fasta \
  --outdir s3://your-bucket/output
```

## 📊 Pipeline Processes

### Process 1: AlphaFold3 Structure Prediction
- **Method**: AWS HealthOmics Ready2Run workflow (ID: 6094971)
- **Output**: High-confidence 3D protein structure
- **Runtime**: ~10-60 minutes depending on sequence length

### Process 2: ESMFold Comparative Prediction
- **Method**: Meta's ESMFold transformer-based folding
- **Purpose**: Structural comparison and validation
- **Container**: Custom ECR image with ESMFold dependencies

### Process 3: ESM-2 Sequence Embeddings
- **Model**: ESM-2 650M parameter model (esm2_t33_650M_UR50D)
- **Output**: 1280-dimensional per-residue embeddings
- **Features**: Contextual protein language representations

### Process 4: Structural Feature Extraction
- **SASA**: Solvent Accessible Surface Area calculation
- **B-factors**: Structural flexibility analysis
- **Comparison**: AlphaFold vs ESMFold structural differences

### Process 5: Protein Graph Construction
- **Method**: Contact-based graph generation
- **Cutoff**: 8Å distance threshold for residue contacts
- **Output**: Adjacency matrix and graph topology

### Process 6: GraphBepi Prediction
- **Architecture**: Graph Neural Network with multimodal fusion
- **Features**: Structure + Sequence + Graph topology
- **Output**: Per-residue epitope probabilities

### Process 7: Result Integration
- **Analysis**: Comprehensive prediction summary
- **Metrics**: Confidence scores and performance statistics

## 📁 Output Structure

```
output/
├── alphafold/
│   ├── alphafold_structure.pdb
│   └── alphafold_confidence.json
├── esmfold/
│   ├── esmfold_structure.pdb
│   └── esmfold_confidence.json
├── embeddings/
│   ├── esm2_embeddings.npy
│   └── sequence_features.csv
├── features/
│   ├── structure_features.csv
│   └── structure_comparison.csv
├── graph/
│   ├── protein_graph.json
│   └── adjacency_matrix.npy
├── predictions/
│   ├── epitope_predictions.csv
│   └── confidence_scores.csv
└── final/
    ├── final_analysis_report.csv
    └── summary_statistics.txt
```

## 🔧 Configuration

### Pipeline Parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--fasta` | S3 path to input FASTA file | - | ✅ |
| `--outdir` | S3 output directory | `s3://my-af-project-bucket/output` | ❌ |
| `--help` | Display help message | `false` | ❌ |

### Advanced Configuration

Create `nextflow.config` for custom settings:

```groovy
// nextflow.config
aws {
    region = 'us-east-1'
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
    }
}

process {
    executor = 'awsbatch'
    queue = 'your-batch-queue'
    
    withName:ALPHAFOLD_PREDICTION {
        memory = '16 GB'
        cpus = 4
    }
    
    withName:ESM2_EMBEDDING {
        memory = '32 GB'
        cpus = 8
        accelerator = 1, type: 'nvidia-tesla-v100'
    }
}
```

## 🧪 Example Usage

### Basic Run
```bash
nextflow run main.nf --fasta s3://my-bucket/Q6PSU2.fasta
```

### With Custom Output Directory
```bash
nextflow run main.nf \
  --fasta s3://my-bucket/input/protein.fasta \
  --outdir s3://my-results-bucket/analysis-2024
```

### Local Testing (with Singularity)
```bash
nextflow run main.nf \
  --fasta /path/to/local/protein.fasta \
  --outdir ./results \
  -profile singularity
```

## 📈 Performance Considerations

### Resource Requirements
- **Memory**: 32GB+ for ESM-2 embeddings
- **CPU**: 8+ cores recommended
- **GPU**: Optional, accelerates ESM-2 and GNN processing
- **Storage**: ~5GB per protein analysis

### Optimization Tips
1. **Batch Processing**: Use `Channel.fromPath()` with glob patterns for multiple proteins
2. **Caching**: Enable Nextflow resume (`-resume`) for interrupted runs
3. **GPU Acceleration**: Use GPU instances for ESM-2 and GNN processes
4. **Parallelization**: Adjust `process.cpus` based on available resources

## 🔬 Scientific Background

### AlphaFold3 Integration
This pipeline leverages AWS HealthOmics' AlphaFold3 Ready2Run workflow, providing:
- **State-of-the-art accuracy** for protein structure prediction
- **Confidence scoring** for prediction reliability assessment
- **Scalable cloud execution** without local GPU requirements

### Multimodal Fusion Approach
The GraphBepi framework combines:
1. **Structural Features**: 3D geometry, surface accessibility, flexibility
2. **Sequence Embeddings**: Evolutionary and physicochemical information
3. **Graph Topology**: Residue contact networks and spatial relationships

### Graph Neural Networks for Proteins
- **Node Features**: Per-residue embeddings and structural properties
- **Edge Features**: Spatial proximity and interaction strength
- **Architecture**: Message passing between spatially connected residues

## 🛠️ Troubleshooting

### Common Issues

**1. AWS Permissions Error**
```bash
# Check IAM permissions
aws sts get-caller-identity
aws omics list-workflows --region us-east-1
```

**2. Container Access Issues**
```bash
# Login to ECR
aws ecr get-login-password --region ap-southeast-2 | \
docker login --username AWS --password-stdin 891612562910.dkr.ecr.ap-southeast-2.amazonaws.com
```

**3. Memory Issues with ESM-2**
- Reduce batch size or use model checkpointing
- Consider using `esm2_t12_35M_UR50D` for smaller sequences

**4. S3 Access Problems**
```bash
# Test S3 access
aws s3 ls s3://your-bucket/
aws s3 cp test.txt s3://your-bucket/test.txt
```

## 📚 References

1. **AlphaFold3**: [Accurate structure prediction of biomolecular interactions](https://www.nature.com/articles/s41586-024-07487-w)
2. **ESM-2**: [Language models of protein sequences at the scale of evolution](https://www.science.org/doi/10.1126/science.ade2574)
3. **ESMFold**: [Evolutionary-scale prediction of atomic level protein structure](https://www.science.org/doi/10.1126/science.ade2574)
4. **GraphBepi**: [Graph neural networks for epitope prediction](https://academic.oup.com/bioinformatics)
5. **AWS HealthOmics**: [Ready2Run Workflows Documentation](https://docs.aws.amazon.com/omics/latest/dev/workflows.html)

## 🤝 Contributing

1. **Fork the repository**
2. **Create feature branch**: `git checkout -b feature/amazing-feature`
3. **Commit changes**: `git commit -m 'Add amazing feature'`
4. **Push to branch**: `git push origin feature/amazing-feature`
5. **Open Pull Request**

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **AWS HealthOmics** team for AlphaFold3 Ready2Run workflows
- **Meta AI** for ESM-2 and ESMFold models
- **Nextflow** community for workflow management framework
- **GraphBepi** authors for the epitope prediction methodology

## 📧 Contact

- **Author**: Your Name
- **Email**: your.email@institution.edu
- **Lab**: [Your Lab Website](https://yourlab.edu)
- **Issues**: [GitHub Issues](https://github.com/your-username/graphbepi-pipeline/issues)

---

**Citation**: If you use this pipeline in your research, please cite:
```bibtex
@software{graphbepi_pipeline,
  title={GraphBepi Multimodal Protein Analysis Pipeline},
  author={Your Name},
  year={2024},
  url={https://github.com/your-username/graphbepi-pipeline}
}
```
# MetaMAG Explorer

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://msatti123.github.io/MetaMAG_Explorer/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Pipeline](https://img.shields.io/badge/pipeline-v1.0.0-orange)](https://github.com/msatti123/MetaMAG_Explorer)
[![HPC Ready](https://img.shields.io/badge/HPC-SLURM-red)](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html)
[![Python](https://img.shields.io/badge/python-3.8%20%7C%203.9-blue)](https://www.python.org/)

## 📚 **Complete Documentation Available**

### **[📖 View Full Documentation](https://msatti123.github.io/MetaMAG_Explorer/)** | **[🚀 Installation Guide](https://msatti123.github.io/MetaMAG_Explorer/installation.html)** | **[📋 User Guide](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html)**

---

## Overview

**MetaMAG Explorer** is a comprehensive, modular pipeline for **Metagenome-Assembled Genome (MAG)** analysis, specializing in novel MAG detection and database integration. The pipeline provides an end-to-end solution from raw metagenomic reads to functionally annotated genomes with phylogenetic placement.

### 🎯 Key Highlights
- **Novel MAG Detection**: Automatically identifies and processes novel genomes not present in existing databases
- **Database Integration**: Seamlessly integrates novel MAGs into Kraken2 databases for improved future analyses
- **Rumen Microbiome Specialization**: Includes dedicated workflows for rumen microbiome studies with reference MAG integration
- **Modular Design**: Run individual steps or the complete pipeline based on your needs
- **HPC Optimized**: Built for high-performance computing environments with SLURM integration

## 📖 Documentation

**For detailed instructions, please visit our comprehensive documentation:**

| Documentation | Description |
|--------------|-------------|
| **[🏠 Main Documentation](https://msatti123.github.io/MetaMAG_Explorer/)** | Complete pipeline overview, features, and FAQ |
| **[⚙️ Installation Guide](https://msatti123.github.io/MetaMAG_Explorer/installation.html)** | Step-by-step installation instructions |
| **[📘 User Guide](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html)** | Detailed execution guide for all pipeline steps |

## 🚀 Quick Start

### Prerequisites
- Linux OS (Ubuntu recommended)
- Conda/Miniconda
- Python 3.8 or 3.9
- 16GB+ RAM (32GB+ recommended)
- 100GB+ storage space

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/msatti123/MetaMAG_Explorer.git
cd MetaMAG_Explorer

# 2. Install tools via setup script
python setup_tools.py --all  # Install all tools
# OR
python setup_tools.py --tools fastqc fastp  # Install specific tools

# 3. Verify installation
python setup_tools.py --verify

# 4. Install the package
pip install -e .
```

**[→ See detailed installation instructions](https://msatti123.github.io/MetaMAG_Explorer/installation.html)**

### Basic Usage

```bash
# Create configuration file (only 3 lines needed!)
cat > project_config.yaml << EOF
input_dir: "/path/to/raw/reads"
output_dir: "/path/to/output"
reference: "/path/to/host/genome.fa"  # Optional
EOF

# Create sample list
echo "Sample_001" > samples.txt
echo "Sample_002" >> samples.txt

# Run quality control
python -m MetaMAG.main \
  --project_config project_config.yaml \
  --samples-file samples.txt \
  --steps qc \
  --batch_size 20 \
  --cpus 4 \
  --memory "8G" \
  --time "2:00:00"
```

**[→ See complete user guide with all steps](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html)**

## 🔧 Pipeline Components

### Core Modules

1. **Quality Control & Preprocessing**
   - FastQC quality assessment
   - Fastp trimming
   - Host genome removal

2. **Assembly**
   - Single-sample assembly (IDBA-UD)
   - Co-assembly (MEGAHIT)
   - Assembly evaluation (MetaQUAST)

3. **Binning & Refinement**
   - Multiple binning algorithms (MetaBAT2, MaxBin2, CONCOCT)
   - Bin refinement (DAS Tool)
   - Quality assessment (CheckM)
   - Dereplication (dRep)

4. **Taxonomic Classification**
   - GTDB-Tk classification
   - Novel MAG identification
   - Kraken2 database integration

5. **Functional Annotation**
   - EggNOG-mapper annotation
   - dbCAN CAZyme identification
   - KEGG pathway analysis
   - Advanced visualizations

6. **Abundance & Phylogeny**
   - Kraken2/Bracken abundance estimation
   - Phylogenetic tree construction
   - Publication-quality visualizations

**[→ View detailed pipeline workflow](https://msatti123.github.io/MetaMAG_Explorer/#pipeline-overview)**

## 🐄 Rumen Microbiome Features

Special workflows for rumen microbiome analysis:
- Automatic download of reference MAGs (RUG, RMGMC, MGnify)
- Rumen-specific MAG dereplication
- Integration with established rumen reference genomes
- Novel rumen MAG detection and cataloging

**[→ Learn about rumen-specific features](https://msatti123.github.io/MetaMAG_Explorer/#rumen-mags)**

## 📊 Key Outputs

```
output_dir/
├── QC/                      # Quality control reports
├── Assembly/                # Assembled contigs
├── Bin_Refinement/         
│   └── drep/               # Dereplicated MAGs
├── Novel_Mags/             
│   ├── gtdbtk/             # Taxonomic classifications
│   └── true_novel_MAGs/    # Identified novel genomes
├── Functional_Annotation/   
│   ├── eggNOG/             # Functional annotations
│   └── dbCAN/              # CAZyme predictions
├── Kraken_Database/        # Updated database with novel MAGs
└── Abundance/              # MAG abundance profiles
```

## 🛠️ Available Pipeline Steps

Run individual steps or combinations:

```bash
# Single step
--steps qc

# Multiple steps
--steps trimming host_removal assembly

# Complete phases
--steps qc trimming host_removal single_assembly single_binning
```

**Available steps:**
`qc`, `multiqc`, `trimming`, `host_removal`, `single_assembly`, `co_assembly`, `single_binning`, `single_bin_refinement`, `evaluation`, `dRep`, `gtdbtk`, `identify_novel_mags`, `process_novel_mags`, `build_kraken_db`, `eggnog_annotation`, `dbcan_annotation`, `functional_analysis`, `advanced_visualizations`, `abundance_estimation`, `mags_tree`, `tree_visualization`

**[→ See step-by-step execution guide](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html)**

## 🖥️ HPC/SLURM Integration

Example SLURM submission:

```bash
#!/bin/bash
#SBATCH -p partition_name
#SBATCH --mem=32G
#SBATCH --cpus=16
#SBATCH -t 24:00:00

source /path/to/conda/activate metamag
python -m MetaMAG.main \
  --project_config config.yaml \
  --samples-file samples.txt \
  --steps qc trimming assembly
```

**[→ View complete SLURM examples](https://msatti123.github.io/MetaMAG_Explorer/user-guide.html#slurm-setup)**

## 📋 Requirements

### Computational Resources (Recommended)
- **Assembly**: 100-200 GB RAM
- **Binning**: 32-64 GB RAM  
- **GTDB-Tk**: 200 GB RAM
- **Storage**: 1-5 TB for typical projects

### Core Dependencies
- FastQC, Fastp, BWA, Samtools
- IDBA-UD, MEGAHIT, MetaQUAST
- MetaBAT2, MaxBin2, CONCOCT, DAS Tool
- CheckM, dRep, GTDB-Tk
- Kraken2, Bracken
- EggNOG-mapper, dbCAN, Prodigal

**[→ View complete tool list and versions](https://msatti123.github.io/MetaMAG_Explorer/installation.html)**

## 📚 Documentation & Support

- **[📖 Full Documentation](https://msatti123.github.io/MetaMAG_Explorer/)** - Complete pipeline documentation
- **[❓ FAQ](https://msatti123.github.io/MetaMAG_Explorer/#faq)** - Frequently asked questions
- **[🔧 Troubleshooting](https://msatti123.github.io/MetaMAG_Explorer/#troubleshooting)** - Common issues and solutions
- **[📧 Contact](mailto:mariasatti12@gmail.com)** - mariasatti12@gmail.com

## 📄 Citation

**Manuscript in preparation**

If you use MetaMAG Explorer in your research, please cite:

```
Satti et al. (2025) MetaMAG Explorer: A comprehensive pipeline for novel MAG 
discovery and metagenomic profiling. [Manuscript in preparation]
```

For now, please cite the repository:
```
MetaMAG Explorer (2025). GitHub repository. 
https://github.com/msatti123/MetaMAG_Explorer
```

## 🤝 Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/NewFeature`)
3. Commit your changes (`git commit -m 'Add NewFeature'`)
4. Push to the branch (`git push origin feature/NewFeature`)
5. Open a Pull Request

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Tools and database developers cited in our [references](https://msatti123.github.io/MetaMAG_Explorer/#references)
- GTDB team for the Genome Taxonomy Database
- Rumen microbiome research community

## 🚦 Project Status

- ✅ **Version 1.0.0** - Stable release
- 🔄 Active development and maintenance
- 📝 Manuscript in preparation

---

<p align="center">
  <b>For comprehensive documentation and detailed instructions, please visit:</b><br>
  <a href="https://msatti123.github.io/MetaMAG_Explorer/">https://msatti123.github.io/MetaMAG_Explorer/</a>
</p>

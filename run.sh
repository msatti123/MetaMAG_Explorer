#!/bin/bash
# ============================================================================
# MetaMAG Explorer - SLURM Submission Script
# Version: 1.0.0
# Description: Template script for running MetaMAG Explorer pipeline on HPC
# ============================================================================

# ============================================================================
# SLURM CONFIGURATION
# Adjust these parameters based on your HPC cluster and job requirements
# ============================================================================
#SBATCH -p normal                     # Partition/queue name (e.g., normal, gpu, highmem)
#SBATCH -N 1                          # Number of nodes
#SBATCH -n 1                          # Number of tasks
#SBATCH --mem=20G                     # Memory allocation (adjust per step)
#SBATCH -t 24:00:00                   # Time limit (HH:MM:SS or D-HH:MM:SS)
#SBATCH --output=metamag_%j.out       # Standard output log (%j = job ID)
#SBATCH --error=metamag_%j.err        # Standard error log
#SBATCH --job-name=MetaMAG            # Job name for queue display

# ============================================================================
# USER CONFIGURATION - MODIFY THESE PATHS
# ============================================================================

# [REQUIRED] Base directories
WORK_DIR="/path/to/MetaMAG-1.0.0"                    # MetaMAG installation directory
OUTPUT_BASE="/path/to/output/directory"               # Where results will be stored
CONDA_ENV="/path/to/miniconda3/bin/activate metamag"  # Full path to conda activate + environment

# [REQUIRED] Input files
SAMPLES="/path/to/samples.txt"                        # List of sample names (one per line)
PROJECT_CONFIG="/path/to/project_config.yaml"         # Project configuration file
LOG_DIR="./logs/pipeline_run"                         # Directory for pipeline logs

# ============================================================================
# STEP SELECTION - CHOOSE WHICH STEP(S) TO RUN
# ============================================================================
# Available steps (uncomment/modify as needed):
# - qc                        : Quality control with FastQC
# - multiqc                   : Aggregate QC reports
# - trimming                  : Read trimming with fastp
# - host_removal              : Remove host contamination
# - single_assembly           : Individual sample assembly (IDBA-UD)
# - co_assembly              : Co-assembly of all samples (MEGAHIT)
# - metaquast                : Assembly quality assessment
# - single_binning           : Bin contigs into MAGs
# - single_bin_refinement    : Refine bins with DAS Tool
# - evaluation               : Evaluate MAG quality with CheckM
# - dRep                     : Dereplicate MAGs
# - gtdbtk                   : Taxonomic classification
# - identify_novel_mags      : Find novel MAGs
# - process_novel_mags       : Process and add novel MAGs to database
# - eggnog_annotation        : Functional annotation with eggNOG
# - dbcan_annotation         : CAZyme annotation
# - functional_analysis      : Integrated functional analysis
# - advanced_visualizations  : Generate comprehensive visualizations
# - kraken_abundance         : Taxonomic profiling
# - abundance_estimation     : Complete abundance analysis with Bracken
# - mags_tree               : Build phylogenetic tree
# - tree_visualization      : Visualize phylogenetic tree

STEPS="qc"  # <-- MODIFY THIS: Select step(s) to run (space-separated for multiple)

# ============================================================================
# RESOURCE ALLOCATION
# Adjust based on step requirements and cluster capabilities
# ============================================================================
BATCH_SIZE=16         # Number of samples to process in parallel
CPUS=8                # Number of CPU cores per job
MEMORY="20G"          # Memory per job (e.g., 8G, 32G, 100G)
TIME="24:00:00"       # Time limit per job (adjust based on step)


# ============================================================================
# ASSEMBLY PARAMETERS
# ============================================================================
ASSEMBLY_TYPE="both"           # Options: idba, megahit, both
BINNING_METHODS="metabat2 maxbin2 concoct"  # Space-separated list

# ============================================================================
# REFERENCE GENOMES & DATABASES
# ============================================================================
# Host genome for decontamination (optional, for host-associated samples)
REFERENCE_GENOME="/path/to/host/reference/genome.fna"

# Kraken2 database for taxonomic profiling
KRAKEN_DB_PATH="/path/to/kraken2/database"

# KEGG database file for functional analysis (optional but recommended)
KEGG_DB_FILE="/path/to/kegg/ko.txt"

# ============================================================================
# RUMEN-SPECIFIC CONFIGURATION (Optional - for rumen microbiome projects)
# ============================================================================
# Uncomment and set these if processing rumen microbiome data

# RUMENREF_DATASET="/path/to/rumen/reference/mags"           # Downloaded reference MAGs
# RUMENREF_DREP="/path/to/rumen/dereplicated"               # Dereplicated reference MAGs
# RUMENREF_MAGS="/path/to/rumen/dereplicated_genomes"       # Final reference MAG set
# RUMEN_REF_MAGS_DIR="/path/to/novel/rumen/mags"           # Novel rumen MAGs
# RUMEN_ADDED_MAGS_DIR="/path/to/previously/added/mags"    # Previously processed MAGs
# IS_RUMEN_DATA=true  # Set to true for rumen data, false otherwise

# ============================================================================
# MAG PROCESSING DIRECTORIES (Optional - for custom paths)
# ============================================================================
# Uncomment and modify if using non-standard paths

# INPUT_MAGS_DIR="/path/to/custom/input/mags"              # Custom MAG input directory
# MAGS_DIR="/path/to/dereplicated/genomes"                 # MAGs for annotation
# EVALUATION_INPUT_DIR="/path/to/bins/for/evaluation"      # Custom bins for CheckM

# ============================================================================
# PHYLOGENETIC ANALYSIS PARAMETERS
# ============================================================================
# For tree construction and visualization
OUTGROUP_TAXON="p__Firmicutes"  # Taxonomic outgroup for tree rooting
                                  # Examples: p__Firmicutes, p__Proteobacteria, c__Bacilli

# Paths for enhanced tree visualization (modify based on your output)
TAXONOMY_SOURCE="$OUTPUT_BASE/Novel_Mags/gtdbtk/gtdbtk.bac120.summary.tsv"
NOVEL_MAGS_FILE="$OUTPUT_BASE/Novel_Mags/novel_mags_list.txt"
CAZYME_ANNOTATIONS="$OUTPUT_BASE/cazyme_annotations/category_counts.xlsx"
COG_ANNOTATIONS="$OUTPUT_BASE/functional_analysis/files/cog/cog_summary_by_category.xlsx"
KEGG_ANNOTATIONS="$OUTPUT_BASE/functional_analysis/files/kegg/ko_abundance_by_level1.xlsx"
ANNOTATION_TYPE="auto"  # Options: eggnog, kegg, dbcan, auto

# ============================================================================
# ABUNDANCE ESTIMATION PARAMETERS
# ============================================================================
TAXONOMIC_LEVELS="S G"  # Taxonomic levels: S (species), G (genus), F (family), 
                        # O (order), C (class), P (phylum), K (kingdom)
READ_LENGTH=150         # Sequencing read length (e.g., 100, 150, 250)
THRESHOLD=10            # Minimum read count threshold for Bracken

# Optional: Metadata file for group comparisons (CSV format)
# METADATA_FILE="/path/to/metadata.csv"  # Format: Sample,Treatment

# ============================================================================
# DATABASE BUILD OPTIONS
# ============================================================================
MERGE_MODE=true  # true: merge with existing database, false: build fresh database

# ============================================================================
# PIPELINE EXECUTION - DO NOT MODIFY BELOW THIS LINE
# ============================================================================

echo "=========================================================="
echo "MetaMAG Explorer Pipeline Execution"
echo "=========================================================="
echo "Start time: $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID:-LOCAL}"
echo "Running on node: ${SLURM_NODELIST:-$(hostname)}"
echo "Selected step(s): $STEPS"
echo "Output directory: $OUTPUT_BASE"
echo "=========================================================="

# Change to working directory
cd $WORK_DIR || { echo "ERROR: Cannot access $WORK_DIR"; exit 1; }

# Activate Conda environment
source $CONDA_ENV || { echo "ERROR: Cannot activate conda environment"; exit 1; }

# Create log directory if it doesn't exist
mkdir -p $LOG_DIR

# ============================================================================
# STEP DISPATCH LOGIC
# ============================================================================

# Common parameters for all steps
COMMON_PARAMS="--project_config $PROJECT_CONFIG --batch_size $BATCH_SIZE --cpus $CPUS --memory $MEMORY --time $TIME --log_dir $LOG_DIR"

# Execute selected step(s)
case "$STEPS" in
    # =========== Quality Control ===========
    qc|multiqc)
        echo "Running $STEPS..."
        if [ "$STEPS" = "multiqc" ]; then
            python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --input_dir $OUTPUT_BASE/QC
        else
            python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        fi
        ;;
    
    # =========== Read Processing ===========
    trimming)
        echo "Running read trimming with fastp..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        ;;
    
    host_removal)
        echo "Running host contamination removal..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --reference $REFERENCE_GENOME
        ;;
    
    # =========== Assembly ===========
    single_assembly|co_assembly)
        echo "Running $STEPS..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        ;;
    
    metaquast*)
        echo "Running assembly quality assessment..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --assembly_type $ASSEMBLY_TYPE
        ;;
    
    # =========== Binning & Refinement ===========
    *binning*)
        echo "Running binning: $STEPS..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --binning_methods $BINNING_METHODS
        ;;
    
    *bin_refinement*)
        echo "Running bin refinement with DAS Tool..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --score_threshold 0.5
        ;;
    
    # =========== Quality Assessment ===========
    evaluation)
        echo "Running MAG quality evaluation with CheckM..."
        if [ -n "$EVALUATION_INPUT_DIR" ]; then
            python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --evaluation_input_dir $EVALUATION_INPUT_DIR
        else
            python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        fi
        ;;
    
    dRep)
        echo "Running MAG dereplication..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        ;;
    
    # =========== Taxonomy ===========
    gtdbtk)
        echo "Running GTDB-Tk taxonomic classification..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
        ;;
    
    # =========== Novel MAG Processing ===========
    identify_novel_mags)
        echo "Identifying novel MAGs..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
        ;;
    
    process_novel_mags)
        echo "Processing novel MAGs..."
        PROCESS_CMD="python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH"
        
        if [ "$MERGE_MODE" = true ]; then
            PROCESS_CMD="$PROCESS_CMD --merge_mode"
        else
            PROCESS_CMD="$PROCESS_CMD --no_merge"
        fi
        
        if [ "$IS_RUMEN_DATA" = true ]; then
            PROCESS_CMD="$PROCESS_CMD --is_rumen_data --rumen_ref_mags_dir $RUMEN_REF_MAGS_DIR --rumen_added_mags_dir $RUMEN_ADDED_MAGS_DIR"
        fi
        
        eval $PROCESS_CMD
        ;;
    
    # =========== Functional Annotation ===========
    eggnog_annotation|dbcan_annotation)
        echo "Running $STEPS..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --mags_dir ${MAGS_DIR:-$OUTPUT_BASE/Bin_Refinement/drep/dRep_output/dereplicated_genomes}
        ;;
    
    functional_analysis)
        echo "Running integrated functional analysis..."
        if [ -n "$KEGG_DB_FILE" ]; then
            python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --kegg_db_file $KEGG_DB_FILE
        else
            python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
        fi
        ;;
    
    advanced_visualizations)
        echo "Generating advanced visualizations..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
        ;;
    
    # =========== Abundance Estimation ===========
    kraken_abundance)
        echo "Running Kraken2 taxonomic profiling..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH
        ;;
    
    abundance_estimation)
        echo "Running comprehensive abundance estimation with Bracken..."
        ABUNDANCE_CMD="python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS"
        ABUNDANCE_CMD="$ABUNDANCE_CMD --kraken_db_path $KRAKEN_DB_PATH"
        ABUNDANCE_CMD="$ABUNDANCE_CMD --taxonomic_levels $TAXONOMIC_LEVELS"
        ABUNDANCE_CMD="$ABUNDANCE_CMD --read_length $READ_LENGTH --threshold $THRESHOLD"
        
        if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
            ABUNDANCE_CMD="$ABUNDANCE_CMD --metadata_file $METADATA_FILE"
        fi
        
        eval $ABUNDANCE_CMD
        ;;
    
    # =========== Phylogenetic Analysis ===========
    mags_tree)
        echo "Building phylogenetic tree..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --mags_dir ${MAGS_DIR:-$OUTPUT_BASE/Bin_Refinement/drep/dRep_output/dereplicated_genomes} --outgroup_taxon "$OUTGROUP_TAXON"
        ;;
    
    tree_visualization)
        echo "Creating tree visualization..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --taxonomy_source $TAXONOMY_SOURCE
        ;;
    
    enhanced_tree_visualization)
        echo "Creating enhanced tree visualization..."
        TREE_CMD="python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --taxonomy_source $TAXONOMY_SOURCE"
        
        [ -f "$NOVEL_MAGS_FILE" ] && TREE_CMD="$TREE_CMD --novel_mags_file $NOVEL_MAGS_FILE"
        [ -f "$CAZYME_ANNOTATIONS" ] && TREE_CMD="$TREE_CMD --cazyme_annotations $CAZYME_ANNOTATIONS"
        [ -f "$COG_ANNOTATIONS" ] && TREE_CMD="$TREE_CMD --cog_annotations $COG_ANNOTATIONS"
        [ -f "$KEGG_ANNOTATIONS" ] && TREE_CMD="$TREE_CMD --kegg_annotations $KEGG_ANNOTATIONS"
        
        TREE_CMD="$TREE_CMD --annotation_type $ANNOTATION_TYPE"
        eval $TREE_CMD
        ;;
    
    # =========== Rumen-specific Steps ===========
    rumen_refmags_download)
        echo "Downloading rumen reference MAGs..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --dataset_dir $RUMENREF_DATASET
        ;;
    
    rumen_refmags_drep)
        echo "Dereplicating rumen reference MAGs..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --dataset_dir $RUMENREF_DATASET --drep_output_dir $RUMENREF_DREP
        ;;
    
    rumen_drep)
        echo "Integrating project MAGs with rumen references..."
        python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --ref_mags_dir $RUMENREF_MAGS
        ;;
    
    # =========== Default Case ===========
    *)
        echo "Running pipeline step: $STEPS..."
        python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
        ;;
esac

# ============================================================================
# JOB COMPLETION
# ============================================================================
echo "=========================================================="
echo "Pipeline execution completed at $(date)"
echo "Job duration: $SECONDS seconds"
echo "=========================================================="

# Display step-specific output locations
case "$STEPS" in
    qc)
        echo " QC reports available in: $OUTPUT_BASE/QC/"
        ;;
    trimming)
        echo " Trimmed reads available in: $OUTPUT_BASE/Trimming/"
        ;;
    host_removal)
        echo " Host-filtered reads available in: $OUTPUT_BASE/Host_Removal/"
        ;;
    *assembly*)
        echo " Assemblies available in: $OUTPUT_BASE/Assembly/"
        ;;
    *binning*)
        echo " Bins available in: $OUTPUT_BASE/Binning/"
        ;;
    evaluation)
        echo " CheckM results available in: $OUTPUT_BASE/Evaluation/"
        ;;
    dRep)
        echo " Dereplicated MAGs available in: $OUTPUT_BASE/Bin_Refinement/drep/"
        ;;
    gtdbtk)
        echo " Taxonomy results available in: $OUTPUT_BASE/Novel_Mags/gtdbtk/"
        ;;
    *annotation*)
        echo " Functional annotations available in: $OUTPUT_BASE/Functional_Annotation/"
        ;;
    abundance_estimation)
        echo " Abundance results available in:"
        echo "  - Kraken outputs: $OUTPUT_BASE/Kraken_Abundance/"
        echo "  - Bracken estimates: $OUTPUT_BASE/Bracken_Abundance_*/"
        echo "  - Merged matrices: $OUTPUT_BASE/Merged_Bracken_Outputs/"
        echo "  - Visualizations: $OUTPUT_BASE/Abundance_Plots/"
        echo "  - Summary report: $OUTPUT_BASE/abundance_summary.html"
        ;;
    *tree*)
        echo " Phylogenetic trees available in: $OUTPUT_BASE/Phylogeny/"
        ;;
    advanced_visualizations)
        echo " Visualizations available in: $OUTPUT_BASE/Advanced_visualizations/"
        ;;
esac

echo "=========================================================="
echo "Check log files for details: $LOG_DIR"
echo "=========================================================="

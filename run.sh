#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 2:30:00
#SBATCH --output=pipeline_meta_%j.out
#SBATCH --error=pipeline_meta_%j.err

# ===================================================================
# MetaG Pipeline Runner Script
# Executes specific steps of the metagenomics pipeline
# ===================================================================

# ===================================================================
# BASIC CONFIGURATION
# ===================================================================
WORK_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0"
CONDA_ENV="/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate metamag"
SAMPLES="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/Data_Set/Rumen_Data/sra_list.txt"
PROJECT_CONFIG="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/rumen_project_config.yaml"
LOG_DIR="./logs/rumen/advanced_visualizations2"

# ===================================================================
# STEP SELECTION
# ===================================================================
STEPS="advanced_visualizations"  # Select step to run

# ===================================================================
# RESOURCE ALLOCATION
# ===================================================================
BATCH_SIZE=16
CPUS=8
MEMORY="20G"
TIME="15-12:00:00"

# ===================================================================
# ASSEMBLY & QC PARAMETERS
# ===================================================================
ASSEMBLY_TYPE="both"  # Options: idba, megahit, both
BINNING_METHODS="maxbin2 concoct"

# ===================================================================
# REFERENCE GENOME & DATABASE PATHS
# ===================================================================
REFERENCE_GENOME="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/Data_Set/Plant_Data/Reference_Genome/GCF_002870075.4_Lsat_Salinas_v11_genomic.fna"
#KRAKEN_DB_PATH="/usr/home/ironbank/researchdata/metagenome/meta_pipeline/db/Struo2/GTDB_kreken2_207/GTDB_release207/kraken2"
#KEGG_DB_FILE="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/database/ko.txt"

# ===================================================================
# RUMEN REFERENCE DATA PATHS
# ===================================================================
#RUMENREF_DATASET="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/rumenref_mags"
#RUMENREF_DREP="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/rumenref_drep"
#RUMENREF_MAGS="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/rumenref_drep/dereplicated_genomes"
#RUMEN_REF_MAGS_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/novel_rumenref_mags" 
#RUMEN_ADDED_MAGS_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/RumenRefMAGs_Added"
#IS_RUMEN_DATA=true  # Set to true for rumen data, false otherwise

# ===================================================================
# MAG PROCESSING DIRECTORIES
# ===================================================================
#INPUT_MAGS_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/input_mags"
#MAGS_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/Bin_Refinement/drep/dRep_output/dereplicated_genomes"
#EVALUATION_INPUT_DIR="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline/custom_mags_dir"

# ===================================================================
# PHYLOGENETIC ANALYSIS PARAMETERS
# ===================================================================
#OUTGROUP_TAXON="p__Verrucomicrobiota"
TAXONOMY_SOURCE="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/Novel_Mags/gtdbtk/gtdbtk.bac120.summary.tsv"

# ===================================================================
# ENHANCED TREE VISUALIZATION PARAMETERS
# ===================================================================
NOVEL_MAGS_FILE="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/novel_mags_list.txt"
CAZYME_ANNOTATIONS="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/cazyme_annotations/category_counts.xlsx"
COG_ANNOTATIONS="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/functional_analysis/files/cog/cog_summary_by_category.xlsx"
KEGG_ANNOTATIONS="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/pipeline_test/MetaMAG-1.0.0/run/Plant_Analysis/functional_analysis/files/kegg/ko_abundance_by_level1.xlsx"
ANNOTATION_TYPE="auto"  # Options: eggnog, kegg, dbcan, auto


# ===================================================================
# ABUNDANCE ESTIMATION PARAMETERS
# ===================================================================
KRAKEN_DB_PATH="/usr/home/qgg/zexi/INCOME/meta_pipeline/pipeline_database/human"
TAXONOMIC_LEVELS="S G"  # Species and Genus (can add F O C P K)
READ_LENGTH=150         # Adjust based on your sequencing
THRESHOLD=10
# OPTIONAL: METADATA FILE (uncomment if available)
#METADATA_FILE="/usr/home/qgg/zexi/INCOME/Exp2/MetaG_Exp2_HSink_analysis/abundance/abundance_new_db/metadata.csv"

# ===================================================================
# DATABASE OPTIONS
# ===================================================================
#NO_MERGE=false  # Set to true to build fresh database, false to merge with existing

# ===================================================================
# PIPELINE EXECUTION
# ===================================================================

echo "=========================================================="
echo "MetaG Pipeline Execution Started"
echo "=========================================================="
echo "Start time: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Selected step: $STEPS"
echo "=========================================================="

# Change to working directory
cd $WORK_DIR

# Activate Conda environment
source $CONDA_ENV

# ===================================================================
# STEP DISPATCH LOGIC
# ===================================================================

# Common parameters
COMMON_PARAMS="--project_config $PROJECT_CONFIG --batch_size $BATCH_SIZE --cpus $CPUS --memory $MEMORY --time $TIME --log_dir $LOG_DIR"

# Split steps into an array
IFS=' ' read -ra STEP_ARRAY <<< "$STEPS"


case "$STEPS" in
	# ===============================================================
	# ASSEMBLY QUALITY ASSESSMENT
	# ===============================================================
	*"metaquast"*)
		echo "Running MetaQUAST for individual sample assemblies..."
		python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --assembly_type $ASSEMBLY_TYPE
		;;
	
	*"metaquast_coassembly"*)
		echo "Running MetaQUAST for co-assembly..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
		;;
	
	# ===============================================================
	# READ PROCESSING
	# ===============================================================
	*"host_removal"*)
		echo "Running host removal..."
		python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --reference $REFERENCE_GENOME
		;;
	
	# ===============================================================
	# BINNING
	# ===============================================================
	*"single_binning"*|*"co_binning"*)
		echo "Running binning step: $STEPS..."
		python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --binning_methods $BINNING_METHODS
		;;
	
	# ===============================================================
	# BIN EVALUATION
	# ===============================================================
	*"evaluation"*)
		echo "Running bin evaluation..."
		if [ -n "$EVALUATION_INPUT_DIR" ]; then
			python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --evaluation_input_dir $EVALUATION_INPUT_DIR
		else
			python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
		fi
		;;
	
	# ===============================================================
	# RUMEN REFERENCE MAG PROCESSING
	# ===============================================================
	*"rumen_refmags_download"*)
		echo "Downloading rumen reference MAGs..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --dataset_dir $RUMENREF_DATASET
		;;
	
	*"rumen_refmags_drep"*)
		echo "Dereplicating rumen reference MAGs..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --dataset_dir $RUMENREF_DATASET --drep_output_dir $RUMENREF_DREP
		;;
	
	*"rumen_drep"*)
		echo "Combining MAGs with rumen reference..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --ref_mags_dir $RUMENREF_MAGS
		;;
	
	# ===============================================================
	# DATABASE CONSTRUCTION
	# ===============================================================
	*"build_kraken_db"*)
		echo "Building Kraken2 database..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH
		;;
	
	# ===============================================================
	# NOVEL MAG PROCESSING
	# ===============================================================
	*"process_novel_mags"*)
		echo "Processing novel MAGs..."
		
		# Build command with common parameters
		PROCESS_CMD="python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH"
		
		# Add input MAGs directory if specified
		if [ -n "$INPUT_MAGS_DIR" ]; then
			PROCESS_CMD="$PROCESS_CMD --input_mags_dir $INPUT_MAGS_DIR"
		fi
		
		# Add merge mode setting
		if [ "$NO_MERGE" = true ]; then
			PROCESS_CMD="$PROCESS_CMD --no_merge"
			echo "  - Building fresh database (merge mode disabled)"
		else
			echo "  - Merging with existing database"
		fi
		
		# Handle rumen-specific processing
		if [ "$IS_RUMEN_DATA" = true ]; then
			if [ -z "$RUMEN_REF_MAGS_DIR" ] || [ -z "$RUMEN_ADDED_MAGS_DIR" ]; then
				echo "[ERROR] For rumen data processing, both RUMEN_REF_MAGS_DIR and RUMEN_ADDED_MAGS_DIR must be set"
				exit 1
			fi
			PROCESS_CMD="$PROCESS_CMD --is_rumen_data --rumen_ref_mags_dir $RUMEN_REF_MAGS_DIR --rumen_added_mags_dir $RUMEN_ADDED_MAGS_DIR"
			echo "  - Processing as rumen data"
		fi
		
		# Execute the command
		eval $PROCESS_CMD
		;;
	
	# ===============================================================
	# FUNCTIONAL ANNOTATION
	# ===============================================================
	*"eggnog_annotation"*|*"dbcan_annotation"*)
		echo "Running $STEPS..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --mags_dir $MAGS_DIR
		;;
	
	*"functional_analysis"*)
		echo "Running functional analysis..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --kegg_db_file $KEGG_DB_FILE
		;;

        *"advanced_visualizations"*)
		echo "Running functional analysis..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS
		;;

	
	# ===============================================================
	# ABUNDANCE ESTIMATION
	# ===============================================================
	*"kraken_abundance"*)
		echo "Running Kraken abundance estimation..."
		python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH
		;;
	
	*"abundance_estimation"*)
		echo "Running comprehensive abundance estimation..."
		echo "  - Taxonomic levels: $TAXONOMIC_LEVELS"
		echo "  - Read length: $READ_LENGTH"
		echo "  - Threshold: $THRESHOLD"
		
		# Build abundance command
		ABUNDANCE_CMD="python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS --kraken_db_path $KRAKEN_DB_PATH --taxonomic_levels $TAXONOMIC_LEVELS --read_length $READ_LENGTH --threshold $THRESHOLD"
		
		# Add metadata file if available
		if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
			echo "  - Using metadata file: $METADATA_FILE"
			ABUNDANCE_CMD="$ABUNDANCE_CMD --metadata_file $METADATA_FILE"
		else
			echo "  - No metadata file provided - will create automatic groups"
		fi
		
		# Execute the command
		eval $ABUNDANCE_CMD
		;;
	
	# ===============================================================
	# PHYLOGENETIC ANALYSIS
	# ===============================================================
	*"mags_tree"*)
		echo "Constructing MAGs phylogenetic tree..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --mags_dir $MAGS_DIR --outgroup_taxon "$OUTGROUP_TAXON"
		;;
	
	*"tree_visualization"*)
		echo "Visualizing phylogenetic tree..."
		python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --taxonomy_source $TAXONOMY_SOURCE
		;;

        *"enhanced_tree_visualization"*)
            echo "Running enhanced tree visualization with novel MAG highlighting and functional annotations..."
    
            # Build enhanced visualization command
            python3 -m MetaMAG.main $COMMON_PARAMS --steps $STEPS --taxonomy_source $TAXONOMY_SOURCE --novel_mags_file $NOVEL_MAGS_FILE --cazyme_annotations $CAZYME_ANNOTATIONS --cog_annotations $COG_ANNOTATIONS --kegg_annotations $KEGG_ANNOTATIONS --annotation_type $ANNOTATION_TYPE
            ;;
	
	# ===============================================================
	# DEFAULT CASE
	# ===============================================================
	*)
		echo "Running default pipeline step: $STEPS..."
		python3 -m MetaMAG.main $COMMON_PARAMS --samples-file $SAMPLES --steps $STEPS
		;;
        
esac

# ===================================================================
# COMPLETION MESSAGE
# ===================================================================
echo "=========================================================="
echo "Pipeline completed at $(date)"
echo "=========================================================="

# Show results location for abundance estimation
if [ "$STEPS" == "abundance_estimation" ]; then
    echo ""
    echo "Abundance estimation results can be found in:"
    echo "  - Kraken_Abundance/: Raw Kraken2 outputs"
    echo "  - Bracken_Abundance_S/: Species-level Bracken estimates"
    echo "  - Bracken_Abundance_G/: Genus-level Bracken estimates"
    echo "  - Merged_Bracken_Outputs/: Combined abundance matrices"
    echo "  - Abundance_Plots/: Publication-quality visualizations"
    echo "  - abundance_summary.html: Interactive summary dashboard"
    echo ""
    if [ -n "$METADATA_FILE" ] && [ -f "$METADATA_FILE" ]; then
        echo "  ?? Plots include group comparisons based on metadata"
    else
        echo "  ?? Plots show automatic grouping based on sample names"
        echo "    For better group comparisons, provide a metadata file"
    fi
    echo "=========================================================="
fi
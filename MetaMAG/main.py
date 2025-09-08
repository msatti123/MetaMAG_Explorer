#!/usr/bin/env python3
"""
Pipeline job submission script for MetaG pipeline.
Handles SLURM job creation and submission with proper argument configuration.
"""

import argparse
import subprocess
import os
from MetaMAG.config import config


# ===================================================================
# PIPELINE CONFIGURATION
# ===================================================================

# Define steps that don't require sample processing
NO_SAMPLES_REQUIRED = {
    "multiqc","metaquast_coassembly", "rumen_refmags_download", "rumen_refmags_drep", 
    "co_bin_refinement", "rumen_refmags_gtdbtk", "gtdbtk", "rumen_drep", 
    "add_mags_to_repo", "build_kraken_db", "identify_novel_mags", 
    "process_novel_mags", "eggnog_annotation", "dbcan_annotation", 
    "functional_analysis", "mags_tree", "tree_visualization",
    # ADD NEW STEP-BY-STEP NOVEL MAG STEPS:
    "repository_dereplication", "rumen_dereplication", "add_to_repository",
    "prepare_taxonomy_input", "create_taxonomy", "rename_headers", 
    "add_kraken_library", "build_kraken_database", "build_bracken_database", 
    "build_database_distribution", "enhanced_tree_visualization", "advanced_visualizations",

}

# Define all valid pipeline steps
VALID_STEPS = {
    # Quality Control & Preprocessing
    "qc", "multiqc", "trimming", "host_removal",
    
    # Assembly & Assessment  
    "single_assembly", "co_assembly", "metaquast", "metaquast_coassembly",
    
    # Binning & Refinement
    "single_binning", "co_binning", "single_bin_refinement", "co_bin_refinement",
    "dRep", "evaluation",
    
    # Taxonomic Classification
    "gtdbtk", "novel_genome_refinement",
    
    # Rumen Reference MAGs
    "rumen_refmags_download", "rumen_refmags_drep", "rumen_refmags_gtdbtk", "rumen_drep",
    
    # Novel MAG Processing
    "identify_novel_mags", "process_novel_mags", "add_mags_to_repo",

    # Novel MAG Processing - Step-by-Step
    "repository_dereplication", "rumen_dereplication", "add_to_repository",
    "prepare_taxonomy_input", "create_taxonomy", "rename_headers",
    "add_kraken_library", "build_kraken_database", "build_bracken_database",
    "build_database_distribution",
    
    # Database Construction
    "build_kraken_db",
    
    # Functional Annotation
    "eggnog_annotation", "dbcan_annotation", "functional_analysis", "advanced_visualizations",
    
    # Abundance Analysis
    "kraken_abundance", "abundance_estimation",
    
    # Phylogenetic Analysis
    "mags_tree", "tree_visualization", "enhanced_tree_visualization", 
}


# ===================================================================
# UTILITY FUNCTIONS
# ===================================================================

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    os.makedirs(directory, exist_ok=True)
    print(f"[DEBUG] Ensured directory exists: {directory}")


def submit_job(script_path):
    """Submit a job to SLURM using sbatch."""
    print(f"Submitting job: {script_path}")
    try:
        result = subprocess.run(
            f"sbatch {script_path}", 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True
        )
        print(f"[DEBUG] SLURM Response: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to submit job: {e}")
        print(f"[ERROR] stderr: {e.stderr}")


def read_samples(samples, samples_file, steps):
    """
    Read samples from either a provided list or a text file.
    Returns empty list for steps that don't require samples.
    """
    # If we only have steps that don't need samples, return empty list
    if all(step in NO_SAMPLES_REQUIRED for step in steps):
        return []
        
    if samples:
        return samples  # Direct list input
    
    if samples_file:
        if os.path.exists(samples_file):
            with open(samples_file, "r") as f:
                return [line.strip() for line in f if line.strip()]
        else:
            print(f"[ERROR] Samples file not found: {samples_file}")
            exit(1)

    # Only show error if we have steps that require samples
    if any(step not in NO_SAMPLES_REQUIRED for step in steps):
        print("[ERROR] You must provide samples using either --samples or --samples-file.")
        exit(1)
    
    return []


def batch_samples(samples, batch_size):
    """Split samples into batches of specified size."""
    # If no samples (e.g., for rumen_refmags steps), create a single empty batch
    if not samples:
        return [[]]
    
    return [samples[i:i + batch_size] for i in range(0, len(samples), batch_size)]


# ===================================================================
# SLURM SCRIPT GENERATION
# ===================================================================

def get_conda_environment():
    """Get Conda environment path from config."""
    conda_activate = config["environment"].get("conda_env", "")
    if not conda_activate:
        print("[ERROR] Conda environment path not found in config.py. Exiting.")
        exit(1)
    return conda_activate


def create_slurm_header(cpus, memory, time, output_dir):
    """Create SLURM header with resource specifications."""
    return f"""#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {cpus}
#SBATCH --mem={memory}
#SBATCH -t {time}
#SBATCH --output={output_dir}/slurm_%j.out
#SBATCH --error={output_dir}/slurm_%j.err

# Activate Conda environment
source {get_conda_environment()}

"""


def build_step_command(step, samples, project_config, cpus, memory, time, **kwargs):
    """Build command for a specific pipeline step with all required arguments."""
    # Basic required arguments
    cmd_args = f"--project_config {project_config} --step {step} --cpus {cpus} --memory {memory} --time {time}"
    
    # Add samples if needed for this step
    if step not in NO_SAMPLES_REQUIRED:
        cmd_args += f" --samples {' '.join(samples)}"
    
    # Step-specific argument handling
    step_args = build_step_specific_args(step, **kwargs)
    if step_args:
        cmd_args += " " + step_args
    
    return f"""
echo "Running step: {step}"
python3 run_step.py {cmd_args}
"""


def build_step_specific_args(step, **kwargs):
    """Build step-specific arguments based on the pipeline step."""
    args = []

    # MultiQC
    if step == "multiqc" and kwargs.get('input_dir'):
        args.append(f"--input_dir {kwargs['input_dir']}") 
 
    # Assembly & Quality Control Arguments
    if step == "metaquast" and kwargs.get('assembly_type'):
        args.append(f"--assembly_type {kwargs['assembly_type']}")
    
    # Binning Arguments
    if step in ["single_binning", "co_binning"] and kwargs.get('binning_methods'):
        args.append(f"--binning_methods {' '.join(kwargs['binning_methods'])}")
    
    # Host Removal Arguments
    if step == "host_removal" and kwargs.get('reference'):
        args.append(f"--reference {kwargs['reference']}")
    
    # DAS Tool Arguments
    if step in ["single_bin_refinement", "co_bin_refinement"]:
        args.append(f"--score_threshold {kwargs.get('score_threshold', 0.5)}")
    
    # Evaluation Arguments
    if step == "evaluation" and kwargs.get('evaluation_input_dir'):
        args.append(f"--evaluation_input_dir {kwargs['evaluation_input_dir']}")
    
    # Rumen Reference MAGs Arguments
    if step in ["rumen_refmags_download", "rumen_refmags_drep"] and kwargs.get('dataset_dir'):
        args.append(f"--dataset_dir {kwargs['dataset_dir']}")
    
    if step == "rumen_refmags_drep" and kwargs.get('drep_output_dir'):
        args.append(f"--drep_output_dir {kwargs['drep_output_dir']}")
    
    if step == "rumen_drep" and kwargs.get('ref_mags_dir'):
        args.append(f"--ref_mags_dir {kwargs['ref_mags_dir']}")
    
    # GTDB-Tk Arguments
    if step in ["gtdbtk", "rumen_refmags_gtdbtk"]:
        if kwargs.get('genome_dir'):
            args.append(f"--genome_dir {kwargs['genome_dir']}")
        if kwargs.get('gtdbtk_out_dir'):
            args.append(f"--gtdbtk_out_dir {kwargs['gtdbtk_out_dir']}")
        args.append(f"--genome_extension {kwargs.get('genome_extension', '.fa')}")
        if kwargs.get('is_rumen_data') and step == "gtdbtk":
            args.append("--is_rumen_data")
    
    if step == "rumen_refmags_gtdbtk":
        if kwargs.get('novel_output_dir'):
            args.append(f"--novel_output_dir {kwargs['novel_output_dir']}")
        if kwargs.get('only_novel'):
            args.append("--only_novel")
    
    # Database Arguments
    if step in ["build_kraken_db", "kraken_abundance", "abundance_estimation"] and kwargs.get('kraken_db_path'):
        args.append(f"--kraken_db_path {kwargs['kraken_db_path']}")
    
    # Novel MAGs Processing Arguments
    if step == "process_novel_mags":
        if kwargs.get('is_rumen_data'):
            args.append("--is_rumen_data")
            if kwargs.get('rumen_ref_mags_dir'):
                args.append(f"--rumen_ref_mags_dir {kwargs['rumen_ref_mags_dir']}")
            if kwargs.get('rumen_added_mags_dir'):
                args.append(f"--rumen_added_mags_dir {kwargs['rumen_added_mags_dir']}")
        if kwargs.get('input_mags_dir'):
            args.append(f"--input_mags_dir {kwargs['input_mags_dir']}")
        if kwargs.get('kraken_db_path'):
            args.append(f"--kraken_db_path {kwargs['kraken_db_path']}")
        if kwargs.get('custom_gtdbtk_dir'):
            args.append(f"--custom_gtdbtk_dir {kwargs['custom_gtdbtk_dir']}")
        if not kwargs.get('no_merge'):
            args.append("--merge_mode")
    
    if step == "identify_novel_mags" and kwargs.get('custom_gtdbtk_dir'):
        args.append(f"--custom_gtdbtk_dir {kwargs['custom_gtdbtk_dir']}")
    
    # Step-by-Step Novel MAG Processing Arguments
    if step == "rumen_dereplication":
        if kwargs.get('rumen_ref_mags_dir'):
            args.append(f"--rumen_ref_mags_dir {kwargs['rumen_ref_mags_dir']}")
        if kwargs.get('rumen_added_mags_dir'):
            args.append(f"--rumen_added_mags_dir {kwargs['rumen_added_mags_dir']}")
    
    if step in ["prepare_taxonomy_input", "create_taxonomy", "rename_headers", 
                "add_kraken_library", "build_kraken_database", "build_bracken_database", 
                "build_database_distribution"]:
        if kwargs.get('kraken_db_path'):
            args.append(f"--kraken_db_path {kwargs['kraken_db_path']}")
    
    if step in ["prepare_taxonomy_input", "create_taxonomy"]:
        if kwargs.get('custom_gtdbtk_dir'):
            args.append(f"--custom_gtdbtk_dir {kwargs['custom_gtdbtk_dir']}")
        if kwargs.get('rumen_ref_mags_dir'):
            args.append(f"--rumen_ref_mags_dir {kwargs['rumen_ref_mags_dir']}")
    
    if step == "create_taxonomy":
        if not kwargs.get('no_merge'):
            args.append("--merge_mode")
    
    if step in ["build_bracken_database", "build_database_distribution"]:
        args.append(f"--read_length {kwargs.get('read_length', 150)}")
    
    # Functional Annotation Arguments
    if step in ["eggnog_annotation", "dbcan_annotation"] and kwargs.get('mags_dir'):
        args.append(f"--mags_dir {kwargs['mags_dir']}")
    
    if step == "functional_analysis" and kwargs.get('kegg_db_file'):
        args.append(f"--kegg_db_file {kwargs['kegg_db_file']}")
    
    if step == "advanced_visualizations":
        if kwargs.get('mag_lists'):
            args.append(f"--mag_lists {' '.join(kwargs['mag_lists'])}")
        if kwargs.get('list_names'):
            args.append(f"--list_names {' '.join(kwargs['list_names'])}")
        if kwargs.get('viz_output_suffix'):
            args.append(f"--viz_output_suffix {kwargs['viz_output_suffix']}")
        if kwargs.get('viz_types'):
            args.append(f"--viz_types {' '.join(kwargs['viz_types'])}")
    
    # Phylogenetic Analysis Arguments
    if step == "mags_tree":
        if kwargs.get('mags_dir'):
            args.append(f"--mags_dir {kwargs['mags_dir']}")
        if kwargs.get('outgroup_taxon'):
            args.append(f"--outgroup_taxon {kwargs['outgroup_taxon']}")
        else:
            print("[WARNING] No outgroup taxon provided for mags_tree step. This is required.")
    
    if step == "tree_visualization":
        if kwargs.get('taxonomy_source'):
            args.append(f"--taxonomy_source {kwargs['taxonomy_source']}")
        if kwargs.get('mags_dir'):
            args.append(f"--mags_dir {kwargs['mags_dir']}")
    
    if step == "enhanced_tree_visualization":
        if kwargs.get('taxonomy_source'):
            args.append(f"--taxonomy_source {kwargs['taxonomy_source']}")
        if kwargs.get('novel_mags_file'):
            args.append(f"--novel_mags_file {kwargs['novel_mags_file']}")
        if kwargs.get('functional_annotations'):
            args.append(f"--functional_annotations {kwargs['functional_annotations']}")
        if kwargs.get('annotation_type'):
            args.append(f"--annotation_type {kwargs['annotation_type']}")
        if kwargs.get('cazyme_annotations'):
            args.append(f"--cazyme_annotations {kwargs['cazyme_annotations']}")
        if kwargs.get('cog_annotations'):
            args.append(f"--cog_annotations {kwargs['cog_annotations']}")
        if kwargs.get('kegg_annotations'):
            args.append(f"--kegg_annotations {kwargs['kegg_annotations']}")

    # Abundance Estimation Arguments
    if step == "abundance_estimation":
        if kwargs.get('taxonomic_levels'):
            args.append(f"--taxonomic_levels {' '.join(kwargs['taxonomic_levels'])}")
        args.append(f"--read_length {kwargs.get('read_length', 150)}")
        args.append(f"--threshold {kwargs.get('threshold', 10)}")
        if kwargs.get('metadata_file'):
            args.append(f"--metadata_file {kwargs['metadata_file']}")
    
    return " ".join(args)


def create_slurm_script(samples, steps, cpus, memory, time, output_dir, batch_id, 
                       project_config, **kwargs):
    """Create a SLURM script to process multiple samples in parallel."""
    ensure_directory_exists(output_dir)
    script_name = f"slurm_batch_{batch_id}.sh"
    script_path = os.path.join(output_dir, script_name)

    print(f"[DEBUG] Creating SLURM script: {script_path}")

    # Start with SLURM header
    job_script = create_slurm_header(cpus, memory, time, output_dir)

    # Add commands for each step
    for step in steps:
        job_script += build_step_command(
            step, samples, project_config, cpus, memory, time, **kwargs
        )

    # Write script to file
    try:
        with open(script_path, 'w') as script_file:
            script_file.write(job_script)
        print(f"[DEBUG] SLURM script successfully created: {script_path}")
    except Exception as e:
        print(f"[ERROR] Failed to create SLURM script: {e}")

    return script_path


# ===================================================================
# ARGUMENT PARSING
# ===================================================================

def parse_arguments():
    """Parse command line arguments for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Pipeline to process metagenomic samples in parallel.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available pipeline steps:
  Quality Control & Preprocessing: qc, multiqc, trimming, host_removal
  Assembly & Assessment: single_assembly, co_assembly, metaquast, metaquast_coassembly
  Binning & Refinement: single_binning, co_binning, single_bin_refinement, co_bin_refinement, dRep, evaluation
  Taxonomic Classification: gtdbtk, novel_genome_refinement
  Rumen Reference MAGs: rumen_refmags_download, rumen_refmags_drep, rumen_refmags_gtdbtk, rumen_drep
  
  Novel MAG Processing (Complete): identify_novel_mags, process_novel_mags, add_mags_to_repo
  
  Novel MAG Processing (Step-by-Step):
    repository_dereplication - Dereplicate against existing MAGs repository
    rumen_dereplication - Dereplicate against rumen reference MAGs (if applicable)
    add_to_repository - Add validated novel MAGs to repository
    prepare_taxonomy_input - Prepare taxonomy input files
    create_taxonomy - Create taxonomy files (nodes.dmp, names.dmp)
    rename_headers - Rename FASTA headers for Kraken compatibility
    add_kraken_library - Add MAGs to Kraken library
    build_kraken_database - Build main Kraken database
    build_bracken_database - Build Bracken database for abundance estimation
    build_database_distribution - Build optimized database distribution
  
  Database Construction: build_kraken_db
  Functional Annotation: eggnog_annotation, dbcan_annotation, functional_analysis
  Abundance Analysis: kraken_abundance, abundance_estimation
  Phylogenetic Analysis: mags_tree, tree_visualization
        """
    )
    
    # ===================================================================
    # BASIC ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--project_config", 
        type=str, 
        required=True,
        help="Path to the project configuration file (YAML)."
    )
    parser.add_argument(
        "--samples", 
        nargs="+", 
        default=None,
        help="List of sample IDs (e.g., --samples ERR2200487 ERR2200488)"
    )
    parser.add_argument(
        "--samples-file", 
        type=str,
        help="Path to a text file containing sample IDs (one per line)"
    )
    parser.add_argument(
        "--steps", 
        nargs="*", 
        type=str, 
        required=True,
        help="Steps to run."
    )
    parser.add_argument(
        "--batch_size", 
        type=int, 
        required=True,
        help="Number of samples per job."
    )
    parser.add_argument(
        "--cpus", 
        type=int, 
        required=True,
        help="Number of CPUs per job."
    )
    parser.add_argument(
        "--memory", 
        required=True,
        help="Memory allocation per job."
    )
    parser.add_argument(
        "--time", 
        required=True,
        help="Time limit per job."
    )
    parser.add_argument(
        "--log_dir", 
        default="./logs",
        help="Directory for SLURM logs."
    )
    
    # ===================================================================
    # ASSEMBLY & QC ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--input_dir", 
        type=str,
        help="Input directory to search for QC reports (for MultiQC)"
    )
    parser.add_argument(
        "--assembly_type", 
        type=str, 
        default="both",
        choices=["idba", "megahit", "both"],
        help="Assembly type to evaluate with metaQUAST: idba, megahit, or both"
    )
    parser.add_argument(
        "--reference", 
        type=str,
        help="Reference genome for host removal (step 3: host_removal)",
        default=None
    )
    
    # ===================================================================
    # BINNING ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--binning_methods", 
        nargs="+", 
        default=[],
        help="Binning methods to use (optional)"
    )
    parser.add_argument(
        "--score_threshold", 
        type=float, 
        default=0.5,
        help="Score threshold for DAS Tool (default: 0.5)"
    )
    
    # ===================================================================
    # EVALUATION ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--evaluation_input_dir", 
        type=str,
        help="Custom input directory containing genomes for CheckM evaluation"
    )
    
    # ===================================================================
    # RUMEN REFERENCE MAG ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--dataset_dir", 
        type=str,
        help="Directory for downloading or accessing reference MAGs data."
    )
    parser.add_argument(
        "--drep_output_dir", 
        type=str,
        help="Directory for dRep output on reference MAGs."
    )
    parser.add_argument(
        "--ref_mags_dir", 
        type=str,
        help="Directory containing reference MAGs for combining with project MAGs."
    )
    
    # ===================================================================
    # GTDB-TK ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--genome_dir", 
        type=str,
        help="Directory containing input genomes for GTDB-Tk"
    )
    parser.add_argument(
        "--gtdbtk_out_dir", 
        type=str,
        help="Output directory for GTDB-Tk results"
    )
    parser.add_argument(
        "--genome_extension", 
        type=str, 
        default=".fa",
        help="File extension of input genomes (default: .fa)"
    )
    parser.add_argument(
        "--is_rumen_data", 
        action="store_true",
        help="Whether the data is rumen-related for GTDB-Tk."
    )
    
    # ===================================================================
    # NOVEL MAG PROCESSING ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--input_mags_dir", 
        type=str,
        help="Directory containing input MAGs for processing"
    )
    parser.add_argument(
        "--custom_gtdbtk_dir", 
        type=str,
        help="Custom path to GTDB-Tk results directory"
    )
    parser.add_argument(
        "--novel_output_dir", 
        type=str,
        help="Output directory for novel MAGs (for rumen_refmags_gtdbtk step)"
    )
    parser.add_argument(
        "--only_novel", 
        action="store_true",
        help="Only identify novel MAGs, skip GTDB-Tk if already run"
    )
    parser.add_argument(
        "--rumen_ref_mags_dir", 
        type=str,
        help="Path to rumen reference MAGs directory (required if is_rumen=True)"
    )
    parser.add_argument(
        "--rumen_added_mags_dir", 
        type=str,
        help="Path to previously added rumen MAGs directory (required if is_rumen=True)"
    )
    parser.add_argument(
        "--no_merge", 
        action="store_true",
        help="Disable merge mode (build fresh database)"
    )
    
    # ===================================================================
    # DATABASE ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--kraken_db_path", 
        type=str,
        help="Path to build Kraken database"
    )
    parser.add_argument(
        "--kegg_db_file", 
        type=str,
        help="Path to KEGG database file (ko.txt) for functional analysis"
    )
    
    # ===================================================================
    # ANNOTATION ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--mags_dir", 
        type=str,
        help="Directory containing MAGs for EggNOG annotation"
    )
    parser.add_argument(
        "--mag_lists", 
        nargs="+",
        help="Paths to MAG list files for visualization"
    )
    parser.add_argument(
        "--list_names", 
        nargs="+",
        help="Names for each MAG list"
    )
    parser.add_argument(
        "--viz_output_suffix", 
        type=str,
        default="",
        help="Suffix for visualization output directory"
    )
    parser.add_argument(
        "--viz_types", 
        nargs="+",
        choices=['functional', 'taxonomic', 'integrated', 'abundance', 'comparative'],
        help="Specific visualization types to create"
    )
    
    # ===================================================================
    # PHYLOGENETIC ANALYSIS ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--outgroup_taxon", 
        type=str,
        help="Outgroup taxon for phylogenetic tree construction (e.g., p__Firmicutes)"
    )
    parser.add_argument(
        "--taxonomy_source", 
        type=str,
        help="Path to GTDB-Tk taxonomy file for tree visualization"
    )
    parser.add_argument(
        "--novel_mags_file", 
        type=str,
        help="Path to file containing list of novel MAG names for highlighting"
    )
    parser.add_argument(
        "--functional_annotations", 
        type=str,
        help="Path to functional annotation results for ring visualization"
    )
    parser.add_argument(
        "--annotation_type", 
        type=str,
        choices=["eggnog", "kegg", "dbcan", "auto"],
        default="auto",
        help="Type of functional annotation (default: auto)"
    )
    parser.add_argument(
        "--cazyme_annotations", 
        type=str,
        help="Path to CAZyme annotation Excel file for enhanced tree visualization"
    )
    parser.add_argument(
        "--cog_annotations", 
        type=str,
        help="Path to COG annotation Excel file for enhanced tree visualization"
    )
    parser.add_argument(
        "--kegg_annotations", 
        type=str,
        help="Path to KEGG annotation Excel file for enhanced tree visualization"
    )
    
    # ===================================================================
    # ABUNDANCE ESTIMATION ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--taxonomic_levels", 
        nargs="+", 
        default=['S', 'G'],
        choices=['S', 'G', 'F', 'O', 'C', 'P', 'K'],
        help="Taxonomic levels for abundance estimation (default: S G)"
    )
    parser.add_argument(
        "--read_length", 
        type=int, 
        default=150,
        help="Read length for Bracken (default: 150)"
    )
    parser.add_argument(
        "--threshold", 
        type=int, 
        default=10,
        help="Minimum reads for Bracken classification (default: 10)"
    )
    parser.add_argument(
        "--metadata_file", 
        type=str,
        help="Path to metadata file for abundance visualization"
    )
    
    return parser.parse_args()


# ===================================================================
# MAIN FUNCTION
# ===================================================================

def main():
    """Main function to orchestrate job submission."""
    # Parse arguments
    args = parse_arguments()
    
    # Validate steps
    invalid_steps = [step for step in args.steps if step not in VALID_STEPS]
    if invalid_steps:
        print(f"[ERROR] Invalid steps provided: {', '.join(invalid_steps)}")
        return
    
    # Read samples
    samples = read_samples(args.samples, args.samples_file, args.steps)
    
    # Create batches
    batches = batch_samples(samples, args.batch_size)
    batch_id = 0
    
    # Process each batch
    for batch in batches:
        batch_id += 1
        if batch:
            print(f"Preparing job for batch {batch_id}: {batch}")
        else:
            print(f"Preparing job for batch {batch_id} (no samples required)")
        
        #Create SLURM script
        slurm_script = create_slurm_script(
            batch,  # samples (positional)
            args.steps,  # steps (positional)
            args.cpus,  # cpus (positional)
            args.memory,  # memory (positional)
            args.time,  # time (positional)
            args.log_dir,  # output_dir (positional)
            batch_id,  # batch_id (positional)
            args.project_config,  # project_config (positional)
            # Pass all arguments as kwargs
            input_dir=args.input_dir,
            assembly_type=args.assembly_type,
            binning_methods=args.binning_methods,
            score_threshold=args.score_threshold,
            evaluation_input_dir=args.evaluation_input_dir,
            reference=args.reference,
            dataset_dir=args.dataset_dir,
            drep_output_dir=args.drep_output_dir,
            genome_dir=args.genome_dir,
            gtdbtk_out_dir=args.gtdbtk_out_dir,
            genome_extension=args.genome_extension,
            custom_gtdbtk_dir=args.custom_gtdbtk_dir,
            novel_output_dir=args.novel_output_dir,
            only_novel=args.only_novel,
            ref_mags_dir=args.ref_mags_dir,
            is_rumen_data=args.is_rumen_data,
            kraken_db_path=args.kraken_db_path,
            no_merge=args.no_merge,
            input_mags_dir=args.input_mags_dir,
            rumen_ref_mags_dir=args.rumen_ref_mags_dir,
            rumen_added_mags_dir=args.rumen_added_mags_dir,
            mags_dir=args.mags_dir,
            kegg_db_file=args.kegg_db_file,
            outgroup_taxon=args.outgroup_taxon,
            taxonomy_source=args.taxonomy_source,
            taxonomic_levels=args.taxonomic_levels,
            read_length=args.read_length,
            threshold=args.threshold,
            metadata_file=args.metadata_file,
            novel_mags_file=args.novel_mags_file,
            functional_annotations=args.functional_annotations,
            annotation_type=args.annotation_type,
            cazyme_annotations=args.cazyme_annotations,
            cog_annotations=args.cog_annotations,
            kegg_annotations=args.kegg_annotations,
            mag_lists=args.mag_lists,
            list_names=args.list_names,
            viz_output_suffix=args.viz_output_suffix,
            viz_types=args.viz_types
        )
        
        # Submit job
        submit_job(slurm_script)
    
    print("All jobs submitted successfully.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Script to run individual steps of the MetaG pipeline.
This script is called by SLURM batch job scripts.
"""

import argparse
import yaml
from MetaMAG import (
    qc, multiqc, trimming, host_removal, assembly, quast, mapping, binning,
    das_tool, dereplicate, evaluation, gtdbtk,
    rumen_refmags_download, rumen_refmags_drep, rumen_refmags_gtdbtk, 
    rumen_drep, mags_repository, kraken_database, novel_mags, nmags_pipeline, 
    eggnog, functional_analysis, mags_tree, phylo_tree_vis, enhanced_phylo_tree_vis, kraken_abundance, 
    dbcan, abundance_estimation, advanced_visualizations
)


def parse_arguments():
    """Parse command line arguments for the pipeline runner."""
    parser = argparse.ArgumentParser(
        description="Run a specific pipeline step.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available pipeline steps:
  Quality Control & Preprocessing:
    - qc, multiqc, trimming, host_removal
  
  Assembly & Assessment:
    - single_assembly, co_assembly, metaquast, metaquast_coassembly
  
  Binning & Refinement:
    - single_binning, co_binning, single_bin_refinement, co_bin_refinement
    - dRep, evaluation
  
  Taxonomic Classification:
    - gtdbtk
  
  Rumen Reference MAGs:
    - rumen_refmags_download, rumen_refmags_drep, rumen_refmags_gtdbtk
    - rumen_drep
  
  Novel MAG Processing:
    - identify_novel_mags, process_novel_mags, add_mags_to_repo
  
  Database Construction:
    - build_kraken_db
  
  Functional Annotation:
    - eggnog_annotation, dbcan_annotation, functional_analysis
  
  Abundance Analysis:
    - kraken_abundance, abundance_estimation
  
  Phylogenetic Analysis:
    - mags_tree, tree_visualization
        """
    )
    
    # ===================================================================
    # BASIC ARGUMENTS (REQUIRED FOR ALL STEPS)
    # ===================================================================
    parser.add_argument(
        "--project_config", 
        type=str, 
        required=True,
        help="Path to the project configuration file (YAML)."
    )
    parser.add_argument(
        "--step", 
        type=str, 
        required=True,
        help="Step to execute."
    )
    parser.add_argument(
        "--cpus", 
        type=int, 
        required=True,
        help="Number of CPUs to allocate."
    )
    parser.add_argument(
        "--memory", 
        required=True,
        help="Memory allocation (e.g., 20G)."
    )
    parser.add_argument(
        "--time", 
        required=True,
        help="Time limit (e.g., 5-12:00:00)."
    )

    # ===================================================================
    # MULTIQC ARGUMENTS
    # ===================================================================

    parser.add_argument(
        "--input_dir", 
        type=str,
        help="Input directory to search for QC reports (for MultiQC)"
    )

    # ===================================================================
    # SAMPLE AND REFERENCE ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--samples", 
        nargs="+",
        help="Samples to process (required for most steps)."
    )
    parser.add_argument(
        "--reference", 
        required=False, 
        default=None,
        help="Path to reference genome index."
    )
    
    # ===================================================================
    # ASSEMBLY ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--assembly_type", 
        type=str, 
        default="both",
        choices=["idba", "megahit", "both"],
        help="Assembly type to evaluate with metaQUAST: idba, megahit, or both"
    )
    
    # ===================================================================
    # BINNING ARGUMENTS
    # ===================================================================
    parser.add_argument(
        "--binning_methods", 
        nargs="+", 
        default=[],
        help="Binning methods to use (optional)."
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
        help="Custom input directory containing genomes for CheckM analysis"
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
        help="Directory containing input genomes for GTDB-Tk."
    )
    parser.add_argument(
        "--gtdbtk_out_dir", 
        type=str,
        help="Output directory for GTDB-Tk results."
    )
    parser.add_argument(
        "--genome_extension", 
        type=str, 
        default=".fa",
        help="File extension of input genomes (default: .fa)."
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
        "--merge_mode", 
        action="store_true",
        help="Whether to merge with existing Kraken database"
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
        help="Path to MAGs directory for EggNOG annotation/phylogenetic tree construction"
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
        help="Path to file containing list of novel MAG names (one per line) for highlighting"
    )
    parser.add_argument(
        "--functional_annotations", 
        type=str,
        help="Path to functional annotation summary file or directory"
    )
    parser.add_argument(
        "--annotation_type", 
        type=str,
        choices=["eggnog", "kegg", "dbcan", "auto"],
        default="auto",
        help="Type of functional annotation (default: auto-detect)"
    )
    parser.add_argument(
        "--cazyme_annotations", 
        type=str,
        help="Path to CAZyme annotation Excel file"
    )
    parser.add_argument(
        "--cog_annotations", 
        type=str,
        help="Path to COG annotation Excel file"
    )
    parser.add_argument(
        "--kegg_annotations", 
        type=str,
        help="Path to KEGG annotation Excel file"
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


def load_project_config(config_path):
    """Load project configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def validate_step_requirements(args):
    """Validate that required arguments are provided for the selected step."""
    # Steps that don't require samples
    no_samples_steps = {
        "multiqc", "metaquast_coassembly", "gtdbtk", "rumen_refmags_download", 
        "rumen_refmags_drep", "co_bin_refinement", "rumen_refmags_gtdbtk", 
        "rumen_drep", "identify_novel_mags", "add_mags_to_repo", 
        "build_kraken_db", "process_novel_mags", "eggnog_annotation", 
        "functional_analysis", "dbcan_annotation", "mags_tree", 
        "tree_visualization","repository_dereplication", "rumen_dereplication", "add_to_repository",
        "prepare_taxonomy_input", "create_taxonomy", "rename_headers",
        "add_kraken_library", "build_kraken_database", "build_bracken_database",
        "build_database_distribution", "enhanced_tree_visualization", "advanced_visualizations"
    }
    
    if args.step not in no_samples_steps and not args.samples:
        raise ValueError(f"The step '{args.step}' requires samples but none were provided.")


def create_step_mapping():
    """Create mapping between step names and their corresponding functions."""
    return {
        # Quality Control & Preprocessing
        "qc": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
              qc.run(sample, cpus, memory, time, project_input, project_output),

        "multiqc": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                   multiqc.run(samples, cpus, memory, time, project_input, project_output),

        "trimming": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                   trimming.run(sample, cpus, memory, time, project_input, project_output),
        
        "host_removal": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                       host_removal.run(sample, cpus, memory, time, kwargs['reference'], project_input, project_output),
        
        # Assembly & Assessment
        "single_assembly": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                          assembly.run_idba(sample, cpus, memory, time, project_input, project_output),
        
        "co_assembly": lambda samples, cpus, memory, time, project_input, project_output, **kwargs: 
                      assembly.run_megahit(samples, cpus, memory, time, project_input, project_output),
        
        "metaquast": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                    quast.run_metaquast(sample, cpus, memory, time, project_input, project_output, kwargs.get('assembly_type', 'both')),
        
        "metaquast_coassembly": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                               quast.run_metaquast_coassembly(cpus, memory, time, project_input, project_output),
        
        # Binning & Refinement
        "single_binning": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                         binning.run_binning(sample, cpus, memory, time, kwargs.get('binning_methods', []), project_input, project_output),
        
        "co_binning": lambda samples, cpus, memory, time, project_input, project_output, **kwargs: 
                     binning.run_binning(samples, cpus, memory, time, kwargs.get('binning_methods', []), project_input, project_output),
        
        "single_bin_refinement": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                                das_tool.run_das_tool_single(sample, cpus, memory, time, project_input, project_output, kwargs.get('score_threshold', 0.5)),
        
        "co_bin_refinement": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                            das_tool.run_das_tool_coassembly(cpus, memory, time, project_input, project_output, kwargs.get('score_threshold', 0.5)),
        
        "dRep": lambda samples, cpus, memory, time, project_input, project_output, **kwargs: 
               dereplicate.run(samples, cpus, memory, time, project_input, project_output),
        
        "evaluation": lambda samples, cpus, memory, time, project_input, project_output, **kwargs: 
                     evaluation.run(samples, cpus, memory, time, project_input, project_output, kwargs.get('evaluation_input_dir')),
        
        # Taxonomic Classification
        "gtdbtk": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                 gtdbtk.run_gtdbtk(cpus, memory, time, project_input, project_output),
        
        # Rumen Reference MAGs
        "rumen_refmags_download": lambda dataset_dir: 
                                 rumen_refmags_download.run(dataset_dir),
        
        "rumen_refmags_drep": lambda dataset_dir, drep_output_dir, cpus, project_config: 
                             rumen_refmags_drep.run(dataset_dir, drep_output_dir, cpus, project_config),
        
        "rumen_refmags_gtdbtk": lambda cpus, genome_dir, gtdbtk_out_dir, genome_extension, novel_output_dir, only_novel: 
                               rumen_refmags_gtdbtk.run_gtdbtk_classify(
                                   genome_dir=genome_dir, 
                                   out_dir=gtdbtk_out_dir, 
                                   extension=genome_extension, 
                                   cpus=cpus, 
                                   identify_novel=True, 
                                   novel_output_dir=novel_output_dir, 
                                   force_rerun=not only_novel
                               ),
        
        "rumen_drep": lambda ref_mags_dir, cpus, memory, time, project_input, project_output, **kwargs: 
                     rumen_drep.run(ref_mags_dir, cpus, memory, time, project_input, project_output),
        
        # Novel MAG Processing
        "process_novel_mags": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                             nmags_pipeline.run(
                                 cpus, memory, time, project_input, project_output, 
                                 kwargs['is_rumen_data'], kwargs.get('input_mags_dir'), kwargs['kraken_db_path'], 
                                 kwargs.get('rumen_ref_mags_dir'), kwargs.get('rumen_added_mags_dir'), 
                                 kwargs.get('merge_mode', False), kwargs.get('custom_gtdbtk_dir')
                             ),
        
        "identify_novel_mags": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                              novel_mags.run(cpus, memory, time, project_input, project_output, kwargs.get('custom_gtdbtk_dir')),
        
        "add_mags_to_repo": lambda cpus, project_input, project_output, **kwargs: 
                           mags_repository.add_mags_to_repo(cpus, project_input, project_output),

        # Novel MAG Processing - Step-by-Step
        "repository_dereplication": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                                   nmags_pipeline.repository_dereplication(cpus, memory, time, project_input, project_output, **kwargs),
        
        "rumen_dereplication": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                              nmags_pipeline.rumen_dereplication(cpus, memory, time, project_input, project_output, **kwargs),
        
        "add_to_repository": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                           nmags_pipeline.add_to_repository_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "prepare_taxonomy_input": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                                nmags_pipeline.prepare_taxonomy_input_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "create_taxonomy": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                         nmags_pipeline.create_taxonomy_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "rename_headers": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                        nmags_pipeline.rename_headers_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "add_kraken_library": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                            nmags_pipeline.add_kraken_library_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "build_kraken_database": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                               nmags_pipeline.build_kraken_database_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "build_bracken_database": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                                nmags_pipeline.build_bracken_database_step(cpus, memory, time, project_input, project_output, **kwargs),
        
        "build_database_distribution": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                                     nmags_pipeline.build_database_distribution_step(cpus, memory, time, project_input, project_output, **kwargs),

        # Database Construction
        "build_kraken_db": lambda cpus, kraken_db_path, project_output, **kwargs: 
                          kraken_database.add_nmags_to_kraken(project_output, cpus, kraken_db_path),
        
        # Functional Annotation
        "eggnog_annotation": lambda cpus, memory, time, project_input, project_output, mags_dir, **kwargs: 
                            eggnog.run(cpus, memory, time, project_input, project_output, mags_dir),
        
        "dbcan_annotation": lambda cpus, memory, time, project_input, project_output, mags_dir, **kwargs: 
                           dbcan.run(cpus, memory, time, project_input, project_output, mags_dir),
        
        "functional_analysis": lambda project_output, **kwargs: 
                              functional_analysis.run_functional_analysis(project_output, kwargs.get('kegg_db_file')),
        "advanced_visualizations": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                              advanced_visualizations.run(cpus, memory, time, project_input, project_output,
                                             kwargs.get('mag_lists'), kwargs.get('list_names'),
                                             kwargs.get('viz_output_suffix', ''), kwargs.get('viz_types')),

        # Abundance Analysis
        "kraken_abundance": lambda sample, cpus, memory, time, project_input, project_output, **kwargs: 
                           kraken_abundance.run(sample, cpus, memory, time, project_input, project_output, kwargs.get('kraken_db_path')),
        
        "abundance_estimation": lambda samples, cpus, memory, time, project_input, project_output, **kwargs: 
                               abundance_estimation.run(
                                   samples, cpus, memory, time, project_input, project_output, 
                                   kwargs['kraken_db_path'], kwargs['taxonomic_levels'], 
                                   kwargs['read_length'], kwargs['threshold'], kwargs.get('metadata_file')
                               ),
        
        # Phylogenetic Analysis
        "mags_tree": lambda cpus, memory, time, project_input, project_output, **kwargs: 
            mags_tree.run_mags_tree(cpus, memory, time, kwargs.get('mags_dir'), project_input, project_output, kwargs.get('outgroup_taxon')),
        
        "tree_visualization": lambda cpus, memory, time, project_input, project_output, taxonomy_source, **kwargs: 
                             phylo_tree_vis.run_tree_visualization(cpus, memory, time, project_input, project_output, taxonomy_source),

        "enhanced_tree_visualization": lambda cpus, memory, time, project_input, project_output, **kwargs: 
                                  enhanced_phylo_tree_vis.run_enhanced_tree_visualization(
                                      cpus, memory, time, project_input, project_output, 
                                      kwargs.get('taxonomy_source'), kwargs.get('novel_mags_file'), 
                                      kwargs.get('functional_annotations'), kwargs.get('annotation_type', 'auto'),
                                      kwargs.get('cazyme_annotations'), kwargs.get('cog_annotations'), 
                                      kwargs.get('kegg_annotations')
                                  ),
    }


def execute_step(args, project_config, step_mapping):
    """Execute the specified pipeline step with appropriate arguments."""
    if args.step not in step_mapping:
        raise ValueError(f"Invalid step: {args.step}. Valid steps: {list(step_mapping.keys())}.")
    
    # Extract project paths
    project_input = project_config["input_dir"]
    project_output = project_config["output_dir"]
    reference_genome = project_config.get("reference", args.reference)
    
    print(f"[INFO] Using project input: {project_input}")
    print(f"[INFO] Using project output: {project_output}")
    print(f"[INFO] Using reference genome: {reference_genome}")
    print(f"[INFO] Selected step: {args.step}")
    
    # Display samples for applicable steps
    if args.samples and args.step not in [
        "metaquast_coassembly", "gtdbtk", "rumen_refmags_download", 
        "rumen_refmags_drep", "co_bin_refinement", "rumen_refmags_gtdbtk", 
        "rumen_drep", "identify_novel_mags", "add_mags_to_repo", 
        "build_kraken_db", "process_novel_mags", "eggnog_annotation", 
        "functional_analysis", "dbcan_annotation", "mags_tree", 
        "tree_visualization",
        "repository_dereplication", "rumen_dereplication", "add_to_repository",
        "prepare_taxonomy_input", "create_taxonomy", "rename_headers",
        "add_kraken_library", "build_kraken_database", "build_bracken_database",
        "build_database_distribution", "advanced_visualizations"
    ]:
        print(f"[INFO] Processing samples: {args.samples}")
    
    # Prepare arguments dictionary
    kwargs = {
        'input_dir': args.input_dir,
        'reference': reference_genome,
        'assembly_type': args.assembly_type,
        'binning_methods': args.binning_methods,
        'score_threshold': args.score_threshold,
        'evaluation_input_dir': args.evaluation_input_dir,
        'dataset_dir': args.dataset_dir,
        'drep_output_dir': args.drep_output_dir,
        'ref_mags_dir': args.ref_mags_dir,
        'genome_dir': args.genome_dir or project_output,
        'gtdbtk_out_dir': args.gtdbtk_out_dir or f"{project_output}/gtdbtk",
        'genome_extension': args.genome_extension,
        'is_rumen_data': args.is_rumen_data,
        'input_mags_dir': args.input_mags_dir,
        'custom_gtdbtk_dir': args.custom_gtdbtk_dir,
        'novel_output_dir': args.novel_output_dir,
        'only_novel': args.only_novel,
        'rumen_ref_mags_dir': args.rumen_ref_mags_dir,
        'rumen_added_mags_dir': args.rumen_added_mags_dir,
        'merge_mode': args.merge_mode,
        'kraken_db_path': args.kraken_db_path,
        'kegg_db_file': args.kegg_db_file,
        'mags_dir': args.mags_dir,
        'outgroup_taxon': args.outgroup_taxon,
        'taxonomy_source': args.taxonomy_source,
        'taxonomic_levels': args.taxonomic_levels,
        'read_length': args.read_length,
        'threshold': args.threshold,
        'metadata_file': args.metadata_file,
        'project_config': args.project_config,
        'novel_mags_file': args.novel_mags_file,
        'functional_annotations': args.functional_annotations, 
        'annotation_type': args.annotation_type,
        'cazyme_annotations': args.cazyme_annotations,
        'cog_annotations': args.cog_annotations,
        'kegg_annotations': args.kegg_annotations,
        'mag_lists': args.mag_lists,
        'list_names': args.list_names,
        'viz_output_suffix': args.viz_output_suffix,
        'viz_types': args.viz_types,
    }
    
    # Get the step function
    step_func = step_mapping[args.step]
    
    # Execute based on step requirements
    if args.step in ["co_assembly", "co_binning", "dRep", "evaluation", "abundance_estimation", "mutliqc"]:
        # Steps that process all samples together
        step_func(args.samples, args.cpus, args.memory, args.time, project_input, project_output, **kwargs)
    
    elif args.step in [
        "metaquast_coassembly", "co_bin_refinement", "gtdbtk", 
        "identify_novel_mags", "process_novel_mags", "eggnog_annotation", 
        "dbcan_annotation", "mags_tree", "tree_visualization","repository_dereplication", "rumen_dereplication", "add_to_repository",
        "prepare_taxonomy_input", "create_taxonomy", "rename_headers",
        "add_kraken_library", "build_kraken_database", "build_bracken_database",
        "build_database_distribution","enhanced_tree_visualization", "advanced_visualizations"
    ]:
        # Steps that don't require samples but need resource parameters
        step_func(args.cpus, args.memory, args.time, project_input, project_output, **kwargs)
    
    elif args.step in ["rumen_refmags_download"]:
        # Steps with minimal arguments
        step_func(args.dataset_dir)
    
    elif args.step in ["rumen_refmags_drep", "rumen_drep"]:
        # Steps with specific argument patterns
        if args.step == "rumen_refmags_drep":
            step_func(args.dataset_dir, args.drep_output_dir, args.cpus, args.project_config)
        else:
            step_func(args.ref_mags_dir, args.cpus, args.memory, args.time, project_input, project_output, **kwargs)

    elif args.step == "rumen_refmags_gtdbtk":
        # GTDB-Tk specific execution
        step_func(args.cpus, args.genome_dir or project_output, 
                 args.gtdbtk_out_dir or f"{project_output}/gtdbtk", 
                 args.genome_extension, args.novel_output_dir, args.only_novel)
    
    elif args.step == "functional_analysis":
        # Functional analysis execution
        step_func(project_output, **kwargs)
    
    elif args.step == "add_mags_to_repo":
        # MAGs repository execution
        step_func(args.cpus, project_input, project_output, **kwargs)
    
    elif args.step == "build_kraken_db":
        # Kraken database building
        step_func(args.cpus, args.kraken_db_path, project_output, **kwargs)
    
    elif args.step in ["metaquast", "single_bin_refinement", "kraken_abundance"]:
        # Steps that process each sample individually
        for sample in args.samples:
            if args.step == "kraken_abundance":
                step_func(sample, args.cpus, args.memory, args.time, project_input, project_output, **kwargs)
            else:
                step_func(sample, args.cpus, args.memory, args.time, project_input, project_output, **kwargs)
    
    else:
        # Default: process each sample individually
        for sample in args.samples:
            step_func(sample, args.cpus, args.memory, args.time, project_input, project_output, **kwargs)


def main():
    """Main function to orchestrate pipeline step execution."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Load project configuration
    project_config = load_project_config(args.project_config)
    
    # Validate step requirements
    validate_step_requirements(args)
    
    # Create step mapping
    step_mapping = create_step_mapping()
    
    # Execute the specified step
    execute_step(args, project_config, step_mapping)


if __name__ == "__main__":
    main()

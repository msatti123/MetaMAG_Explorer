#!/usr/bin/env python3
"""
Script for running GTDB-Tk classification workflow and identifying novel MAGs.
"""

import os
import sys
import argparse
import subprocess
import shutil
import pandas as pd
import yaml

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    os.makedirs(directory, exist_ok=True)

def get_tool_path(tool_name):
    """Get tool path from config.py - handles multiple import scenarios"""
    try:
        # Method 1: Try to import config from MetaMAG package
        from MetaMAG.config import config
        return config["tools"].get(tool_name)
    except ImportError:
        try:
            # Method 2: Fallback to import from current directory
            from config import config
            return config["tools"].get(tool_name)
        except ImportError:
            try:
                # Method 3: Try importing from parent directory
                sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
                from config import config
                return config["tools"].get(tool_name)
            except ImportError:
                print(f"[WARNING] Could not import config module for tool: {tool_name}")
                return None
    except (KeyError, TypeError) as e:
        print(f"[WARNING] Error accessing config for {tool_name}: {e}")
        return None

def yaml_get_tool_path(yaml_config, tool_name):
    """Get tool path from YAML config - safe extraction"""
    try:
        if yaml_config and isinstance(yaml_config, dict):
            tools = yaml_config.get("tools", {})
            if isinstance(tools, dict):
                return tools.get(tool_name)
        return None
    except (KeyError, TypeError, AttributeError) as e:
        print(f"[WARNING] Error accessing YAML config for {tool_name}: {e}")
        return None

def validate_gtdbtk_availability(gtdbtk_activate=None):
    """Validate that GTDB-Tk is available either via activation script or system PATH"""
    if gtdbtk_activate:
        # Split the activation command if it contains spaces
        parts = gtdbtk_activate.split()
        if len(parts) >= 1:
            activation_script = parts[0]
            env_name = parts[1] if len(parts) > 1 else "base"
            
            # Check if activation script exists
            if os.path.exists(activation_script):
                print(f"[INFO] Conda activation script found: {activation_script}")
                return f"{activation_script} {env_name}", True
            else:
                print(f"[WARNING] Activation script not found: {activation_script}")
    
    # Check if GTDB-Tk is available in system PATH
    try:
        result = subprocess.run(
            ["gtdbtk", "--version"], 
            check=True, 
            capture_output=True, 
            text=True,
            timeout=10
        )
        print(f"[INFO] System GTDB-Tk found: {result.stdout.strip()}")
        return None, True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        print(f"[ERROR] GTDB-Tk not found in system PATH")
        return None, False

def check_gtdbtk_completed(out_dir):
    """
    Check if GTDB-Tk has already been run successfully.
    
    Args:
        out_dir (str): GTDB-Tk output directory
        
    Returns:
        bool: True if GTDB-Tk has been completed, False otherwise
    """
    # Check for key output files
    bac_summary = os.path.join(out_dir, "classify", "gtdbtk.bac120.summary.tsv")
    ar_summary = os.path.join(out_dir, "classify", "gtdbtk.ar53.summary.tsv")
    
    # Consider GTDB-Tk completed if at least one of the summary files exists
    if os.path.exists(bac_summary) or os.path.exists(ar_summary):
        if os.path.exists(bac_summary):
            print(f"[INFO] Found existing bacterial classification results: {bac_summary}")
        if os.path.exists(ar_summary):
            print(f"[INFO] Found existing archaeal classification results: {ar_summary}")
        return True
    
    return False

def identify_novel_mags(taxonomy_dir, output_dir, source_mags_dir, extension=".fa"):
    """
    Identify novel MAGs from GTDB-Tk classification results and organize them.
    
    Args:
        taxonomy_dir (str): Directory containing GTDB-Tk classification results
        output_dir (str): Base directory for novel MAGs output
        source_mags_dir (str): Directory containing source MAG files
        extension (str): File extension of MAG files
        
    Returns:
        int: 0 for success, 1 for failure
    """
    print(f"[INFO] Identifying novel MAGs from GTDB-Tk results...")
    print(f"[INFO] Taxonomy directory: {taxonomy_dir}")
    print(f"[INFO] Output directory: {output_dir}")
    print(f"[INFO] Source MAGs directory: {source_mags_dir}")
    
    # Create output directories
    mags_dir = os.path.join(output_dir, "mags")
    taxonomy_dir_output = os.path.join(output_dir, "taxonomy")
    
    ensure_directory_exists(mags_dir)
    ensure_directory_exists(taxonomy_dir_output)
    
    # Find GTDB-Tk summary files
    bac_summary = os.path.join(taxonomy_dir, "classify", "gtdbtk.bac120.summary.tsv")
    ar_summary = os.path.join(taxonomy_dir, "classify", "gtdbtk.ar53.summary.tsv")
    
    novel_bac_mags = []
    novel_ar_mags = []
    
    # Process bacterial summary file if it exists
    if os.path.exists(bac_summary):
        print(f"[INFO] Processing bacterial classification summary: {bac_summary}")
        try:
            bac_df = pd.read_csv(bac_summary, sep='\t')
            print(f"[INFO] Found {len(bac_df)} bacterial MAGs in classification results")
            
            # Find MAGs with missing species designation (s__)
            novel_bac_mags = bac_df[bac_df['classification'].str.split(';').str[-1].str.strip() == 's__']['user_genome'].tolist()
            
            if novel_bac_mags:
                # Filter dataframe to only include novel MAGs
                novel_bac_df = bac_df[bac_df['user_genome'].isin(novel_bac_mags)]
                
                # Save novel bacterial taxonomy to output directory
                bac_output_file = os.path.join(taxonomy_dir_output, "novel_bac_taxonomy.tsv")
                novel_bac_df.to_csv(bac_output_file, sep='\t', index=False)
                print(f"[INFO] Found {len(novel_bac_mags)} novel bacterial MAGs")
                print(f"[INFO] Saved novel bacterial taxonomy to: {bac_output_file}")
            else:
                print("[INFO] No novel bacterial MAGs found")
        except Exception as e:
            print(f"[ERROR] Failed to process bacterial summary: {e}")
            return 1
    else:
        print(f"[WARNING] Bacterial summary file not found: {bac_summary}")
    
    # Process archaeal summary file if it exists
    if os.path.exists(ar_summary):
        print(f"[INFO] Processing archaeal classification summary: {ar_summary}")
        try:
            ar_df = pd.read_csv(ar_summary, sep='\t')
            print(f"[INFO] Found {len(ar_df)} archaeal MAGs in classification results")
            
            # Find MAGs with missing species designation (s__)
            novel_ar_mags = ar_df[ar_df['classification'].str.split(';').str[-1].str.strip() == 's__']['user_genome'].tolist()
            
            if novel_ar_mags:
                # Filter dataframe to only include novel MAGs
                novel_ar_df = ar_df[ar_df['user_genome'].isin(novel_ar_mags)]
                
                # Save novel archaeal taxonomy to output directory
                ar_output_file = os.path.join(taxonomy_dir_output, "novel_ar_taxonomy.tsv")
                novel_ar_df.to_csv(ar_output_file, sep='\t', index=False)
                print(f"[INFO] Found {len(novel_ar_mags)} novel archaeal MAGs")
                print(f"[INFO] Saved novel archaeal taxonomy to: {ar_output_file}")
            else:
                print("[INFO] No novel archaeal MAGs found")
        except Exception as e:
            print(f"[ERROR] Failed to process archaeal summary: {e}")
            return 1
    else:
        print(f"[INFO] Archaeal summary file not found: {ar_summary}")
    
    # Combine novel MAGs from both domains
    novel_mags = novel_bac_mags + novel_ar_mags
    
    if not novel_mags:
        print("[WARNING] No novel MAGs found in either domain")
        return 0
    
    print(f"[INFO] Total novel reference MAGs identified: {len(novel_mags)}")
    
    # Create a combined list of all novel MAGs
    try:
        novel_list_file = os.path.join(taxonomy_dir_output, "novel_mags_list.txt")
        with open(novel_list_file, 'w') as f:
            for mag in novel_mags:
                f.write(f"{mag}\n")
        print(f"[INFO] Saved list of all novel MAGs to: {novel_list_file}")
    except Exception as e:
        print(f"[WARNING] Failed to save combined MAG list: {e}")
    
    # Copy MAG files to output directory
    copied_count = 0
    for mag in novel_mags:
        # Try different extensions if not found
        extensions_to_try = [extension]
        if extension != '.fa':
            extensions_to_try.extend(['.fa', '.fasta', '.fna'])
        
        source_file = None
        for ext in extensions_to_try:
            potential_file = os.path.join(source_mags_dir, f"{mag}{ext}")
            if os.path.exists(potential_file):
                source_file = potential_file
                break
        
        if source_file:
            dest_file = os.path.join(mags_dir, f"{mag}{extension}")
            try:
                shutil.copy2(source_file, dest_file)
                copied_count += 1
            except Exception as e:
                print(f"[ERROR] Failed to copy {mag}: {e}")
        else:
            print(f"[WARNING] Source MAG file not found for {mag} with any extension")
    
    print(f"[INFO] Successfully copied {copied_count}/{len(novel_mags)} novel MAGs to {mags_dir}")
    
    return 0

def run_gtdbtk_classify(genome_dir, out_dir, extension='.fa', cpus=64, config_file=None, identify_novel=True, novel_output_dir=None, force_rerun=False):
    """
    Main function to run GTDB-Tk classification workflow and optionally identify novel MAGs.
    
    Args:
        genome_dir (str): Directory containing input genome files
        out_dir (str): Directory to store GTDB-Tk output
        extension (str, optional): File extension of genomes. Defaults to '.fa'
        cpus (int, optional): Number of CPUs to use. Defaults to 64
        config_file (str, optional): Path to additional configuration file
        identify_novel (bool, optional): Whether to identify novel MAGs. Defaults to True.
        novel_output_dir (str, optional): Output directory for novel MAGs. Defaults to None.
        force_rerun (bool, optional): Force rerun of GTDB-Tk even if output exists. Defaults to False.
    
    Returns:
        int: 0 for success, 1 for failure
    """
    print("=" * 70)
    print("Starting GTDB-Tk Classification and Novel MAG Identification")
    print("=" * 70)
    print(f"Genome directory: {genome_dir}")
    print(f"Output directory: {out_dir}")
    print(f"File extension: {extension}")
    print(f"CPUs: {cpus}")
    print(f"Config file: {config_file if config_file else 'None'}")
    print(f"Identify novel: {identify_novel}")
    print(f"Force rerun: {force_rerun}")
    print("=" * 70)
    
    # Ensure output directory exists
    ensure_directory_exists(out_dir)
    
    # Check if GTDB-Tk has already been run
    gtdbtk_completed = check_gtdbtk_completed(out_dir)
    
    # Skip GTDB-Tk if already completed and not forced to rerun
    if gtdbtk_completed and not force_rerun:
        print("[INFO] GTDB-Tk has already been run. Skipping classification step.")
    else:
        print("[INFO] Running GTDB-Tk classification...")
        
        # Get GTDB-Tk activation path
        gtdbtk_activate = None
        
        # Check config file if provided
        if config_file and os.path.exists(config_file):
            print(f"[INFO] Loading configuration from: {config_file}")
            try:
                with open(config_file, 'r') as f:
                    yaml_config = yaml.safe_load(f)
                    if yaml_config and 'tools' in yaml_config:
                        gtdbtk_activate = yaml_get_tool_path(yaml_config, "gtdbtk")
                        if gtdbtk_activate:
                            print(f"[INFO] Found gtdbtk_activate in YAML config: {gtdbtk_activate}")
            except Exception as e:
                print(f"[WARNING] Error loading config file {config_file}: {e}")
        
        # Fallback to config.py if no activation path found
        if not gtdbtk_activate:
            print("[INFO] Trying to get gtdbtk_activate from config.py")
            gtdbtk_activate = get_tool_path("gtdbtk")
            if gtdbtk_activate:
                print(f"[INFO] Found gtdbtk_activate in config.py: {gtdbtk_activate}")
        
        # Validate GTDB-Tk availability
        validated_activate, gtdbtk_available = validate_gtdbtk_availability(gtdbtk_activate)
        if not gtdbtk_available:
            print("[ERROR] GTDB-Tk is not available via activation script or system PATH")
            print("[ERROR] Please install GTDB-Tk or configure the activation script path")
            return 1
        
        # Verify genome directory exists
        if not os.path.exists(genome_dir):
            print(f"[ERROR] Genome directory does not exist: {genome_dir}")
            return 1
        
        # Find genome files
        genome_files = [f for f in os.listdir(genome_dir) if f.endswith(extension)]
        
        if not genome_files:
            print(f"[ERROR] No {extension} files found in {genome_dir}")
            return 1
        
        print(f"[INFO] Found {len(genome_files)} genome files with extension {extension}")
        
        # Prepare GTDB-Tk command
        gtdbtk_cmd = (
            f"gtdbtk classify_wf "
            f"--genome_dir {genome_dir} "
            f"--out_dir {out_dir} "
            f"--extension {extension} "
            f"--skip_ani_screen "
            f"--cpus {cpus}"
        )
        
        # If activation path is provided, wrap command in activation
        if validated_activate:
            print(f"[INFO] Using GTDB-Tk activation script: {validated_activate}")
            full_cmd = f"bash -c 'source {validated_activate} && {gtdbtk_cmd}'"
        else:
            print("[INFO] Using system GTDB-Tk installation")
            full_cmd = gtdbtk_cmd
        
        print(f"[INFO] GTDB-Tk command: {full_cmd}")
        
        # Run GTDB-Tk
        try:
            print("[INFO] Starting GTDB-Tk classification workflow...")
            result = subprocess.run(
                full_cmd, 
                shell=True, 
                check=True, 
                stderr=subprocess.PIPE, 
                stdout=subprocess.PIPE,
                text=True,
                timeout=14400  # 4 hour timeout
            )
            print("[INFO] GTDB-Tk classification completed successfully")
            if result.stdout:
                print(f"[INFO] GTDB-Tk stdout: {result.stdout}")
            gtdbtk_completed = True
            
        except subprocess.TimeoutExpired:
            print("[ERROR] GTDB-Tk execution timed out (4 hours)")
            return 1
        except subprocess.CalledProcessError as e:
            print("[ERROR] GTDB-Tk execution failed")
            print(f"[ERROR] Return code: {e.returncode}")
            if e.stdout:
                print(f"[ERROR] Command output: {e.stdout}")
            if e.stderr:
                print(f"[ERROR] Command error: {e.stderr}")
            return 1
        except Exception as e:
            print(f"[ERROR] Unexpected error during GTDB-Tk execution: {e}")
            return 1
    
    # Identify novel MAGs if requested and GTDB-Tk has completed
    if identify_novel and gtdbtk_completed:
        if not novel_output_dir:
            # Default to 'novel_rumenref_mags' in the parent directory of out_dir
            parent_dir = os.path.dirname(out_dir)
            novel_output_dir = os.path.join(parent_dir, "novel_rumenref_mags")
        
        print(f"[INFO] Identifying novel MAGs and saving to {novel_output_dir}")
        ensure_directory_exists(novel_output_dir)
        novel_exit_code = identify_novel_mags(out_dir, novel_output_dir, genome_dir, extension)
        
        if novel_exit_code != 0:
            print("[ERROR] Failed to identify novel MAGs")
            return novel_exit_code
    
    print("=" * 70)
    print("? GTDB-Tk workflow completed successfully")
    print("=" * 70)
    
    return 0

def main():
    """Command line interface for GTDB-Tk classification and novel MAG identification"""
    parser = argparse.ArgumentParser(
        description="Run GTDB-Tk classification workflow and identify novel MAGs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full workflow
  python rumen_refmags_gtdbtk.py --genome_dir /path/to/genomes --gtdbtk_out_dir /path/to/output --cpus 32
  
  # Only identify novel MAGs (skip GTDB-Tk if already run)
  python rumen_refmags_gtdbtk.py --genome_dir /path/to/genomes --gtdbtk_out_dir /path/to/output --only_novel
  
  # Force rerun GTDB-Tk even if output exists
  python rumen_refmags_gtdbtk.py --genome_dir /path/to/genomes --gtdbtk_out_dir /path/to/output --force_rerun
        """
    )
    
    parser.add_argument("--genome_dir", type=str, required=True, 
                        help="Directory containing input genomes")
    parser.add_argument("--gtdbtk_out_dir", type=str, required=True, 
                        help="Output directory for GTDB-Tk results")
    parser.add_argument("--genome_extension", type=str, default=".fa", 
                        help="File extension of input genomes (default: .fa)")
    parser.add_argument("--cpus", type=int, default=64, 
                        help="Number of CPUs to use (default: 64)")
    parser.add_argument("--config", type=str, help="Path to additional configuration file")
    parser.add_argument("--novel_output_dir", type=str, 
                        help="Output directory for novel MAGs (default: novel_rumenref_mags in the parent directory of gtdbtk_out_dir)")
    parser.add_argument("--skip_novel", action="store_true", 
                        help="Skip novel MAG identification")
    parser.add_argument("--force_rerun", action="store_true", 
                        help="Force rerun of GTDB-Tk even if output exists")
    parser.add_argument("--only_novel", action="store_true", 
                        help="Only identify novel MAGs (skip GTDB-Tk if already run)")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.cpus < 1:
        print("[ERROR] Number of CPUs must be at least 1")
        sys.exit(1)
    
    # If only_novel is specified, don't force rerun of GTDB-Tk
    if args.only_novel:
        args.force_rerun = False
    
    # Call run function with parsed arguments
    exit_code = run_gtdbtk_classify(
        genome_dir=args.genome_dir, 
        out_dir=args.gtdbtk_out_dir, 
        extension=args.genome_extension,
        cpus=args.cpus,
        config_file=args.config,
        identify_novel=not args.skip_novel,
        novel_output_dir=args.novel_output_dir,
        force_rerun=args.force_rerun
    )
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
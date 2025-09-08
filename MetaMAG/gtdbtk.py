import os
import subprocess
import argparse
import sys
from MetaMAG.config import config, get_tool_path, get_tool_command
from MetaMAG.utils import ensure_directory_exists, run_command

def run_gtdbtk_classify(genome_dir, out_dir, cpus, extension=".fa"):
    """
    Run GTDB-Tk classify workflow on genomes.
    
    Args:
        genome_dir (str): Directory containing input genomes
        out_dir (str): Output directory for GTDB-Tk results
        cpus (int): Number of CPUs to use
        extension (str): File extension of input genomes
    
    Returns:
        bool: Success or failure
    """
    # Create output directory
    ensure_directory_exists(out_dir)
    
    # Check if genome directory exists
    if not os.path.exists(genome_dir):
        print(f"[ERROR] Genome directory not found: {genome_dir}")
        return False
    
    # Find genome files
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(extension)]
    if not genome_files:
        print(f"[ERROR] No files with extension {extension} found in {genome_dir}")
        return False
    
    print(f"[INFO] Found {len(genome_files)} genome files with extension {extension}")
    
    # Check if output already exists
    bac_summary = os.path.join(out_dir, "classify", "gtdbtk.bac120.summary.tsv")
    ar_summary = os.path.join(out_dir, "classify", "gtdbtk.ar53.summary.tsv")
    
    if os.path.exists(bac_summary) or os.path.exists(ar_summary):
        print("[INFO] GTDB-Tk classification results already exist. Skipping classification.")
        if os.path.exists(bac_summary):
            print(f"[INFO] Bacterial summary exists: {bac_summary}")
        if os.path.exists(ar_summary):
            print(f"[INFO] Archaeal summary exists: {ar_summary}")
        return True
    
    # Get GTDB-Tk activation command from config
    gtdbtk_activate = get_tool_path("gtdbtk")
    if not gtdbtk_activate:
        print(f"[ERROR] Gtdbtk_Activate tool not found in configuration")
        return
    
    # Create the log file path
    log_file = os.path.join(out_dir, "gtdbtk_output.log")
    
    # Build the GTDB-Tk command - only the classify workflow
    gtdbtk_command = (
        f"bash -c 'source {gtdbtk_activate} && "
        f"gtdbtk classify_wf --genome_dir {genome_dir} --out_dir {out_dir} "
        f"--extension {extension} --skip_ani_screen --cpus {cpus} > {log_file} 2>&1'"
    )
    
    print(f"[INFO] Running GTDB-Tk classify_wf...")
    print(f"[DEBUG] Command: {gtdbtk_command}")
    
    # Use the run_command utility from your existing code
    try:
        run_command(gtdbtk_command)
        print(f"[INFO] GTDB-Tk classification completed successfully")
        
        # Check if classification produced expected output files
        if os.path.exists(bac_summary):
            print(f"[INFO] Bacterial classification summary created: {bac_summary}")
        if os.path.exists(ar_summary):
            print(f"[INFO] Archaeal classification summary created: {ar_summary}")
            
        return True
    except Exception as e:
        print(f"[ERROR] GTDB-Tk execution failed: {e}")
        return False

def run_gtdbtk(cpus, memory, time, project_input, project_output):
    """
    Main function to run GTDB-Tk.
    
    Args:
        cpus (int): Number of CPUs to allocate
        memory (str): Memory allocation
        time (str): Time allocation
        project_input (str): Project input directory
        project_output (str): Project output directory
    
    Returns:
        bool: Success or failure
    """
    # Define fixed paths for input and output
    genome_dir = os.path.join(project_output, "Bin_Refinement", "drep", "dRep_output", "dereplicated_genomes")
    gtdbtk_out_dir = os.path.join(project_output, "Novel_Mags", "gtdbtk")
    extension = ".fa"
    
    print(f"[INFO] Using MAGs from: {genome_dir}")
    print(f"[INFO] Using output directory: {gtdbtk_out_dir}")
    print(f"[INFO] Using {cpus} CPUs for GTDB-Tk classification")
    
    return run_gtdbtk_classify(
        genome_dir=genome_dir,
        out_dir=gtdbtk_out_dir,
        cpus=cpus,
        extension=extension
    )

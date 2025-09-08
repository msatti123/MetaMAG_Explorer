#!/usr/bin/env python3
"""
Module for running Kraken2 abundance estimation on metagenomic samples.
This module handles taxonomic classification using Kraken2.
"""

import os
import subprocess
import sys
from MetaMAG.config import config, get_tool_path, get_tool_command

def run(sample, cpus, memory, time, project_input, project_output, kraken_db_path):
    """
    Run Kraken2 abundance estimation on a sample.
    
    Args:
        sample (str): Sample ID
        cpus (int): Number of CPUs to use
        memory (str): Memory allocation
        time (str): Time allocation
        project_input (str): Path to project input directory
        project_output (str): Path to project output directory
        kraken_db_path (str): Path to Kraken2 database
    
    Returns:
        bool: True if successful, False otherwise
    """
    print(f"[INFO] Running Kraken2 abundance estimation for sample: {sample}")
    
    # Get paths from config
    perl_path = config["environment"].get("PERL5LIB", "")
    kraken2_path = get_tool_path("kraken2")
    if not kraken2_path:
        print(f"[ERROR] Kraken2 tool not found in configuration")
        return
    
    if not perl_path:
        print("[ERROR] PERL5LIB not found in config. Please add it to the config file.")
        return False
    
    if not kraken2_path:
        print("[ERROR] Kraken2 path not found in config. Please add it to the config file.")
        return False
    
    if not kraken_db_path:
        print("[ERROR] Kraken2 database path not provided.")
        return False
    
    # Set up paths for input fastq files (after host removal)
    host_removal_dir = os.path.join(project_output, "Host_Removal")
    r1_file = os.path.join(host_removal_dir, f"{sample}_unmapped_R1.fastq.gz")
    r2_file = os.path.join(host_removal_dir, f"{sample}_unmapped_R2.fastq.gz")
    
    # Check if input files exist
    if not os.path.exists(r1_file) or not os.path.exists(r2_file):
        print(f"[ERROR] Input files not found for sample {sample}.")
        print(f"Expected: {r1_file} and {r2_file}")
        return False
    
    # Create output directory
    output_dir = os.path.join(project_output, "Kraken_Abundance")
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output files
    kraken_output = os.path.join(output_dir, f"{sample}_kraken.txt")
    kraken_report = os.path.join(output_dir, f"{sample}_kreport.txt")
    
    # Construct Kraken2 command
    cmd = f"PERL5LIB={perl_path} {kraken2_path} --paired --gzip-compressed --db {kraken_db_path} " \
          f"--threads {cpus} --report {kraken_report} --output {kraken_output} {r1_file} {r2_file}"
    
    print(f"[INFO] Running command: {cmd}")
    
    try:
        # Execute Kraken2 command
        subprocess.run(cmd, shell=True, check=True)
        print(f"[INFO] Kraken2 abundance estimation completed for sample {sample}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to run Kraken2 for sample {sample}: {e}")
        return False

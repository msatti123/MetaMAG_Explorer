#!/usr/bin/env python3
"""
Module for running EggNOG-mapper on MAGs.
This module handles the annotation of MAGs using EggNOG-mapper.
"""

import os
import subprocess
import yaml
import shutil
import time
import glob
from MetaMAG.config import config, get_tool_path, get_tool_command
from MetaMAG.functional_analysis import run_functional_analysis

def run(cpus, memory, time_limit, project_input, project_output, mags_dir=None):
    """
    Run EggNOG-mapper on the provided MAGs directory.
    
    Parameters:
    -----------
    cpus : int
        Number of CPUs to use for annotation
    memory : str
        Memory allocation (e.g., "200G")
    time_limit : str
        Time allocation (e.g., "30-12:00:00")
    project_input : str
        Path to project input directory
    project_output : str
        Path to project output directory
    mags_dir : str, optional
        Path to MAGs directory. If None, uses default from project configuration.
    """
    print(f"[INFO] Starting EggNOG-mapper annotation")
    print(f"[INFO] Using {cpus} CPUs")
    
    # Set up input and output directories
    if not mags_dir:
        mags_dir = os.path.join(project_output, "final_mags")
    
    if not os.path.exists(mags_dir):
        print(f"[ERROR] MAGs directory not found: {mags_dir}")
        return

    # Create output directory
    eggnog_output_dir = os.path.join(project_output, "eggnog_output")
    os.makedirs(eggnog_output_dir, exist_ok=True)
    
    # Get EggNOG-mapper paths from config
    eggnog_mapper_path = get_tool_path("eggnog_mapper")
    if not eggnog_mapper_path:
        print(f"[ERROR] Eggnog_Mapper tool not found in configuration")
        return
    eggnog_db_dir = get_tool_path("eggnog_db_dir")
    if not eggnog_db_dir:
        print(f"[ERROR] Eggnog_Db_Dir tool not found in configuration")
        return
    prodigal_path = get_tool_path("prodigal")
    if not prodigal_path:
        print(f"[ERROR] Prodigal tool not found in configuration")
        return  # Default to "prodigal" if not specified
    
    if not eggnog_mapper_path:
        print(f"[ERROR] EggNOG-mapper path not found in config. Make sure 'eggnog_mapper' is defined in the tools section.")
        return
    
    if not os.path.exists(eggnog_mapper_path):
        print(f"[ERROR] EggNOG-mapper path does not exist: {eggnog_mapper_path}")
        return
    
    if not eggnog_db_dir:
        print(f"[ERROR] EggNOG database directory not found in config. Make sure 'eggnog_db_dir' is defined in the tools section.")
        return
    
    if not os.path.exists(eggnog_db_dir):
        print(f"[ERROR] EggNOG database directory does not exist: {eggnog_db_dir}")
        return

    # Find all MAG files
    mag_files = []
    for root, dirs, files in os.walk(mags_dir):
        for file in files:
            if file.endswith((".fa", ".fasta", ".fna")):
                mag_files.append(os.path.join(root, file))
    
    # Check if all MAGs have already been processed by checking annotation files
    need_processing = False
    mags_to_process = []
    
    for mag_file in mag_files:
        mag_base_name = os.path.basename(mag_file).replace(".", "_")
        
        # Check if annotation file exists
        annotation_pattern = os.path.join(eggnog_output_dir, f"{mag_base_name}*.emapper.annotations")
        annotation_files = glob.glob(annotation_pattern)
        
        # Also check for protein base names
        protein_base_name = os.path.basename(mag_file) + ".faa"
        protein_base_name = protein_base_name.replace(".", "_")
        protein_annotation_pattern = os.path.join(eggnog_output_dir, f"{protein_base_name}*.emapper.annotations")
        protein_annotation_files = glob.glob(protein_annotation_pattern)
        
        # Combine both potential matches
        all_annotation_files = annotation_files + protein_annotation_files
        
        if not all_annotation_files:
            need_processing = True
            mags_to_process.append(mag_file)
    
    if not need_processing:
        print(f"[INFO] All {len(mag_files)} MAGs already have annotation files. Skipping annotation step.")
        print(f"[INFO] Extracting functional annotations from existing files.")
        extract_functional_annotations(eggnog_output_dir, project_output)
        
        # Run the functional analysis
        run_functional_analysis(project_output)
        return
    else:
        print(f"[INFO] Found {len(mags_to_process)} MAGs that need processing out of {len(mag_files)} total MAGs.")
    
    # Create temporary directory for EggNOG database
    tmp_db_dir = "/dev/shm/eggnogdb"
    os.makedirs(tmp_db_dir, exist_ok=True)
    
    try:
        # Copy the pre-prepared EggNOG database files
        print(f"[INFO] Copying EggNOG database files from {eggnog_db_dir} to {tmp_db_dir}")
        copy_cmd = f"cp -r {eggnog_db_dir}/* {tmp_db_dir}/"
        subprocess.run(copy_cmd, shell=True, check=True)
        
        # Check if the diamond database exists in the copied files
        diamond_db = os.path.join(tmp_db_dir, "eggnog_proteins.dmnd")
        if not os.path.exists(diamond_db):
            print(f"[WARNING] Diamond database not found in copied files: {diamond_db}")
            print(f"[WARNING] Will use HMMER instead of Diamond")
            search_mode = "hmmer"
        else:
            search_mode = "diamond"
            print(f"[INFO] Diamond database found: {diamond_db}")
        
        # Step 1: First look for protein files for MAGs that need processing
        protein_files_map = {}  # Map MAG file to its protein file

        # Check for existing protein files in MAGs directory
        for mag_file in mags_to_process:
            mag_base = os.path.basename(mag_file)
            mag_dir = os.path.dirname(mag_file)
            
            # Common naming patterns for protein files
            potential_protein_names = [
                f"{mag_base}.faa",
                f"{mag_base}.proteins.fa",
                f"{mag_base}.pep",
                f"{os.path.splitext(mag_base)[0]}.faa",
                f"{os.path.splitext(mag_base)[0]}.proteins.fa"
            ]
            
            # Check MAG directory for protein files
            protein_found = False
            for protein_name in potential_protein_names:
                protein_path = os.path.join(mag_dir, protein_name)
                if os.path.exists(protein_path):
                    protein_files_map[mag_file] = protein_path
                    protein_found = True
                    break
            
            # Check output directory for existing protein files
            if not protein_found:
                for protein_name in potential_protein_names:
                    protein_path = os.path.join(eggnog_output_dir, protein_name)
                    if os.path.exists(protein_path):
                        protein_files_map[mag_file] = protein_path
                        protein_found = True
                        break
        
        # Step 2: For MAGs without protein files, predict them
        for mag_file in mags_to_process:
            if mag_file not in protein_files_map:
                protein_file = os.path.join(eggnog_output_dir, os.path.basename(mag_file) + ".faa")
                
                # Skip if protein prediction already exists
                if os.path.exists(protein_file):
                    print(f"[INFO] Protein file already exists: {protein_file}")
                    protein_files_map[mag_file] = protein_file
                    continue
                    
                print(f"[INFO] Predicting proteins from {mag_file}")
                prodigal_cmd = f"{prodigal_path} -i {mag_file} -a {protein_file} -p meta"
                subprocess.run(prodigal_cmd, shell=True, check=True)
                protein_files_map[mag_file] = protein_file
        
        # Step 3: Run EggNOG-mapper on each protein file
        for mag_file, protein_file in protein_files_map.items():
            # Get bases for checking existing annotations
            mag_base = os.path.basename(mag_file).replace(".", "_")
            protein_base = os.path.basename(protein_file).replace(".", "_")
            
            # Check if annotation already exists using direct or derived name
            annotation_exists = False
            
            # Check for direct protein file annotation
            annotation_pattern = os.path.join(eggnog_output_dir, f"{protein_base}*.emapper.annotations")
            if glob.glob(annotation_pattern):
                print(f"[INFO] Annotation already exists for {protein_file}. Skipping.")
                annotation_exists = True
            
            # Check for MAG-derived annotation 
            mag_annotation_pattern = os.path.join(eggnog_output_dir, f"{mag_base}*.emapper.annotations")
            if glob.glob(mag_annotation_pattern):
                print(f"[INFO] Annotation already exists for {mag_file}. Skipping.")
                annotation_exists = True
                
            if annotation_exists:
                continue
                
            # Run EggNOG-mapper
            print(f"[INFO] Processing {protein_file} with EggNOG-mapper")
            
            # Set output prefix based on protein filename
            output_prefix = os.path.join(eggnog_output_dir, protein_base)
            
            # Run EggNOG-mapper with usemem option for faster processing
            eggnog_cmd = f"{eggnog_mapper_path} -i {protein_file} --data_dir {tmp_db_dir} --dbmem --usemem -o {output_prefix} --output_dir {eggnog_output_dir} --cpu {cpus} -m {search_mode}"
            
            print(f"[INFO] Running command: {eggnog_cmd}")
            subprocess.run(eggnog_cmd, shell=True, check=True)
            
        print(f"[INFO] EggNOG annotation completed. Results available in {eggnog_output_dir}")
        
        # Extract functional annotations (KOs, COGs, CAZy)
        print(f"[INFO] Extracting functional annotations (KOs, COGs, CAZy)")
        extract_functional_annotations(eggnog_output_dir, project_output)
        
        # Run the functional analysis
        run_functional_analysis(project_output)
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] EggNOG-mapper process failed: {e}")
    except Exception as e:
        print(f"[ERROR] An error occurred during EggNOG annotation: {e}")
    finally:
        # Clean up temporary database files
        if os.path.exists(tmp_db_dir):
            print(f"[INFO] Cleaning up temporary EggNOG database files")
            shutil.rmtree(tmp_db_dir)

def extract_functional_annotations(eggnog_output_dir, project_output):
    """
    Extract KOs, COGs, and CAZy annotations from EggNOG-mapper output.
    
    Parameters:
    -----------
    eggnog_output_dir : str
        Path to directory containing EggNOG-mapper output files
    project_output : str
        Path to project output directory
    """
    # Create output directories
    functional_dir = os.path.join(project_output, "functional_annotations")
    os.makedirs(functional_dir, exist_ok=True)
    
    ko_dir = os.path.join(functional_dir, "KOs")
    cog_dir = os.path.join(functional_dir, "COGs")
    cazy_dir = os.path.join(functional_dir, "CAZy")
    
    os.makedirs(ko_dir, exist_ok=True)
    os.makedirs(cog_dir, exist_ok=True)
    os.makedirs(cazy_dir, exist_ok=True)
    
    # Count how many new annotations we're extracting
    new_ko_count = 0
    new_cog_count = 0
    new_cazy_count = 0
    
    # Process each EggNOG annotation file
    for file in os.listdir(eggnog_output_dir):
        if file.endswith(".emapper.annotations"):
            input_file = os.path.join(eggnog_output_dir, file)
            base_name = os.path.splitext(os.path.splitext(file)[0])[0]
            
            # Check if the outputs already exist
            ko_file = os.path.join(ko_dir, f"{base_name}_KOs.tsv")
            cog_file = os.path.join(cog_dir, f"{base_name}_COGs.tsv")
            cazy_file = os.path.join(cazy_dir, f"{base_name}_CAZy.tsv")
            
            # Extract KOs if needed
            if not os.path.exists(ko_file) or os.path.getsize(ko_file) == 0:
                extract_cmd = f"awk -F'\\t' '{{if ($10 ~ /^K[0-9]/) print $1\"\\t\"$10}}' {input_file} > {ko_file}"
                subprocess.run(extract_cmd, shell=True)
                new_ko_count += 1
            
            # Extract COGs if needed
            if not os.path.exists(cog_file) or os.path.getsize(cog_file) == 0:
                extract_cmd = f"awk -F'\\t' '{{if ($8 != \"-\") print $1\"\\t\"$8}}' {input_file} > {cog_file}"
                subprocess.run(extract_cmd, shell=True)
                new_cog_count += 1
            
            # Extract CAZy if needed
            if not os.path.exists(cazy_file) or os.path.getsize(cazy_file) == 0:
                extract_cmd = f"awk -F'\\t' '{{if ($21 != \"-\" && $21 != \"\") print $1\"\\t\"$21}}' {input_file} > {cazy_file}"
                subprocess.run(extract_cmd, shell=True)
                new_cazy_count += 1
    
    print(f"[INFO] Functional annotations extracted to {functional_dir}")
    print(f"[INFO] Extracted {new_ko_count} new KO files to: {ko_dir}")
    print(f"[INFO] Extracted {new_cog_count} new COG files to: {cog_dir}")
    print(f"[INFO] Extracted {new_cazy_count} new CAZy files to: {cazy_dir}")

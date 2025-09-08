import os
import shutil
import pandas as pd
from MetaMAG.utils import ensure_directory_exists

def identify_novel_mags(gtdbtk_dir, output_dir, source_mags_dir, extension=".fa"):
    """
    Identify novel MAGs from GTDB-Tk classification results and copy them to output directory.
    
    Args:
        gtdbtk_dir (str): Directory containing GTDB-Tk classification results
        output_dir (str): Directory to place novel MAGs
        source_mags_dir (str): Directory containing source MAG files
        extension (str): File extension of MAG files
        
    Returns:
        bool: Success or failure
    """
    # Ensure output directory exists
    ensure_directory_exists(output_dir)
    
    # Find GTDB-Tk summary files - check multiple possible locations
    possible_locations = [
        # Standard classify subdirectory
        os.path.join(gtdbtk_dir, "classify", "gtdbtk.bac120.summary.tsv"),
        # Directly in the GTDB-Tk directory
        os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv"),
        # Alternative naming conventions
        os.path.join(gtdbtk_dir, "summary", "gtdbtk.bac120.summary.tsv"),
    ]
    
    bac_summary = None
    for location in possible_locations:
        if os.path.exists(location):
            bac_summary = location
            print(f"[INFO] Found bacterial summary at: {bac_summary}")
            break
    
    # Similarly for archaeal summary
    ar_possible_locations = [
        os.path.join(gtdbtk_dir, "classify", "gtdbtk.ar53.summary.tsv"),
        os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv"),
        os.path.join(gtdbtk_dir, "summary", "gtdbtk.ar53.summary.tsv"),
    ]
    
    ar_summary = None
    for location in ar_possible_locations:
        if os.path.exists(location):
            ar_summary = location
            print(f"[INFO] Found archaeal summary at: {ar_summary}")
            break
    
    novel_mags = []
    
    # Process bacterial summary file
    if bac_summary:
        print(f"[INFO] Processing bacterial classification summary: {bac_summary}")
        try:
            bac_df = pd.read_csv(bac_summary, sep='\t')
            # Find MAGs with missing species designation
            novel_bac_mags = bac_df[bac_df['classification'].str.split(';').str[-1].str.strip() == 's__']['user_genome'].tolist()
            novel_mags.extend(novel_bac_mags)
            print(f"[INFO] Found {len(novel_bac_mags)} novel bacterial MAGs")
        except Exception as e:
            print(f"[ERROR] Failed to process bacterial summary: {e}")
    else:
        print(f"[WARNING] Bacterial summary file not found in any of the expected locations")
        print(f"[WARNING] Checked: {possible_locations}")
    
    # Process archaeal summary file if it exists
    if ar_summary:
        print(f"[INFO] Processing archaeal classification summary: {ar_summary}")
        try:
            ar_df = pd.read_csv(ar_summary, sep='\t')
            # Find MAGs with missing species designation
            novel_ar_mags = ar_df[ar_df['classification'].str.split(';').str[-1].str.strip() == 's__']['user_genome'].tolist()
            novel_mags.extend(novel_ar_mags)
            print(f"[INFO] Found {len(novel_ar_mags)} novel archaeal MAGs")
        except Exception as e:
            print(f"[ERROR] Failed to process archaeal summary: {e}")
    else:
        print(f"[INFO] Archaeal summary file not found (this is normal if no archaea are present)")
    
    # Exit if no novel MAGs found
    if not novel_mags:
        print("[WARNING] No novel MAGs identified")
        return False
    
    print(f"[INFO] Total novel MAGs identified: {len(novel_mags)}")
    
    # Copy novel MAGs to output directory
    copied_count = 0
    for mag in novel_mags:
        source_file = os.path.join(source_mags_dir, f"{mag}{extension}")
        dest_file = os.path.join(output_dir, f"{mag}{extension}")
        
        if os.path.exists(source_file):
            try:
                shutil.copy2(source_file, dest_file)
                copied_count += 1
            except Exception as e:
                print(f"[ERROR] Failed to copy {mag}: {e}")
        else:
            print(f"[WARNING] Source MAG file not found: {source_file}")
    
    print(f"[INFO] Successfully copied {copied_count} novel MAGs to {output_dir}")
    return True

def run(cpus, memory, time, project_input, project_output, custom_gtdbtk_dir=None):
    """
    Main function to identify novel MAGs.
    
    Args:
        cpus (int): Number of CPUs allocated
        memory (str): Memory allocation 
        time (str): Time allocation
        project_input (str): Project input directory
        project_output (str): Project output directory
        custom_gtdbtk_dir (str): Custom path to GTDB-Tk results (optional)
        
    Returns:
        bool: Success or failure
    """
    print("[INFO] Starting novel MAGs identification")
    
    # Use custom GTDB-Tk directory if provided, otherwise use default
    if custom_gtdbtk_dir:
        gtdbtk_dir = custom_gtdbtk_dir
        print(f"[INFO] Using custom GTDB-Tk directory: {gtdbtk_dir}")
    else:
        gtdbtk_dir = os.path.join(project_output, "Novel_Mags", "gtdbtk")
        print(f"[INFO] Using default GTDB-Tk directory: {gtdbtk_dir}")
    
    # Source MAGs could be from dRep output or a previous step
    source_mags_dir = os.path.join(project_output, "Bin_Refinement", "drep", "dRep_output", "dereplicated_genomes")
    
    # Output directory for novel MAGs
    output_dir = os.path.join(project_output, "Novel_Mags", "UniqueMags")
    
    # Check if GTDB-Tk results exist
    if not os.path.exists(gtdbtk_dir):
        print(f"[ERROR] GTDB-Tk directory not found: {gtdbtk_dir}")
        return False
    
    # Check if source MAGs directory exists
    if not os.path.exists(source_mags_dir):
        print(f"[ERROR] Source MAGs directory not found: {source_mags_dir}")
        return False
    
    # Process GTDB-Tk results and identify novel MAGs
    return identify_novel_mags(gtdbtk_dir, output_dir, source_mags_dir)

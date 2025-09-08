import os
import shutil
import pandas as pd
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.novel_mags import identify_novel_mags
from MetaMAG.mags_repository import add_mags_to_repo, drep_nmags_with_repo
from MetaMAG.kraken_database import add_nmags_to_kraken
from MetaMAG.config import config

def process_novel_mags_pipeline(cpus, memory, time, project_input, project_output, 
                               is_rumen=False, input_mags_dir=None, kraken_db_path=None,
                               rumen_ref_mags_dir=None, rumen_added_mags_dir=None, merge_mode=True,custom_gtdbtk_dir=None):
    """
    Comprehensive pipeline for processing novel MAGs:
    1. Identify novel MAGs using GTDB-Tk results
    2. Dereplicate against existing MAGs repository
    3. If rumen data, dereplicate against rumen reference MAGs (explicitly provided)
    4. Add to MAGs repository
    5. Build Kraken database with merge support
    
    Args:
        cpus (int): Number of CPUs to use
        memory (str): Memory allocation
        time (str): Time allocation
        project_input (str): Project input directory
        project_output (str): Project output directory
        is_rumen (bool): Whether the data is from rumen source
        input_mags_dir (str): Directory containing input MAGs (optional)
        kraken_db_path (str): Path to build Kraken database (optional)
        rumen_ref_mags_dir (str): Path to rumen reference MAGs directory (required if is_rumen=True)
        rumen_added_mags_dir (str): Path to previously added rumen MAGs directory (required if is_rumen=True)
        merge_mode (bool): Whether to merge with existing Kraken database (default: True)
    
    Returns:
        bool: Success or failure
    """
    print("[INFO] Starting comprehensive novel MAGs pipeline")
    
    # Define directories
    if not input_mags_dir:
        input_mags_dir = os.path.join(project_output, "Bin_Refinement", "drep", "dRep_output", "dereplicated_genomes")
    
   # Use custom GTDB-Tk directory if provided, otherwise use default
    if custom_gtdbtk_dir:
        gtdbtk_dir = custom_gtdbtk_dir
        print(f"[INFO] Using custom GTDB-Tk directory: {gtdbtk_dir}")
    else:
        gtdbtk_dir = os.path.join(project_output, "Novel_Mags", "gtdbtk")
        print(f"[INFO] Using default GTDB-Tk directory: {gtdbtk_dir}")


    candidate_nmags_dir = os.path.join(project_output, "Novel_Mags", "UniqueMags")
    filtered_nmags_dir = os.path.join(project_output, "Novel_Mags", "filtered_NMAGs")
    true_novel_mags_dir = os.path.join(project_output, "Novel_Mags", "true_novel_MAGs")
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    
    # Check if rumen paths are provided when is_rumen=True
    if is_rumen:
        if not rumen_ref_mags_dir:
            print("[ERROR] Rumen reference MAGs directory (rumen_ref_mags_dir) must be provided when is_rumen=True")
            return False
        if not rumen_added_mags_dir:
            print("[ERROR] Rumen added MAGs directory (rumen_added_mags_dir) must be provided when is_rumen=True")
            return False
            
        # Check if the directories exist
        if not os.path.exists(rumen_ref_mags_dir):
            print(f"[ERROR] Rumen reference MAGs directory not found: {rumen_ref_mags_dir}")
            return False
        
        # Check if the rumen reference directory has a 'mags' subdirectory
        rumen_mags_dir = os.path.join(rumen_ref_mags_dir, "mags")
        if os.path.exists(rumen_mags_dir):
            print(f"[INFO] Found mags subdirectory in rumen reference directory: {rumen_mags_dir}")
            # Use this as the actual reference directory for dereplication
            rumen_ref_fasta_dir = rumen_mags_dir
        else:
            # Use the main directory if no mags subdirectory exists
            rumen_ref_fasta_dir = rumen_ref_mags_dir
            
        # Ensure rumen_added_mags_dir exists (it might be empty for first run)
        ensure_directory_exists(rumen_added_mags_dir)
    
    # Create necessary directories
    for directory in [candidate_nmags_dir, filtered_nmags_dir, true_novel_mags_dir, mags_repo_dir]:
        ensure_directory_exists(directory)
    
    # Step 1: Run GTDB-Tk to classify MAGs (assuming GTDB-Tk has been run)
    # This would be handled by a separate step in the pipeline
    
    # Step 2: Identify novel MAGs from GTDB-Tk results
    print("[INFO] Step 1: Identifying novel MAGs from GTDB-Tk results")
    success = identify_novel_mags(gtdbtk_dir, candidate_nmags_dir, input_mags_dir)
    
    if not success:
        print("[WARNING] No novel MAGs identified. Pipeline will exit.")
        return False
    
    # Step 3: Dereplicate candidate novel MAGs against existing MAGs repository
    print("[INFO] Step 2: Dereplicating novel MAGs against existing repository")
    
    # Check if MAGs repository is empty
    repo_empty = len([f for f in os.listdir(mags_repo_dir) if f.endswith(('.fa', '.fasta'))]) == 0
    
    if repo_empty:
        # If repository is empty, copy all candidate MAGs to filtered MAGs
        print("[INFO] MAGs repository is empty, skipping dereplication")
        for mag_file in os.listdir(candidate_nmags_dir):
            if mag_file.endswith(('.fa', '.fasta')):
                shutil.copy2(
                    os.path.join(candidate_nmags_dir, mag_file),
                    os.path.join(filtered_nmags_dir, mag_file)
                )
    else:
        # Dereplicate candidate MAGs against repository
        success = drep_nmags_with_repo(project_output, cpus, candidate_nmags_dir, mags_repo_dir, filtered_nmags_dir)
        
        if not success:
            print("[WARNING] Dereplication against repository failed.")
            return False
    
    # Step 4: For rumen data, dereplicate against rumen reference MAGs
    if is_rumen:
        print("[INFO] Step 3: Dereplicating against rumen reference MAGs")
        
        # Check if rumen reference MAGs exist in the appropriate directory
        rumen_refs_exist = False
        if rumen_ref_fasta_dir and os.path.exists(rumen_ref_fasta_dir):
            rumen_refs_exist = len([f for f in os.listdir(rumen_ref_fasta_dir) if f.endswith(('.fa', '.fasta'))]) > 0
        
        rumen_added_exist = False
        if os.path.exists(rumen_added_mags_dir):
            rumen_added_exist = len([f for f in os.listdir(rumen_added_mags_dir) if f.endswith(('.fa', '.fasta'))]) > 0
        
        if not (rumen_refs_exist or rumen_added_exist):
            print("[INFO] No rumen reference or added MAGs found, skipping rumen-specific dereplication")
            # Copy all filtered MAGs to true novel MAGs
            for mag_file in os.listdir(filtered_nmags_dir):
                if mag_file.endswith(('.fa', '.fasta')):
                    shutil.copy2(
                        os.path.join(filtered_nmags_dir, mag_file),
                        os.path.join(true_novel_mags_dir, mag_file)
                    )
        else:
            # Create combined rumen reference directory
            combined_rumen_dir = os.path.join(project_output, "Novel_Mags", "combined_rumen_refs")
            ensure_directory_exists(combined_rumen_dir)
            
            # Copy rumen reference MAGs to combined directory
            if rumen_refs_exist:
                print(f"[INFO] Copying reference MAGs from: {rumen_ref_fasta_dir}")
                for mag_file in os.listdir(rumen_ref_fasta_dir):
                    if mag_file.endswith(('.fa', '.fasta')):
                        # Make sure to prefix reference MAGs for later identification
                        shutil.copy2(
                            os.path.join(rumen_ref_fasta_dir, mag_file),
                            os.path.join(combined_rumen_dir, f"ref_{mag_file}")
                        )
            
            # Copy rumen added MAGs to combined directory
            if rumen_added_exist:
                print(f"[INFO] Copying previously added MAGs from: {rumen_added_mags_dir}")
                for mag_file in os.listdir(rumen_added_mags_dir):
                    if mag_file.endswith(('.fa', '.fasta')):
                        shutil.copy2(
                            os.path.join(rumen_added_mags_dir, mag_file),
                            os.path.join(combined_rumen_dir, f"added_{mag_file}")
                        )
            
            # Dereplicate filtered MAGs against combined rumen references
            success = drep_nmags_with_repo(project_output, cpus, filtered_nmags_dir, combined_rumen_dir, true_novel_mags_dir)
            
            if not success:
                print("[WARNING] Dereplication against rumen references failed.")
                return False
    else:
        # For non-rumen data, filtered MAGs are the true novel MAGs
        print("[INFO] Non-rumen data, skipping rumen-specific dereplication")
        # Copy all filtered MAGs to true novel MAGs
        for mag_file in os.listdir(filtered_nmags_dir):
            if mag_file.endswith(('.fa', '.fasta')):
                shutil.copy2(
                    os.path.join(filtered_nmags_dir, mag_file),
                    os.path.join(true_novel_mags_dir, mag_file)
                )
    
    # Step 5: Add true novel MAGs to repository
    print("[INFO] Step 4: Adding true novel MAGs to repository")
    success = add_mags_to_repo(cpus, project_input, project_output, source_dir=true_novel_mags_dir)
    
    if not success:
        print("[WARNING] Failed to add MAGs to repository.")
        return False
    
    # If rumen data, also copy to rumen_added_mags_dir
    if is_rumen:
        print(f"[INFO] Copying novel MAGs to rumen added MAGs directory: {rumen_added_mags_dir}")
        for mag_file in os.listdir(true_novel_mags_dir):
            if mag_file.endswith(('.fa', '.fasta')):
                shutil.copy2(
                    os.path.join(true_novel_mags_dir, mag_file),
                    os.path.join(rumen_added_mags_dir, mag_file)
                )
    
    # Step 6: Build Kraken database with novel MAGs
    print("[INFO] Step 5: Building Kraken database with novel MAGs")
    # Pass the merge_mode parameter to add_nmags_to_kraken
    success = add_nmags_to_kraken(project_output, cpus, kraken_db_path, rumen_ref_mags_dir, read_length=150, merge_mode=merge_mode)
    
    if not success:
        print("[WARNING] Failed to build Kraken database.")
        return False
    
    print("[INFO] Comprehensive novel MAGs pipeline completed successfully")
    return True

def run(cpus, memory, time, project_input, project_output, 
       is_rumen=False, input_mags_dir=None, kraken_db_path=None,
       rumen_ref_mags_dir=None, rumen_added_mags_dir=None, merge_mode=True,custom_gtdbtk_dir=None):
    """
    Main run function to process novel MAGs pipeline.
    
    Args:
        cpus (int): Number of CPUs to use
        memory (str): Memory allocation
        time (str): Time allocation
        project_input (str): Project input directory
        project_output (str): Project output directory
        is_rumen (bool): Whether the data is from rumen source
        input_mags_dir (str): Directory containing input MAGs (optional)
        kraken_db_path (str): Path to build Kraken database (optional)
        rumen_ref_mags_dir (str): Path to rumen reference MAGs directory (required if is_rumen=True)
        rumen_added_mags_dir (str): Path to previously added rumen MAGs directory (required if is_rumen=True)
        merge_mode (bool): Whether to merge with existing Kraken database (default: True)
        custom_gtdbtk_dir (str): Custom path to GTDB-Tk results directory (optional)
    
    Returns:
        bool: Success or failure
    """
    return process_novel_mags_pipeline(
        cpus, memory, time, project_input, project_output, 
        is_rumen, input_mags_dir, kraken_db_path,
        rumen_ref_mags_dir, rumen_added_mags_dir, merge_mode, custom_gtdbtk_dir
    )

def repository_dereplication(cpus, memory, time, project_input, project_output, **kwargs):
    """Step 2: Dereplicate against existing MAGs repository"""
    print("[INFO] Starting repository dereplication step")
    
    candidate_nmags_dir = os.path.join(project_output, "Novel_Mags", "UniqueMags")
    filtered_nmags_dir = os.path.join(project_output, "Novel_Mags", "filtered_NMAGs")
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    
    ensure_directory_exists(filtered_nmags_dir)
    
    # Check if MAGs repository is empty
    repo_empty = len([f for f in os.listdir(mags_repo_dir) if f.endswith(('.fa', '.fasta'))]) == 0
    
    if repo_empty:
        print("[INFO] MAGs repository is empty, skipping dereplication")
        for mag_file in os.listdir(candidate_nmags_dir):
            if mag_file.endswith(('.fa', '.fasta')):
                shutil.copy2(
                    os.path.join(candidate_nmags_dir, mag_file),
                    os.path.join(filtered_nmags_dir, mag_file)
                )
    else:
        success = drep_nmags_with_repo(project_output, cpus, candidate_nmags_dir, mags_repo_dir, filtered_nmags_dir)
        if not success:
            print("[ERROR] Repository dereplication failed")
            return False
    
    print("[INFO] Repository dereplication completed successfully")
    return True

def rumen_dereplication(cpus, memory, time, project_input, project_output, 
                       rumen_ref_mags_dir=None, rumen_added_mags_dir=None, **kwargs):
    """Step 3: Dereplicate against rumen reference MAGs"""
    print("[INFO] Starting rumen-specific dereplication step")
    
    if not rumen_ref_mags_dir or not rumen_added_mags_dir:
        print("[ERROR] Rumen directories must be provided for rumen dereplication")
        return False
    
    filtered_nmags_dir = os.path.join(project_output, "Novel_Mags", "filtered_NMAGs")
    true_novel_mags_dir = os.path.join(project_output, "Novel_Mags", "true_novel_MAGs")
    
    ensure_directory_exists(true_novel_mags_dir)
    
    # Check if rumen reference MAGs exist
    rumen_mags_dir = os.path.join(rumen_ref_mags_dir, "mags")
    if os.path.exists(rumen_mags_dir):
        rumen_ref_fasta_dir = rumen_mags_dir
    else:
        rumen_ref_fasta_dir = rumen_ref_mags_dir
    
    rumen_refs_exist = len([f for f in os.listdir(rumen_ref_fasta_dir) if f.endswith(('.fa', '.fasta'))]) > 0
    rumen_added_exist = os.path.exists(rumen_added_mags_dir) and len([f for f in os.listdir(rumen_added_mags_dir) if f.endswith(('.fa', '.fasta'))]) > 0
    
    if not (rumen_refs_exist or rumen_added_exist):
        print("[INFO] No rumen reference MAGs found, copying filtered MAGs as final")
        for mag_file in os.listdir(filtered_nmags_dir):
            if mag_file.endswith(('.fa', '.fasta')):
                shutil.copy2(
                    os.path.join(filtered_nmags_dir, mag_file),
                    os.path.join(true_novel_mags_dir, mag_file)
                )
    else:
        # Create combined rumen reference directory and dereplicate
        combined_rumen_dir = os.path.join(project_output, "Novel_Mags", "combined_rumen_refs")
        ensure_directory_exists(combined_rumen_dir)
        
        # Copy reference and added MAGs with prefixes
        if rumen_refs_exist:
            for mag_file in os.listdir(rumen_ref_fasta_dir):
                if mag_file.endswith(('.fa', '.fasta')):
                    shutil.copy2(
                        os.path.join(rumen_ref_fasta_dir, mag_file),
                        os.path.join(combined_rumen_dir, f"ref_{mag_file}")
                    )
        
        if rumen_added_exist:
            for mag_file in os.listdir(rumen_added_mags_dir):
                if mag_file.endswith(('.fa', '.fasta')):
                    shutil.copy2(
                        os.path.join(rumen_added_mags_dir, mag_file),
                        os.path.join(combined_rumen_dir, f"added_{mag_file}")
                    )
        
        success = drep_nmags_with_repo(project_output, cpus, filtered_nmags_dir, combined_rumen_dir, true_novel_mags_dir)
        if not success:
            print("[ERROR] Rumen dereplication failed")
            return False
    
    print("[INFO] Rumen dereplication completed successfully")
    return True

def add_to_repository_step(cpus, memory, time, project_input, project_output, **kwargs):
    """Step 4: Add novel MAGs to repository"""
    print("[INFO] Starting add to repository step")
    
    true_novel_mags_dir = os.path.join(project_output, "Novel_Mags", "true_novel_MAGs")
    
    success = add_mags_to_repo(cpus, project_input, project_output, source_dir=true_novel_mags_dir)
    
    if not success:
        print("[ERROR] Failed to add MAGs to repository")
        return False
    
    print("[INFO] Add to repository completed successfully")
    return True

def prepare_taxonomy_input_step(cpus, memory, time, project_input, project_output, 
                               kraken_db_path=None, custom_gtdbtk_dir=None, rumen_ref_mags_dir=None, **kwargs):
    """Step 5: Prepare taxonomy input files"""
    print("[INFO] Starting taxonomy input preparation step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import prepare_taxonomy_input
    
    success = prepare_taxonomy_input(
        project_output=project_output,
        kraken_db_path=kraken_db_path,
        custom_gtdbtk_dir=custom_gtdbtk_dir,
        rumen_ref_mags_dir=rumen_ref_mags_dir
    )
    
    if not success:
        print("[ERROR] Failed to prepare taxonomy input")
        return False
    
    print("[INFO] Taxonomy input preparation completed successfully")
    return True

def create_taxonomy_step(cpus, memory, time, project_input, project_output, 
                        kraken_db_path=None, merge_mode=True, custom_gtdbtk_dir=None, **kwargs):
    """Step 6: Create taxonomy files (nodes.dmp, names.dmp)"""
    print("[INFO] Starting taxonomy creation step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import create_taxonomy_files
    
    success = create_taxonomy_files(
        project_output=project_output,
        kraken_db_path=kraken_db_path,
        merge_mode=merge_mode,
        custom_gtdbtk_dir=custom_gtdbtk_dir
    )
    
    if not success:
        print("[ERROR] Failed to create taxonomy files")
        return False
    
    print("[INFO] Taxonomy creation completed successfully")
    return True

def rename_headers_step(cpus, memory, time, project_input, project_output, 
                       kraken_db_path=None, **kwargs):
    """Step 7: Rename FASTA headers for Kraken compatibility"""
    print("[INFO] Starting header renaming step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import rename_fasta_headers
    
    success = rename_fasta_headers(
        project_output=project_output,
        kraken_db_path=kraken_db_path,
        cpus=cpus
    )
    
    if not success:
        print("[ERROR] Failed to rename headers")
        return False
    
    print("[INFO] Header renaming completed successfully")
    return True

def add_kraken_library_step(cpus, memory, time, project_input, project_output, 
                           kraken_db_path=None, **kwargs):
    """Step 8: Add MAGs to Kraken library"""
    print("[INFO] Starting add to Kraken library step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import add_to_kraken_library
    
    success = add_to_kraken_library(
        kraken_db_path=kraken_db_path,
        cpus=cpus
    )
    
    if not success:
        print("[ERROR] Failed to add to Kraken library")
        return False
    
    print("[INFO] Add to Kraken library completed successfully")
    return True

def build_kraken_database_step(cpus, memory, time, project_input, project_output, 
                              kraken_db_path=None, **kwargs):
    """Step 9: Build main Kraken database"""
    print("[INFO] Starting Kraken database building step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import build_kraken_db
    
    success = build_kraken_db(
        kraken_db_path=kraken_db_path,
        cpus=cpus
    )
    
    if not success:
        print("[ERROR] Failed to build Kraken database")
        return False
    
    print("[INFO] Kraken database building completed successfully")
    return True

def build_bracken_database_step(cpus, memory, time, project_input, project_output, 
                               kraken_db_path=None, read_length=150, **kwargs):
    """Step 10: Build Bracken database"""
    print("[INFO] Starting Bracken database building step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import build_bracken_db
    
    success = build_bracken_db(
        kraken_db_path=kraken_db_path,
        cpus=cpus,
        read_length=read_length
    )
    
    if not success:
        print("[ERROR] Failed to build Bracken database")
        return False
    
    print("[INFO] Bracken database building completed successfully")
    return True

def build_database_distribution_step(cpus, memory, time, project_input, project_output, 
                                    kraken_db_path=None, read_length=150, **kwargs):
    """Step 11: Build database distribution"""
    print("[INFO] Starting database distribution building step")
    
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    from MetaMAG.kraken_database import build_db_distribution
    
    success = build_db_distribution(
        kraken_db_path=kraken_db_path,
        cpus=cpus,
        read_length=read_length
    )
    
    if not success:
        print("[ERROR] Failed to build database distribution")
        return False
    
    print("[INFO] Database distribution building completed successfully")
    return True

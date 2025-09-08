import os
import shutil
from MetaMAG.utils import ensure_directory_exists, run_command

def add_mags_to_repo(cpus, project_input, project_output, source_dir=None):
    """
    Add MAGs to repository for future reuse.
    
    Args:
        cpus (int): Number of CPUs to use
        project_input (str): Project input directory
        project_output (str): Project output directory
        source_dir (str): Directory containing MAGs to add (default: project_output/Novel_Mags/true_novel_MAGs)
    
    Returns:
        bool: Success or failure
    """
    print("[INFO] Starting MAGs repository processing")
    
    # Set paths
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    
    if not source_dir:
        source_dir = os.path.join(project_output, "Novel_Mags", "true_novel_MAGs")
    
    # Check if source directory exists
    if not os.path.exists(source_dir):
        print(f"[ERROR] Source directory not found: {source_dir}")
        return False
    
    # Check if repository exists and has MAGs
    repo_has_mags = False
    if os.path.exists(mags_repo_dir):
        mags_count = len([f for f in os.listdir(mags_repo_dir) if f.endswith(('.fa', '.fasta'))])
        if mags_count > 0:
            repo_has_mags = True
            print(f"[INFO] MAGs Repository has {mags_count} existing MAGs")
        else:
            print(f"[INFO] MAGs Repository is empty")
    else:
        print(f"[INFO] MAGs Repository is empty")
    
    # Ensure repository directory exists
    ensure_directory_exists(mags_repo_dir)
    
    # Get list of MAGs to add
    mags_to_add = [f for f in os.listdir(source_dir) if f.endswith(('.fa', '.fasta'))]
    
    if not mags_to_add:
        print(f"[INFO] No MAGs found in source directory: {source_dir}")
        return True  # Not an error, just no MAGs to add
    
    if not repo_has_mags:
        # Repository is empty, simply copy all MAGs
        print(f"[INFO] MAGs Repository is empty, adding Novel MAGs directly")
        
        print(f"[INFO] Adding {len(mags_to_add)} Novel MAGs to repository")
        for mag_file in mags_to_add:
            shutil.copy2(
                os.path.join(source_dir, mag_file),
                os.path.join(mags_repo_dir, mag_file)
            )
        
        print(f"[INFO] Added {len(mags_to_add)} new MAGs to repository")
        return True
    else:
        # For now, simply add all MAGs (we've already done dereplication)
        print(f"[INFO] Adding {len(mags_to_add)} Novel MAGs to repository")
        for mag_file in mags_to_add:
            shutil.copy2(
                os.path.join(source_dir, mag_file),
                os.path.join(mags_repo_dir, mag_file)
            )
        
        print(f"[INFO] Added {len(mags_to_add)} new MAGs to repository")
        return True

def drep_nmags_with_repo(project_output, cpus, input_dir, reference_dir, output_dir):
    """
    Dereplicate Novel MAGs against reference MAGs.
    
    Args:
        project_output (str): Project output directory
        cpus (int): Number of CPUs to use
        input_dir (str): Directory containing novel MAGs to dereplicate
        reference_dir (str): Directory containing reference MAGs
        output_dir (str): Directory to store results (MAGs that pass dereplication)
    
    Returns:
        bool: Success or failure
    """
    # Create temporary directory for combined MAGs
    temp_dir = os.path.join(project_output, "MAGs_Repository", "temp_combined")
    ensure_directory_exists(temp_dir)
    
    # Create drep output directory
    drep_output_dir = os.path.join(project_output, "MAGs_Repository", "drep_output")
    ensure_directory_exists(drep_output_dir)
    
    # Count MAGs in both directories
    input_mags = [f for f in os.listdir(input_dir) if f.endswith(('.fa', '.fasta'))]
    ref_mags = [f for f in os.listdir(reference_dir) if f.endswith(('.fa', '.fasta'))]
    
    if not input_mags:
        print(f"[ERROR] No MAGs found in input directory: {input_dir}")
        return False
    
    if not ref_mags:
        print(f"[INFO] No reference MAGs found in: {reference_dir}")
        
        # Simply copy all input MAGs to output
        for mag_file in input_mags:
            shutil.copy2(
                os.path.join(input_dir, mag_file),
                os.path.join(output_dir, mag_file)
            )
        
        return True
    
    # Create genome list file for dRep
    genome_list_file = os.path.join(temp_dir, "genome_list.txt")
    with open(genome_list_file, 'w') as f:
        # Write input MAGs with "new_" prefix
        for mag_file in input_mags:
            mag_path = os.path.abspath(os.path.join(input_dir, mag_file))
            new_mag_path = os.path.join(temp_dir, f"new_{mag_file}")
            shutil.copy2(mag_path, new_mag_path)
            f.write(f"{new_mag_path}\n")
        
        # Write reference MAGs with "ref_" prefix if they don't already have it
        for mag_file in ref_mags:
            mag_path = os.path.abspath(os.path.join(reference_dir, mag_file))
            if mag_file.startswith('ref_'):
                ref_mag_path = os.path.join(temp_dir, mag_file)
            else:
                ref_mag_path = os.path.join(temp_dir, f"ref_{mag_file}")
            shutil.copy2(mag_path, ref_mag_path)
            f.write(f"{ref_mag_path}\n")
    
    print(f"[INFO] dRep {len(input_mags)} Novel MAGs with {len(ref_mags)} reference MAGs")
    
    # Run dRep with stderr redirected to stdout
    print(f"[INFO] Running dRep on combined MAGs")
    cmd = f"bash -c 'source /usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate drep && dRep dereplicate {drep_output_dir} -g {genome_list_file} -p {cpus} --S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95 2>&1'"
    
    # Try to run command, but we'll check for success ourselves
    run_command(cmd)
    
    # Always check if the output directory exists, regardless of command result
    derep_genomes_dir = os.path.join(drep_output_dir, "dereplicated_genomes")
    
    # Determine success based on directory existence, not command result
    if os.path.exists(derep_genomes_dir) and os.listdir(derep_genomes_dir):
        print(f"[INFO] Found dRep output directory with {len(os.listdir(derep_genomes_dir))} files: {derep_genomes_dir}")
        success = True
    else:
        print(f"[ERROR] dRep failed - output directory missing or empty")
        return False
    
    # Copy dereplicated MAGs to output directory
    # Ensure output directory exists
    ensure_directory_exists(output_dir)
    
    # Count both novel and reference MAGs that passed dereplication
    derep_novel_mags = [f for f in os.listdir(derep_genomes_dir) if f.startswith('new_') and f.endswith(('.fa', '.fasta'))]
    derep_ref_mags = [f for f in os.listdir(derep_genomes_dir) if f.startswith('ref_') and f.endswith(('.fa', '.fasta'))]
    
    # Copy all novel MAGs that passed dereplication to output directory, removing the "new_" prefix
    for mag_file in derep_novel_mags:
        orig_name = mag_file[4:]  # Remove 'new_' prefix
        shutil.copy2(
            os.path.join(derep_genomes_dir, mag_file),
            os.path.join(output_dir, orig_name)
        )
    
    # Copy all reference MAGs that passed dereplication to output directory, keeping the "ref_" prefix
    for mag_file in derep_ref_mags:
        shutil.copy2(
            os.path.join(derep_genomes_dir, mag_file),
            os.path.join(output_dir, mag_file)
        )
    
    print(f"[INFO] After dereplication, {len(derep_novel_mags)} novel MAGs and {len(derep_ref_mags)} reference MAGs will be kept")
    
    return True

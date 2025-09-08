import os
import pandas as pd
import shutil
import subprocess
import json
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def get_conda_base_path():
    """
    Dynamically detect conda installation path.
    Returns the conda base path or None if not found.
    """
    # Method 1: Try to get from existing config conda_env setting
    conda_env_setting = config.get("environment", {}).get("conda_env", "")
    if conda_env_setting and "miniconda3" in conda_env_setting:
        # Extract path from something like "source /path/to/miniconda3/bin/activate metamag"
        parts = conda_env_setting.split()
        for part in parts:
            if "miniconda3" in part or "anaconda3" in part:
                # Remove /bin/activate from the end if present
                if part.endswith("/bin/activate"):
                    return part.replace("/bin/activate", "")
                elif "/bin/" in part:
                    return part.split("/bin/")[0]
    
    # Method 2: Try conda info command
    try:
        result = subprocess.run(['conda', 'info', '--json'], 
                               capture_output=True, text=True, check=True)
        conda_info = json.loads(result.stdout)
        conda_root = conda_info.get('root_prefix')
        if conda_root:
            return conda_root
    except Exception as e:
        print(f"[DEBUG] Could not get conda info: {e}")
    
    # Method 3: Try to find from CONDA_PREFIX environment variable
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        # If we're in a conda environment, get the base
        base_path = conda_prefix
        while base_path and base_path != '/' and 'envs' in base_path:
            base_path = os.path.dirname(base_path)
        if base_path and (base_path.endswith('miniconda3') or base_path.endswith('anaconda3')):
            return base_path
    
    # Method 4: Check common conda installation locations
    common_paths = [
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/anaconda3"),
        "/opt/miniconda3",
        "/opt/anaconda3",
        "/usr/local/miniconda3",
        "/usr/local/anaconda3"
    ]
    
    for path in common_paths:
        if os.path.exists(os.path.join(path, "etc", "profile.d", "conda.sh")):
            return path
    
    print("[WARNING] Could not automatically detect conda installation path")
    return None

def run_mags_tree(cpus, memory, time, mags_dir, project_input, project_output, outgroup_taxon=None):
    """
    Runs GTDB-Tk taxonomy2 for creating phylogenetic tree from MAGs.
    This function specifically targets the tree-building functionality.

    Inputs:
    - mags_dir: Directory containing MAG files
    - outgroup_taxon: User-specified outgroup taxon (e.g., "p__Firmicutes")

    Outputs:
    - {base_dir}/Phylogeny/gtdbtk/taxonomy2/  (User-only tree)
    - {base_dir}/Phylogeny/gtdbtk/taxonomy2/trees/  (Copy of all tree files)
    """
    print(f"[INFO] Starting MAGs tree generation with GTDB-Tk...")
    
    # Base directory
    base_dir = project_output
    
    # Input directory with MAGs
    input_genomes_dir = mags_dir
    
    if not os.path.exists(input_genomes_dir) or not os.listdir(input_genomes_dir):
        print(f"[ERROR] MAGs directory is missing or empty: {input_genomes_dir}")
        return
    
    # Define GTDB-Tk output directories
    taxonomy_output = os.path.join(base_dir, "Novel_Mags", "gtdbtk")
    taxonomy2_output = os.path.join(base_dir, "Phylogeny", "gtdbtk", "taxonomy2")
    
    # Ensure directories exist
    ensure_directory_exists(taxonomy_output)
    ensure_directory_exists(taxonomy2_output)
    
    # Check for classification files (we need these to build the tree)
    bac120_file = os.path.join(taxonomy_output, "gtdbtk.bac120.summary.tsv")
    
    if not os.path.exists(bac120_file):
        print(f"[ERROR] GTDB-Tk classification file not found: {bac120_file}")
        print("[ERROR] Run the gtdbtk step first to generate classification.")
        return
    
    # Generate the custom taxonomy file for tree construction
    custom_taxonomy_file = create_custom_taxonomy_file(taxonomy_output, bac120_file)
    
    if not custom_taxonomy_file:
        print("[ERROR] Failed to create custom taxonomy file. Exiting.")
        return
    
    # Check if outgroup taxon is provided
    if not outgroup_taxon:
        print("[ERROR] No outgroup taxon provided. Please specify an outgroup taxon (e.g., p__Firmicutes).")
        return
    
    print(f"[INFO] Using user-specified outgroup taxon: {outgroup_taxon}")
    
    # Run GTDB-Tk custom tree workflow
    run_gtdbtk_custom_tree(input_genomes_dir, taxonomy2_output, cpus, custom_taxonomy_file, outgroup_taxon)
    
    # Create a trees directory and copy all tree files there
    collect_tree_files(taxonomy2_output)

def create_custom_taxonomy_file(taxonomy_output, bac120_file):
    """
    Creates the custom taxonomy file required for GTDB-Tk de_novo_wf using taxonomy information
    from the classification summary (either BAC or AR).
    The format is: GENOME_ID<TAB>d__;p__;c__;o__;f__;g__;s__
    """

    print("[INFO] Creating custom taxonomy file for GTDB-Tk de_novo_wf")

    # Use bac120_file for bacteria (if exists)
    if os.path.exists(bac120_file):
        df_bac = pd.read_csv(bac120_file, sep='\t')
        # Select the first two columns and remove the header
        df_bac = df_bac.iloc[:, :2].dropna()
        custom_taxonomy_file = os.path.join(taxonomy_output, "custom_taxonomy.txt")
        df_bac.to_csv(custom_taxonomy_file, sep='\t', header=False, index=False)
        print(f"[INFO] Custom taxonomy file created: {custom_taxonomy_file}")
        return custom_taxonomy_file
    else:
        print(f"[ERROR] {bac120_file} not found for custom taxonomy creation.")
        return None

def run_gtdbtk_custom_tree(input_genomes_dir, taxonomy2_output, cpus, custom_taxonomy_file, outgroup_taxon):
    """
    Builds a custom phylogenetic tree using only user genomes (no GTDB genomes).
    """
    print(f"[INFO] Running GTDB-Tk custom tree workflow (User genomes only)")

    # Define expected output files
    final_tree_file = os.path.join(taxonomy2_output, "gtdbtk.tree")
    
    # Skip if tree already exists
    if os.path.exists(final_tree_file):
        print(f"[INFO] GTDB-Tk tree file already exists at {final_tree_file}. Skipping tree construction.")
        return  

    # Get GTDB-Tk tool configuration from config
    gtdbtk_tool_path = get_tool_path("gtdbtk")
    
    if not gtdbtk_tool_path:
        print(f"[ERROR] GTDB-Tk tool not found in configuration")
        return
    
    print(f"[INFO] Running GTDB-Tk de_novo_wf with outgroup: {outgroup_taxon}")
    
    # Check if the tool path is a conda activation command
    if "activate" in gtdbtk_tool_path and "gtdbtk" in gtdbtk_tool_path:
        # Extract environment name from the activation path
        env_name = gtdbtk_tool_path.split()[-1]  # Get the last part (environment name)
        
        # Dynamically get conda base path
        conda_base = get_conda_base_path()
        
        if not conda_base:
            print("[ERROR] Could not detect conda installation path. Please ensure conda is properly installed.")
            return
        
        print(f"[DEBUG] Using conda base path: {conda_base}")
        print(f"[DEBUG] Activating environment: {env_name}")
        
        # Create command that properly activates environment and runs gtdbtk
        cmd_tree = f"""bash -c "source {conda_base}/etc/profile.d/conda.sh && conda activate {env_name} && gtdbtk de_novo_wf --genome_dir {input_genomes_dir} --bacteria --outgroup_taxon '{outgroup_taxon}' --extension .fa --skip_gtdb_refs --custom_taxonomy_file {custom_taxonomy_file} --out_dir {taxonomy2_output} --cpus {cpus}" """
    else:
        # Assume it's a direct path to gtdbtk executable
        cmd_tree = f"{gtdbtk_tool_path} de_novo_wf --genome_dir {input_genomes_dir} --bacteria --outgroup_taxon '{outgroup_taxon}' --extension .fa --skip_gtdb_refs --custom_taxonomy_file {custom_taxonomy_file} --out_dir {taxonomy2_output} --cpus {cpus}"
    
    print(f"[DEBUG] Running command: {cmd_tree}")
    run_command(cmd_tree)

    print(f"[INFO] GTDB-Tk tree construction completed. Output in {taxonomy2_output}")

def collect_tree_files(taxonomy2_output):
    """
    Collects all tree files from the GTDB-Tk output directory and copies them
    to a dedicated 'trees' directory.
    """
    print("[INFO] Collecting tree files and copying to a centralized trees directory...")
    
    # Create a 'trees' directory
    trees_dir = os.path.join(taxonomy2_output, "trees")
    ensure_directory_exists(trees_dir)
    
    # Common tree file extensions
    tree_extensions = [".tree", ".nwk", ".newick", ".treefile"]
    
    # Count tree files
    tree_count = 0
    
    # Walk through the taxonomy2_output directory and find all tree files
    for root, dirs, files in os.walk(taxonomy2_output):
        for file in files:
            # Check if the file ends with any tree extension
            if any(file.endswith(ext) for ext in tree_extensions):
                source_path = os.path.join(root, file)
                dest_path = os.path.join(trees_dir, file)
                
                # Copy the tree file to the trees directory
                try:
                    shutil.copy2(source_path, dest_path)
                    print(f"[INFO] Copied tree file: {file}")
                    tree_count += 1
                except Exception as e:
                    print(f"[ERROR] Failed to copy {file}: {e}")
    
    # Also check for any .tree files that don't have these extensions
    for root, dirs, files in os.walk(taxonomy2_output):
        for file in files:
            # Check if 'tree' is in the filename but not already copied
            if "tree" in file.lower() and not any(file.endswith(ext) for ext in tree_extensions):
                source_path = os.path.join(root, file)
                dest_path = os.path.join(trees_dir, file)
                
                # Skip if the source and destination are the same
                if os.path.samefile(os.path.dirname(source_path), trees_dir):
                    continue
                
                # Copy the tree file to the trees directory
                try:
                    shutil.copy2(source_path, dest_path)
                    print(f"[INFO] Copied tree file: {file}")
                    tree_count += 1
                except Exception as e:
                    print(f"[ERROR] Failed to copy {file}: {e}")
    
    if tree_count > 0:
        print(f"[INFO] Successfully copied {tree_count} tree files to {trees_dir}")
    else:
        print("[WARNING] No tree files found to copy")
import os
import shutil
import pandas as pd
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def prepare_drep_input(samples, project_output):
    """
    Collects bins from multiple samples and co-assembly into a single `dRep_input` directory.

     Copies **single-sample** bins from:
      {project_output}/Bin_Refinement/dastool/Single_Sample/{sample}/{sample}_DAS_DASTool_bins
     Copies **co-assembly** bins (if they exist) from:
      {project_output}/Bin_Refinement/dastool/Coassembly/Coassembly_DAS_DASTool_bins
     Renames bins as:
      - `{sample}_{bin_name}.fa` (for single-sample)
      - `coassembly_{bin_name}.fa` (for co-assembly)
    
    Returns:
      - `drep_input_dir`: Path to the directory containing all renamed bins.
      - `copied_files`: List of copied file names.
    """
    base_dir = project_output
    drep_input_dir = os.path.join(base_dir, "Bin_Refinement", "drep", "dRep_input")
    ensure_directory_exists(drep_input_dir)
    copied_files = []

    # Process single-sample bins
    for sample in samples:
        das_tool_bins_dir = os.path.join(base_dir, "Bin_Refinement", "dastool", "Single_Sample", sample, f"{sample}_DAS_DASTool_bins")

        if os.path.exists(das_tool_bins_dir) and os.listdir(das_tool_bins_dir):
            for filename in os.listdir(das_tool_bins_dir):
                if filename.endswith(".fa"):
                    old_path = os.path.join(das_tool_bins_dir, filename)
                    new_filename = f"{sample}_{filename}"  # Rename bin file
                    new_path = os.path.join(drep_input_dir, new_filename)
                    shutil.copy(old_path, new_path)
                    copied_files.append(new_filename)
        else:
            print(f"[WARNING] No bins found for sample {sample} at {das_tool_bins_dir}")

    # Process co-assembly bins (if available)
    coassembly_bins_dir = os.path.join(base_dir, "Bin_Refinement", "dastool", "Coassembly", "Coassembly_DAS_DASTool_bins")

    if os.path.exists(coassembly_bins_dir) and os.listdir(coassembly_bins_dir):
        for filename in os.listdir(coassembly_bins_dir):
            if filename.endswith(".fa"):
                old_path = os.path.join(coassembly_bins_dir, filename)
                new_filename = f"coassembly_{filename}"  # Rename bin file
                new_path = os.path.join(drep_input_dir, new_filename)
                shutil.copy(old_path, new_path)
                copied_files.append(new_filename)
    else:
        print(f"[INFO] No co-assembly bins found at {coassembly_bins_dir}")

    print(f"[INFO] Copied {len(copied_files)} bin files to {drep_input_dir}")
    return drep_input_dir, copied_files

def handle_drep_failure(drep_output_dir, drep_input_dir, comp_threshold=80, con_threshold=10):
    """
    Handles the case when dRep fails due to too few genomes passing quality thresholds.
    
    This function:
    1. Checks if genomeInfo.csv exists
    2. Identifies genomes that pass quality thresholds
    3. Creates a 'dereplicated_genomes' directory
    4. Copies passing genomes to that directory
    
    Args:
        drep_output_dir: Path to dRep output directory
        drep_input_dir: Path to dRep input directory
        comp_threshold: Completeness threshold (default: 80)
        con_threshold: Contamination threshold (default: 10)
        
    Returns:
        bool: True if recovery was successful, False otherwise
    """
    genome_info_path = os.path.join(drep_output_dir, "data_tables", "genomeInfo.csv")
    
    if not os.path.exists(genome_info_path):
        print(f"[ERROR] genomeInfo.csv not found at {genome_info_path}")
        print("[INFO] This suggests the dRep run failed before quality assessment. Cannot recover.")
        return False
    
    try:
        # Read genome info
        genome_info = pd.read_csv(genome_info_path)
        
        # Apply thresholds
        passing_genomes = genome_info[
            (genome_info['completeness'] >= comp_threshold) & 
            (genome_info['contamination'] <= con_threshold)
        ]
        
        total_genomes = len(genome_info)
        num_passing = len(passing_genomes)
        
        print(f"[INFO] Found {total_genomes} total genomes in genomeInfo.csv")
        print(f"[INFO] {num_passing} genomes ({(num_passing/total_genomes)*100:.2f}%) passed quality thresholds:")
        print(f"       Completeness >= {comp_threshold}%, Contamination <= {con_threshold}%")
        
        if num_passing == 0:
            print("\n[WARNING] No genomes passed quality thresholds!")
            print("You may need to adjust thresholds or check the quality of your input genomes.")
            return False
        
        # Display passing genomes
        print("\nPassing genomes:")
        for idx, row in passing_genomes.iterrows():
            genome_name = row['genome']
            comp = row['completeness']
            con = row['contamination']
            print(f"  - {genome_name} (Completeness: {comp:.2f}%, Contamination: {con:.2f}%)")
        
        # Create dereplicated_genomes directory
        derep_genomes_dir = os.path.join(drep_output_dir, "dereplicated_genomes")
        ensure_directory_exists(derep_genomes_dir)
        
        # Copy passing genomes to the output folder
        copy_count = 0
        for genome_name in passing_genomes['genome']:
            genome_path = os.path.join(drep_input_dir, genome_name)
            if os.path.exists(genome_path):
                shutil.copy(genome_path, derep_genomes_dir)
                copy_count += 1
            else:
                print(f"[WARNING] Could not find genome file: {genome_path}")
        
        print(f"\n[SUCCESS] Copied {copy_count} out of {num_passing} passing genomes to:")
        print(f"          {derep_genomes_dir}")
        
        # Write the passing genomes list to a file
        passing_genomes_output = os.path.join(drep_output_dir, "passing_genomes.csv")
        passing_genomes.to_csv(passing_genomes_output, index=False)
        print(f"\n[INFO] Wrote passing genomes info to: {passing_genomes_output}")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] Failed to process genome information: {e}")
        return False

def run_drep(samples, cpus, memory, time, project_input, project_output):
    """
    Runs `dRep` on all collected bins from multiple samples and co-assembly.

     Input: `{project_output}/Bin_Refinement/drep/dRep_input/`
     Output: `{project_output}/Bin_Refinement/drep/dRep_output/`
    """
    base_dir = project_output
    drep_input_dir, copied_files = prepare_drep_input(samples, project_output)

    if not copied_files:
        print("[ERROR] No bins available for dRep dereplication. Exiting.")
        return
    
    # Define output directory
    drep_output_dir = os.path.join(base_dir, "Bin_Refinement", "drep", "dRep_output")
    ensure_directory_exists(drep_output_dir)

    # Use glob pattern to select all .fa files
    input_pattern = os.path.join(drep_input_dir, "*.fa")

    # First, let's check if we have the drep_activate path in the config (for backward compatibility)
    drep_activate = config.get('tools', {}).get('drep_activate')
    
    if drep_activate:
        # Use the old style activation if available
        cmd = (
            f"bash -c 'source {drep_activate} && "
            f"dRep dereplicate {drep_output_dir} -g {input_pattern} -p {cpus} "
            f"--S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95'"    
        )
    else:
        # Get the drep tool path from config
        drep_tool = get_tool_path("drep")
        if not drep_tool:
            print(f"[ERROR] dRep tool not found in configuration")
            return
        
        # Check if the drep tool path looks like an activation script
        if "activate" in drep_tool:
            # This is actually an activation script, use it properly
            cmd = (
                f"bash -c 'source {drep_tool} && "
                f"dRep dereplicate {drep_output_dir} -g {input_pattern} -p {cpus} "
                f"--S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95'"
            )
        else:
            # For dRep to work, we need CheckM and prodigal in PATH
            # Let's build a PATH that includes all necessary tools
            tool_paths = []
            
            # Add paths for required tools
            for tool in ["prodigal", "checkm"]:
                tool_path = get_tool_path(tool)
                if tool_path and os.path.exists(tool_path):
                    tool_dir = os.path.dirname(tool_path)
                    if tool_dir not in tool_paths:
                        tool_paths.append(tool_dir)
            
            # Add dRep's own directory
            drep_dir = os.path.dirname(drep_tool)
            if drep_dir not in tool_paths:
                tool_paths.append(drep_dir)
            
            # Build PATH export
            if tool_paths:
                path_export = "export PATH=" + ":".join(tool_paths) + ":$PATH && "
            else:
                path_export = ""
            
            cmd = (
                f"bash -c '{path_export}"
                f"{drep_tool} dereplicate {drep_output_dir} -g {input_pattern} -p {cpus} "
                f"--S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95'"
            )

    print(f"[INFO] Running dRep with:")
    print(f"       Input: {input_pattern}")
    print(f"       Output: {drep_output_dir}")
    print(f"       CPUs: {cpus}")
    print(f"[DEBUG] Command: {cmd}")

    try:
        run_command(cmd)
        
        # Check if dRep completed successfully
        derep_genomes_dir = os.path.join(drep_output_dir, "dereplicated_genomes")
        if os.path.exists(derep_genomes_dir) and os.listdir(derep_genomes_dir):
            print(f"[INFO] dRep dereplication completed successfully!")
            print(f"[INFO] Results saved in {drep_output_dir}")
            print(f"[INFO] Dereplicated genomes: {len(os.listdir(derep_genomes_dir))} files")
        else:
            print(f"[WARNING] dRep completed but no dereplicated genomes found.")
            print("[INFO] This might be due to strict quality thresholds.")
            
    except Exception as e:
        print(f"[WARNING] dRep command failed with error: {e}")
        print("[INFO] Attempting to recover and process passing genomes...")
        
        if handle_drep_failure(drep_output_dir, drep_input_dir):
            print("[INFO] Recovery successful. Genomes that passed quality thresholds were copied to 'dereplicated_genomes'")
        else:
            print("[ERROR] Recovery failed. Please check the dRep output directory for more information.")
            print(f"[INFO] Check logs in: {drep_output_dir}")

def run(samples, cpus, memory, time, project_input, project_output):
    """
    Main function to run dRep dereplication.
    This is called by the pipeline step runner.
    """
    print(f"[INFO] Starting dRep dereplication for {len(samples)} samples")
    run_drep(samples, cpus, memory, time, project_input, project_output)
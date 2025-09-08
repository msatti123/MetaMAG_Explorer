import os
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run_metaquast(sample, cpus, memory, time, project_input, project_output, assembly_type="both"):
    """
    Runs metaQUAST on assembly outputs for quality assessment.
    
    Parameters:
    -----------
    sample : str
        Sample ID to process
    cpus : int
        Number of CPUs to use
    memory : str
        Memory allocation
    time : str
        Time allocation
    project_input : str
        Path to project input directory
    project_output : str
        Path to project output directory
    assembly_type : str
        Type of assembly to evaluate: "idba", "megahit", or "both" (default)
    """
    print(f"Running metaQUAST for sample: {sample}")
    
    # Define output directory
    quast_output_dir = f"{project_output}/Assembly_QA/QUAST"
    ensure_directory_exists(quast_output_dir)
    
    # Get the metaQUAST binary from config
    metaquast_bin = get_tool_path("metaquast")
    if not metaquast_bin:
        print(f"[ERROR] Metaquast tool not found in configuration")
        return
    if not metaquast_bin:
        print("[ERROR] metaQUAST binary not found in config.py. Please add it to the config file.")
        return
    
    # Process IDBA assembly if requested
    if assembly_type in ["idba", "both"]:
        idba_scaffold_file = f"{project_output}/Assembly/IDBA/{sample}_idba/scaffold.fa"
        
        if os.path.exists(idba_scaffold_file):
            idba_quast_output = f"{quast_output_dir}/{sample}_idba"
            ensure_directory_exists(idba_quast_output)
            
            # Check if QUAST has already been run
            if os.path.exists(f"{idba_quast_output}/report.html"):
                print(f"[INFO] QUAST already completed for IDBA assembly of sample {sample}. Skipping.")
            else:
                # Run metaQUAST on IDBA assembly
                quast_cmd = (
                    f"{metaquast_bin} {idba_scaffold_file} -o {idba_quast_output} "
                    f"--max-ref-number 0 --threads {cpus} --silent"
                )
                try:
                    run_command(quast_cmd, log_prefix=f"QUAST (IDBA {sample})")
                    print(f"[INFO] metaQUAST completed for IDBA assembly of sample {sample}.")
                except Exception as e:
                    print(f"[ERROR] metaQUAST failed for IDBA assembly of sample {sample}: {e}")
        else:
            print(f"[WARNING] IDBA scaffold file not found for sample {sample}: {idba_scaffold_file}")

    # Process MEGAHIT co-assembly if requested
    if assembly_type in ["megahit", "both"]:
        megahit_contigs_file = f"{project_output}/Assembly/MEGAHIT/coassembly_contigs.fa"
        
        if os.path.exists(megahit_contigs_file):
            megahit_quast_output = f"{quast_output_dir}/megahit_coassembly"
            ensure_directory_exists(megahit_quast_output)
            
            # Check if QUAST has already been run
            if os.path.exists(f"{megahit_quast_output}/report.html"):
                print(f"[INFO] QUAST already completed for MEGAHIT co-assembly. Skipping.")
            else:
                # Run metaQUAST on MEGAHIT co-assembly
                quast_cmd = (
                    f"{metaquast_bin} {megahit_contigs_file} -o {megahit_quast_output} "
                    f"--max-ref-number 0 --threads {cpus} --silent"
                )
                try:
                    run_command(quast_cmd, log_prefix="QUAST (MEGAHIT)")
                    print("[INFO] metaQUAST completed for MEGAHIT co-assembly.")
                except Exception as e:
                    print(f"[ERROR] metaQUAST failed for MEGAHIT co-assembly: {e}")
        else:
            print(f"[WARNING] MEGAHIT contigs file not found: {megahit_contigs_file}")

def run_metaquast_coassembly(cpus, memory, time, project_input, project_output):
    """
    Runs metaQUAST specifically on co-assembly outputs.
    This function is designed to be called for steps that don't require sample IDs.
    """
    print("Running metaQUAST for co-assembly")
    
    # Define output directory
    quast_output_dir = f"{project_output}/Assembly_QA/QUAST"
    ensure_directory_exists(quast_output_dir)
    
    # Get the metaQUAST binary from config
    metaquast_bin = get_tool_path("metaquast")
    if not metaquast_bin:
        print("[ERROR] metaQUAST binary not found in config.py. Please add it to the config file.")
        return
    
    # Process MEGAHIT co-assembly
    megahit_contigs_file = f"{project_output}/Assembly/MEGAHIT/coassembly.contigs.fa"
    
    if os.path.exists(megahit_contigs_file):
        megahit_quast_output = f"{quast_output_dir}/megahit_coassembly"
        ensure_directory_exists(megahit_quast_output)
        
        # Check if QUAST has already been run
        if os.path.exists(f"{megahit_quast_output}/report.html"):
            print(f"[INFO] QUAST already completed for MEGAHIT co-assembly. Skipping.")
        else:
            # Run metaQUAST on MEGAHIT co-assembly
            quast_cmd = (
                f"{metaquast_bin} {megahit_contigs_file} -o {megahit_quast_output} "
                f"--max-ref-number 0 --threads {cpus} --silent"
            )
            try:
                run_command(quast_cmd, log_prefix="QUAST (MEGAHIT)")
                print("[INFO] metaQUAST completed for MEGAHIT co-assembly.")
            except Exception as e:
                print(f"[ERROR] metaQUAST failed for MEGAHIT co-assembly: {e}")
    else:
        print(f"[WARNING] MEGAHIT contigs file not found: {megahit_contigs_file}")

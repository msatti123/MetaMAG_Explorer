import os
import shutil
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run_idba(sample, cpus, memory, time, project_input, project_output):
    """
    Runs IDBA Assembly for a single sample.
    """
    print(f"Running IDBA Assembly for sample: {sample}")

    # Define output directory
    output_dir = f"{project_output}/Assembly/IDBA"
    ensure_directory_exists(output_dir)
    idba_dir = f"{output_dir}/{sample}_idba"

    # Ensure sample-specific subdirectory exists
    ensure_directory_exists(idba_dir)

    # Input files
    host_removal_dir = f"{project_output}/Host_Removal"
    unmapped_r1 = f"{host_removal_dir}/{sample}_unmapped_R1.fastq.gz"
    unmapped_r2 = f"{host_removal_dir}/{sample}_unmapped_R2.fastq.gz"

    # Validate input files
    if not os.path.exists(unmapped_r1) or not os.path.exists(unmapped_r2):
        print(f"[ERROR] Unmapped files not found for sample {sample}: {unmapped_r1}, {unmapped_r2}")
        return

    # Check if IDBA assembly is already completed
    scaffold_file = f"{idba_dir}/scaffold.fa"
    if os.path.exists(scaffold_file):
        print(f"[INFO] IDBA assembly already completed for sample {sample}. Skipping.")
        return

    # Commands for IDBA assembly
    idba_fq2fa = get_tool_path("idba_fq2fa")
    if not idba_fq2fa:
        print(f"[ERROR] Idba_Fq2Fa tool not found in configuration")
        return
    idba = get_tool_path("idba")
    if not idba:
        print(f"[ERROR] Idba tool not found in configuration")
        return
    command = (
        f"zcat {unmapped_r1} > {idba_dir}/{sample}_R1.fastq && "
        f"zcat {unmapped_r2} > {idba_dir}/{sample}_R2.fastq && "
        f"{idba_fq2fa} --merge {idba_dir}/{sample}_R1.fastq {idba_dir}/{sample}_R2.fastq {idba_dir}/{sample}.fas && "
        f"{idba} -r {idba_dir}/{sample}.fas --num_threads {cpus} -o {idba_dir}"
    )

    # Run the IDBA command
    try:
        run_command(command)
        print(f"[INFO] IDBA assembly completed for sample {sample}. Results saved in {idba_dir}.")
    except Exception as e:
        print(f"[ERROR] IDBA assembly failed for sample {sample}: {e}")

def run_megahit(samples, cpus, memory, time, project_input, project_output):
    """
    Runs MEGAHIT Co-Assembly for multiple samples.
    """
    print(f"Running MEGAHIT Co-Assembly for samples: {', '.join(samples)}")

    # Define output directory (one shared directory for all samples)
    output_dir = f"{project_output}/Assembly/MEGAHIT"
    ensure_directory_exists(output_dir)

    # Input files (ensure correct handling of sample names)
    host_removal_dir = f"{project_output}/Host_Removal"
    
    # Ensure correct formatting of sample names
    unmapped_r1_files = [f"{host_removal_dir}/{sample.strip()}_unmapped_R1.fastq.gz" for sample in samples]
    unmapped_r2_files = [f"{host_removal_dir}/{sample.strip()}_unmapped_R2.fastq.gz" for sample in samples]

    # Validate input files
    missing_files = [f for f in unmapped_r1_files + unmapped_r2_files if not os.path.exists(f)]
    if missing_files:
        print(f"[ERROR] Some unmapped reads are missing: {missing_files}")
        return

    # **?? Check if MEGAHIT output directory exists**
    contigs_file = f"{output_dir}/coassembly_contigs.fa"
    if os.path.exists(contigs_file):
        print(f"[INFO] MEGAHIT Co-Assembly already completed. Skipping.")
        return

    if os.path.exists(output_dir):
        print(f"[WARNING] MEGAHIT output directory already exists: {output_dir}")
        
        # Option 1: Remove existing directory
        shutil.rmtree(output_dir)
        print(f"[INFO] Removed existing MEGAHIT output directory: {output_dir}")

        # Option 2: Rename existing directory (instead of deleting)
        # backup_dir = f"{output_dir}_backup_{int(time.time())}"
        # shutil.move(output_dir, backup_dir)
        # print(f"[INFO] Renamed existing MEGAHIT directory to: {backup_dir}")

    # Debugging: Print sample names used in MEGAHIT
    print(f"[DEBUG] Running MEGAHIT on samples: {', '.join(samples)}")

    # Format input files for MEGAHIT
    r1_string = ",".join(unmapped_r1_files)
    r2_string = ",".join(unmapped_r2_files)

    # MEGAHIT command
    megahit = get_tool_path("megahit")
    if not megahit:
        print(f"[ERROR] Megahit tool not found in configuration")
        return
    command = (
        f"{megahit} -1 {r1_string} "
        f"-2 {r2_string} --k-min 21 --k-max 141 --k-step 12 "
        f"-t {cpus} --out-dir {output_dir} --out-prefix coassembly"
    )

    # Run the MEGAHIT command
    try:
        run_command(command, log_prefix=f"MEGAHIT ({', '.join(samples)})")
        print(f"[INFO] MEGAHIT Co-Assembly completed. Results saved in {output_dir}.")
    except Exception as e:
        print(f"[ERROR] MEGAHIT failed: {e}")
        print(f"[INFO] Ensure the directory {output_dir} is cleaned up if retrying.")

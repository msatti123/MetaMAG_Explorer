import os
import argparse
import subprocess
from MetaMAG.utils import ensure_directory_exists
from MetaMAG.config import config, get_tool_path, get_tool_command

def run(sample, cpus, memory, time, reference, project_input, project_output):
    """
    Performs host removal using BWA-MEM and Samtools with better error handling.
    """
    print(f"Running Host Removal for sample: {sample} using reference: {reference}")

    # Define output directory and files
    output_dir = f"{project_output}/Host_Removal"
    ensure_directory_exists(output_dir)
    unmapped_r1 = f"{output_dir}/{sample}_unmapped_R1.fastq.gz"
    unmapped_r2 = f"{output_dir}/{sample}_unmapped_R2.fastq.gz"

    # Check if host removal is already done
    if os.path.exists(unmapped_r1) and os.path.exists(unmapped_r2):
        print(f"[INFO] Host removal already completed for sample {sample}. Skipping.")
        return

    # Input files from the trimming step
    trimming_dir = f"{project_output}/Trimming"
    trimmed_r1 = f"{trimming_dir}/{sample}_trimmed_1.fq.gz"
    trimmed_r2 = f"{trimming_dir}/{sample}_trimmed_2.fq.gz"

    # Validate input files
    if not os.path.exists(trimmed_r1) or not os.path.exists(trimmed_r2):
        print(f"[ERROR] Trimmed files not found for sample {sample}: {trimmed_r1}, {trimmed_r2}")
        return

    # Commands for Host Removal
    bwa = get_tool_path("bwa")
    if not bwa:
        print(f"[ERROR] Bwa tool not found in configuration")
        return
    samtools = get_tool_path("samtools")
    if not samtools:
        print(f"[ERROR] Samtools tool not found in configuration")
        return
    host_removed_bam = f"{output_dir}/{sample}_host_removed.bam"
    sorted_bam = f"{output_dir}/{sample}_sorted.bam"
    temp_sam = f"{output_dir}/{sample}_host_removed.sam"  # Temporary file for debugging

    def run_command(command):
        """ Run a shell command and check for errors. """
        print(f"[INFO] Running command: {command}")
        result = subprocess.run(command, shell=True)
        if result.returncode != 0:
            print(f"[ERROR] Command failed: {command}")
            exit(1)

    # Step 1: Run BWA-MEM (Write to SAM instead of piping)
    run_command(f"{bwa} mem {reference} -t {cpus} {trimmed_r1} {trimmed_r2} > {temp_sam}")

    # Step 2: Convert SAM to BAM
    run_command(f"{samtools} view -bS {temp_sam} > {host_removed_bam}")

    # Step 3: Sort BAM (Use lower memory to avoid crashes)
    run_command(f"{samtools} sort -m 2G -o {sorted_bam} {host_removed_bam}")

    # Step 4: Extract Unmapped Reads
    run_command(f"{samtools} fastq {sorted_bam} -1 {unmapped_r1} -2 {unmapped_r2} -0 /dev/null -s /dev/null -n")

    # Cleanup temporary files
    os.remove(temp_sam)
    print(f"[INFO] Host removal completed for sample {sample}.")

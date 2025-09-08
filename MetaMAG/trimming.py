import os
import glob
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run(sample, cpus, memory, time, project_input, project_output):
    print(f"Running Trimming for sample: {sample}")

    # Define input/output directories
    input_dir = project_input
    output_dir = os.path.join(project_output, "Trimming")
    ensure_directory_exists(output_dir)

    # Detect input files (support various naming conventions)
    possible_r1_names = [f"{sample}_1", f"{sample}_R1", f"{sample}_forward"]
    possible_r2_names = [f"{sample}_2", f"{sample}_R2", f"{sample}_reverse"]

    input_r1 = None
    input_r2 = None

    for ext in [".fastq.gz", ".fastq"]:
        for r1_base, r2_base in zip(possible_r1_names, possible_r2_names):
            r1_path = os.path.join(input_dir, f"{r1_base}{ext}")
            r2_path = os.path.join(input_dir, f"{r2_base}{ext}")
            if os.path.exists(r1_path) and os.path.exists(r2_path):
                input_r1, input_r2 = r1_path, r2_path
                break
        if input_r1 and input_r2:
            break  # Stop searching if valid files are found

    # Validate input files
    if not input_r1 or not input_r2:
        print(f"[ERROR] Paired files not found for sample {sample} in {input_dir}")
        return

    # Define output files
    trimmed_r1 = os.path.join(output_dir, f"{sample}_trimmed_1.fq.gz")
    trimmed_r2 = os.path.join(output_dir, f"{sample}_trimmed_2.fq.gz")
    json_report = os.path.join(output_dir, f"{sample}_fastp.json")
    html_report = os.path.join(output_dir, f"{sample}_fastp.html")

    # Check if trimming is already completed
    if os.path.exists(trimmed_r1) and os.path.exists(trimmed_r2):
        print(f"[INFO] Trimming already completed for sample {sample}. Skipping.")
        return

    print(f"[DEBUG] Trimming output files NOT found, running trimming.")

    # Run trimming command
    fastp = get_tool_path("fastp")
    if not fastp:
        print(f"[ERROR] Fastp tool not found in configuration")
        return
    command = f"{fastp} -i {input_r1} -I {input_r2} " \
              f"-o {trimmed_r1} -O {trimmed_r2} " \
              f"-w {cpus} --detect_adapter_for_pe --cut_mean_quality 20 " \
              f"--length_required 30 -j {json_report} -h {html_report}"

    run_command(command)

    # Confirm output files were generated
    if os.path.exists(trimmed_r1) and os.path.exists(trimmed_r2):
        print(f"[INFO] Trimming completed successfully for sample {sample}.")
    else:
        print(f"[WARNING] Trimming step executed but expected output files are missing for sample {sample}.")

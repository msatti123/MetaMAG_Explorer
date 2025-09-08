import os
import glob
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

# Get PERL5LIB from config
# Set up PERL5LIB - prioritize conda environment over system
perl_lib_path = config.get("environment", {}).get("PERL5LIB", "")

# If no PERL5LIB in config, try to detect conda environment perl path
if not perl_lib_path or not os.path.exists(perl_lib_path):
    print(f"[WARNING] Configured PERL5LIB not found: {perl_lib_path}")
    
    # Try to detect conda environment PERL5LIB
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', '')
    if conda_env:
        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if conda_prefix:
            potential_perl_paths = [
                os.path.join(conda_prefix, "lib", "perl5"),
                os.path.join(conda_prefix, "lib", "site_perl"),
                os.path.join(conda_prefix, "lib", f"perl5/site_perl")
            ]
            
            for perl_path in potential_perl_paths:
                if os.path.exists(perl_path):
                    perl_lib_path = perl_path
                    print(f"[INFO] Using detected conda PERL5LIB: {perl_lib_path}")
                    break

# Set PERL5LIB environment variable
if perl_lib_path and os.path.exists(perl_lib_path):
    print(f"[INFO] Setting PERL5LIB to: {perl_lib_path}")
    perl_env = f"export PERL5LIB={perl_lib_path}; "
else:
    print(f"[WARNING] No valid PERL5LIB path found. FastQC may fail if it requires Perl modules.")
    perl_env = ""

def run(sample, cpus, memory, time, project_input, project_output):
    """
    Runs FastQC on the provided sample.

    Args:
        sample (str): Sample ID
        cpus (int): Number of CPU threads to use
        memory (str): Memory allocation
        time (str): Time limit for the step
        project_input (str): Path to the input directory from project_config.yaml
        project_output (str): Path to the output directory from project_config.yaml
    """

    # Ensure `sample` is a string (extract if list)
    if isinstance(sample, list):
        sample = sample[0]  

    print(f"[INFO] Running QC for sample: {sample}")

    # Define paths
    input_dir = project_input
    output_dir = os.path.join(project_output, "QC")

    # Find input files (handle both .fastq and .fastq.gz)
    input_r1, input_r2 = None, None
    for ext in [".fastq.gz", ".fastq"]:
        r1_path = os.path.join(input_dir, f"{sample}_1{ext}")
        r2_path = os.path.join(input_dir, f"{sample}_2{ext}")
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            input_r1, input_r2 = r1_path, r2_path
            break  # Stop searching once valid files are found

    # Validate input files
    if not input_r1 or not input_r2:
        print(f"[ERROR] Paired files not found for sample {sample} in {input_dir}. Skipping QC.")
        return

    # Define expected output filenames
    qc_html_1 = os.path.join(output_dir, f"{sample}_1_fastqc.html")
    qc_zip_1 = os.path.join(output_dir, f"{sample}_1_fastqc.zip")
    qc_html_2 = os.path.join(output_dir, f"{sample}_2_fastqc.html")
    qc_zip_2 = os.path.join(output_dir, f"{sample}_2_fastqc.zip")

    # Check if QC is already completed
    if all(os.path.exists(f) for f in [qc_html_1, qc_zip_1, qc_html_2, qc_zip_2]):
        print(f"[INFO] QC already completed for sample {sample}. Skipping.")
        return

    print(f"[DEBUG] QC output files NOT found, running QC.")

    # Ensure the output directory exists before running the command
    ensure_directory_exists(output_dir)

    # Run QC command
    fastqc = get_tool_path("fastqc")
    if not fastqc:
        print(f"[ERROR] Fastqc tool not found in configuration")
        return

    # Only export PERL5LIB if it exists
    # Construct and run the FastQC command
    command = (
        f"bash -c '{perl_env} {fastqc} -f fastq -t {cpus} -o {output_dir} {input_r1} {input_r2}'"
    )

    print(f"[DEBUG] Running command: {command}")
    
    run_command(command)

    # Confirm output files were generated
    if all(os.path.exists(f) for f in [qc_html_1, qc_zip_1, qc_html_2, qc_zip_2]):
        print(f"[INFO] QC completed successfully for sample {sample}.")
    else:
        print(f"[WARNING] QC step executed but expected output files are missing for sample {sample}.")
        print(f"[DEBUG] Expected files:")
        for f in [qc_html_1, qc_zip_1, qc_html_2, qc_zip_2]:
            exists = "?" if os.path.exists(f) else "?"
            print(f"[DEBUG]   {exists} {f}")
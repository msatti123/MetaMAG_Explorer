import os
import subprocess
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

# Global sample mapping dictionary
SAMPLE_MAPPING = {}

def get_metawrap_perl_path():
    """
    Dynamically find the correct PERL5LIB path for metaWRAP environment
    """
    metawrap_env_path = get_tool_path('metawrap')
    if metawrap_env_path:
        # Extract the conda environment path
        env_path = metawrap_env_path.split('bin/activate')[0].strip()
        if env_path.startswith('source '):
            env_path = env_path[7:].strip()
        
        # Check multiple possible perl paths
        possible_paths = [
            f"{env_path}/lib/perl5/site_perl/5.22.0",
            f"{env_path}/lib/perl5",
            f"{env_path}/lib/perl5/site_perl",
            f"{env_path}/lib/perl5/5.22.0",
            f"{env_path}/lib/perl5/5.26.2",  # Sometimes different perl versions
            f"{env_path}/lib/perl5/5.30.0",
        ]
        
        existing_paths = []
        for path in possible_paths:
            if os.path.exists(path):
                existing_paths.append(path)
                print(f"[INFO] Found perl path: {path}")
        
        if existing_paths:
            return ":".join(existing_paths)
    
    return None

def check_maxbin2_installation(env_path):
    """
    Check if MaxBin2 is properly installed and find the correct path
    """
    # Common locations for run_MaxBin.pl
    possible_locations = [
        f"{env_path}/bin/run_MaxBin.pl",
        f"{env_path}/scripts/run_MaxBin.pl",
        "/usr/local/bin/run_MaxBin.pl",
        # Add the actual MaxBin2 installation path if known
    ]
    
    for loc in possible_locations:
        if os.path.exists(loc):
            print(f"[INFO] Found run_MaxBin.pl at: {loc}")
            return os.path.dirname(loc)
    
    # Try to find it using which in the environment
    try:
        result = subprocess.run(
            f"source {get_tool_path('metawrap')} && which run_MaxBin.pl",
            shell=True, capture_output=True, text=True, executable='/bin/bash'
        )
        if result.returncode == 0 and result.stdout.strip():
            print(f"[INFO] Found run_MaxBin.pl via which: {result.stdout.strip()}")
            return os.path.dirname(result.stdout.strip())
    except:
        pass
    
    return None

def assign_sample_ids(sample_list, project_output):
    """
    Assigns numeric IDs to samples and saves/updates a mapping file.
    Handles parallel processing by updating the mapping file with new samples.
    """
    reads_dir = f"{project_output}/reads"
    ensure_directory_exists(reads_dir)

    mapping_file = f"{reads_dir}/sample_mapping.txt"
    
    # Load existing mapping if it exists
    if os.path.exists(mapping_file):
        print(f"[INFO] Loading existing sample mapping from: {mapping_file}")
        with open(mapping_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    SAMPLE_MAPPING[parts[1]] = int(parts[0])
    
    # Find the highest existing ID
    max_id = max(SAMPLE_MAPPING.values()) if SAMPLE_MAPPING else 0
    
    # Check for new samples that need IDs
    new_samples = []
    for sample in sample_list:
        if sample not in SAMPLE_MAPPING:
            new_samples.append(sample)
    
    # If there are new samples, assign IDs and update the file
    if new_samples:
        print(f"[INFO] Found {len(new_samples)} new samples to add to mapping")
        
        # Use file locking to prevent race conditions in parallel execution
        import fcntl
        
        # Re-read the file with lock to ensure we have the latest data
        with open(mapping_file, "a+") as f:
            # Lock the file
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            
            # Move to beginning and re-read to get any updates from other processes
            f.seek(0)
            temp_mapping = {}
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    temp_mapping[parts[1]] = int(parts[0])
            
            # Update our global mapping
            SAMPLE_MAPPING.update(temp_mapping)
            
            # Recalculate max_id after re-reading
            max_id = max(SAMPLE_MAPPING.values()) if SAMPLE_MAPPING else 0
            
            # Move to end of file for appending
            f.seek(0, 2)
            
            # Assign IDs to truly new samples
            for sample in sample_list:
                if sample not in SAMPLE_MAPPING:
                    max_id += 1
                    SAMPLE_MAPPING[sample] = max_id
                    f.write(f"{max_id}\t{sample}\n")
                    print(f"[INFO] Assigned ID {max_id} to sample {sample}")
            
            # Release the lock (happens automatically when file closes)
            
    print(f"[INFO] Total samples in mapping: {len(SAMPLE_MAPPING)}")
    
    # Verify all requested samples have mappings
    missing_samples = [s for s in sample_list if s not in SAMPLE_MAPPING]
    if missing_samples:
        print(f"[ERROR] Failed to assign IDs to samples: {', '.join(missing_samples)}")
        raise ValueError(f"Sample ID assignment failed for: {', '.join(missing_samples)}")

def unzip_and_rename_reads(sample, project_output, co_binning=False):
    """
    Copies, unzips, and renames fastq.gz files before binning.
    Assigns numeric IDs dynamically and stores a mapping file.
    """
    host_removal_dir = f"{project_output}/Host_Removal"
    reads_dir = f"{project_output}/reads"
    ensure_directory_exists(reads_dir)

    if sample not in SAMPLE_MAPPING:
        print(f"[ERROR] Sample {sample} is missing from mapping. Cannot process.")
        return None, None

    sample_id = SAMPLE_MAPPING[sample]  # Get numeric ID from mapping
    gz_r1 = f"{host_removal_dir}/{sample}_unmapped_R1.fastq.gz"
    gz_r2 = f"{host_removal_dir}/{sample}_unmapped_R2.fastq.gz"

    # Correctly renamed output files
    renamed_r1 = f"{reads_dir}/reads{sample_id}_1.fastq"
    renamed_r2 = f"{reads_dir}/reads{sample_id}_2.fastq"

    print(f"[DEBUG] Co-binning: {co_binning}, Checking files: {gz_r1}, {gz_r2}")

    # Copy, Unzip & Rename if files do not exist
    if not os.path.exists(renamed_r1) or not os.path.exists(renamed_r2):
        if os.path.exists(gz_r1) and os.path.exists(gz_r2):
            print(f"[INFO] Copying and unzipping reads for sample {sample} (ID {sample_id})")
            run_command(f"gunzip -c {gz_r1} > {renamed_r1}")
            run_command(f"gunzip -c {gz_r2} > {renamed_r2}")
        else:
            print(f"[ERROR] Missing input files for sample {sample}: {gz_r1}, {gz_r2}")
            return None, None
    else:
        print(f"[INFO] Read files already exist for sample {sample} (ID {sample_id}), skipping unzip")

    return renamed_r1, renamed_r2

def check_existing_alignments(output_dir, samples, sample_mapping):
    """
    Check which alignment files already exist to avoid re-processing.
    Returns a list of samples that need alignment.
    """
    work_files_dir = f"{output_dir}/work_files"
    samples_needing_alignment = []
    
    for sample in samples:
        if sample not in sample_mapping:
            continue
            
        sample_id = sample_mapping[sample]
        # Check for sorted BAM file which indicates completed alignment
        sorted_bam = f"{work_files_dir}/reads{sample_id}.bam"
        
        if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
            print(f"[INFO] Alignment already exists for sample {sample} (ID {sample_id}), will skip")
        else:
            samples_needing_alignment.append(sample)
            
    return samples_needing_alignment

def setup_environment_vars():
    """
    Set up all necessary environment variables for metaWRAP
    """
    env_vars = []
    
    # Get dynamic PERL5LIB path
    perl_lib = get_metawrap_perl_path()
    if perl_lib:
        env_vars.append(f"export PERL5LIB={perl_lib}")
        print(f"[INFO] Setting PERL5LIB to: {perl_lib}")
    
    # Get metawrap environment path for PATH updates
    metawrap_env_path = get_tool_path('metawrap')
    if metawrap_env_path:
        env_path = metawrap_env_path.split('bin/activate')[0].strip()
        if env_path.startswith('source '):
            env_path = env_path[7:].strip()
        
        # Check for MaxBin2
        maxbin_path = check_maxbin2_installation(env_path)
        if maxbin_path:
            env_vars.append(f"export PATH={maxbin_path}:$PATH")
            print(f"[INFO] Adding MaxBin2 to PATH: {maxbin_path}")
    
    return " && ".join(env_vars)

def run_binning(sample_or_samples, cpus, memory, time, methods, project_input, project_output):
    """
    Runs binning methods using metaWRAP. Ensures all required read files exist before execution.
    """
    methods = methods or []
    methods_str = " ".join([f"--{method}" for method in methods])

    # Set up environment variables
    env_setup = setup_environment_vars()

    if isinstance(sample_or_samples, list):  # Co-assembly mode
        output_dir = f"{project_output}/Binning/Coassembly/Selected"
        ensure_directory_exists(output_dir)
        contigs = f"{project_output}/Assembly/MEGAHIT/coassembly.contigs.fa"

        # Assign numeric IDs to all samples before processing
        assign_sample_ids(sample_or_samples, project_output)

        # Check which samples need alignment
        samples_needing_alignment = check_existing_alignments(output_dir, sample_or_samples, SAMPLE_MAPPING)
        
        print(f"[INFO] Total samples: {len(sample_or_samples)}")
        print(f"[INFO] Samples needing alignment: {len(samples_needing_alignment)}")
        
        # Copy, Unzip & Rename Reads for Co-Binning
        read_files = []
        missing_samples = []
        
        for sample in sample_or_samples:
            r1, r2 = unzip_and_rename_reads(sample, project_output, co_binning=True)
            if r1 is None or r2 is None:
                print(f"[ERROR] Missing read files for sample {sample}.")
                missing_samples.append(sample)
            else:
                read_files.extend([r1, r2])

        # If ANY sample is missing, stop execution
        if missing_samples:
            print(f"[CRITICAL ERROR] The following samples have missing read files: {', '.join(missing_samples)}")
            print("[CRITICAL ERROR] Co-binning cannot proceed. Please check the input data.")
            exit(1)  # Force exit with error

        read_files_str = " ".join(read_files)
        
        # Create a file listing samples that need alignment
        if samples_needing_alignment:
            alignment_list_file = f"{output_dir}/samples_to_align.txt"
            with open(alignment_list_file, "w") as f:
                for sample in samples_needing_alignment:
                    f.write(f"{sample}\n")

    else:  # Single-sample mode
        output_dir = f"{project_output}/Binning/Single_Sample/Selected/{sample_or_samples}"
        ensure_directory_exists(output_dir)
        contigs = f"{project_output}/Assembly/IDBA/{sample_or_samples}_idba/scaffold.fa"

        assign_sample_ids([sample_or_samples], project_output)

        r1, r2 = unzip_and_rename_reads(sample_or_samples, project_output, co_binning=False)
        if r1 is None or r2 is None:
            print(f"[ERROR] Missing input files for sample {sample_or_samples}. Skipping binning.")
            return

        read_files_str = f"{r1} {r2}"

    # Build the binning command with proper environment setup
    if isinstance(sample_or_samples, list) and len(sample_or_samples) > 10:
        print(f"[INFO] Large dataset detected ({len(sample_or_samples)} samples). Using optimized binning approach.")
        
        # For large datasets, run all methods together with timeout
        binning_command = (
            f"bash -c '"
            f"source {get_tool_path('metawrap')} && "
            f"{env_setup} && "
            f"echo \"[INFO] Environment setup complete\" && "
            f"echo \"[INFO] PATH=$PATH\" && "
            f"echo \"[INFO] PERL5LIB=$PERL5LIB\" && "
            f"# Test MaxBin2 availability\n"
            f"if [[ \"{methods_str}\" == *\"--maxbin2\"* ]]; then "
            f"  which run_MaxBin.pl && echo \"[INFO] MaxBin2 found\" || echo \"[WARNING] MaxBin2 not found\"; "
            f"fi && "
            f"timeout 48h metaWRAP binning {methods_str} "
            f"-t {cpus} -m {memory} "
            f"-a {contigs} "
            f"-o {output_dir} "
            f"{read_files_str} "
            f"2>&1 | tee {output_dir}/binning_progress.log"
            f"'"
        )
    else:
        # For smaller datasets, run methods individually for better error handling
        binning_command = (
            f"bash -c '"
            f"source {get_tool_path('metawrap')} && "
            f"{env_setup} && "
            f"echo \"[INFO] Environment setup complete\" && "
            f"echo \"[INFO] PATH=$PATH\" && "
            f"echo \"[INFO] PERL5LIB=$PERL5LIB\" && "
        )
        
        # Add each method individually
        if "metabat2" in methods:
            binning_command += (
                f"echo \"[INFO] Running MetaBAT2...\"; "
                f"metaWRAP binning --metabat2 -t {cpus} -a {contigs} -o {output_dir} {read_files_str} || "
                f"echo \"[WARNING] MetaBAT2 failed, continuing...\"; "
            )
        
        if "maxbin2" in methods:
            binning_command += (
                f"echo \"[INFO] Running MaxBin2...\"; "
                f"# First check if MaxBin2 is available\n"
                f"if which run_MaxBin.pl > /dev/null 2>&1; then "
                f"  metaWRAP binning --maxbin2 -t {cpus} -a {contigs} -o {output_dir} {read_files_str} || "
                f"  echo \"[WARNING] MaxBin2 failed, continuing...\"; "
                f"else "
                f"  echo \"[ERROR] MaxBin2 not found. Skipping MaxBin2 binning.\"; "
                f"  echo \"[INFO] To fix this, ensure MaxBin2 is installed in the metawrap environment:\"; "
                f"  echo \"[INFO] conda activate metawrap_env && conda install -c bioconda maxbin2\"; "
                f"fi; "
            )
        
        if "concoct" in methods:
            binning_command += (
                f"echo \"[INFO] Running CONCOCT...\"; "
                f"metaWRAP binning --concoct -t {cpus} -a {contigs} -o {output_dir} {read_files_str} || "
                f"echo \"[WARNING] CONCOCT failed, continuing...\"; "
            )
        
        binning_command += (
            f"echo \"[INFO] Binning completed.\" "
            f"'"
        )

    print(f"[DEBUG] Running binning command: {binning_command}")
    
    # Create a status file to track progress
    status_file = f"{output_dir}/binning_status.txt"
    with open(status_file, "w") as f:
        f.write(f"Started: {os.popen('date').read().strip()}\n")
        f.write(f"Samples: {len(sample_or_samples) if isinstance(sample_or_samples, list) else 1}\n")
        f.write(f"Methods: {methods_str}\n")
    
    # Run the command
    run_command(binning_command)
    
    # Update status
    with open(status_file, "a") as f:
        f.write(f"Completed: {os.popen('date').read().strip()}\n")
        
    # Check results
    print(f"[INFO] Checking binning results in {output_dir}")
    for method in methods:
        method_dir = f"{output_dir}/{method.upper()}_bins"
        if os.path.exists(method_dir):
            bin_count = len([f for f in os.listdir(method_dir) if f.endswith('.fa')])
            print(f"[INFO] {method.upper()}: Found {bin_count} bins")
        else:
            print(f"[WARNING] {method.upper()}: No output directory found")
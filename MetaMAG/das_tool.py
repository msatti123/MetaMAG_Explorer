import os
import re
import sys
import shutil
from MetaMAG.utils import ensure_directory_exists, run_command
import subprocess
from MetaMAG.config import config, get_tool_path, get_tool_command

# Read paths from config.
fasta_to_contig2bin = get_tool_path("DAS_Tool_Fasta2Contig")
if not fasta_to_contig2bin:
    print(f"[ERROR] Das_Tool_Fasta2Contig tool not found in configuration")
    sys.exit(1)
das_activate = get_tool_path("das_tool")
if not das_activate:
    print(f"[ERROR] Das_Activate tool not found in configuration")
    sys.exit(1)

r"""
Copies files from source_dir to dest_dir that match the pattern:
  ^bin\.[0-9]+\.fa$
Only these FASTA files (e.g. bin.62.fa, bin.63.fa, etc.) will be copied.
"""
def filter_valid_bins(source_dir, dest_dir):
    ensure_directory_exists(dest_dir)
    pattern = re.compile(r'^bin\.[0-9]+\.fa$')
    valid_files = []
    for f in os.listdir(source_dir):
        full_path = os.path.join(source_dir, f)
        if os.path.isfile(full_path) and pattern.match(f):
            shutil.copy(full_path, dest_dir)
            valid_files.append(f)
    return valid_files

def run_das_tool_single(sample, cpus, memory, time, project_input, project_output, score_threshold=0.5):
    """
    Runs DAS Tool for single-sample binning.
    For each binning method, it first checks if bins exist, then filters the bin folder
    so that only FASTA files with numeric bin names (e.g. bin.62.fa) are retained.
    It then generates a TSV file using Fasta_to_Contig2Bin.sh on the filtered directory
    and runs DAS Tool with only the methods that have valid bins.
    All outputs are saved under:
      {project_output}/Bin_Refinement/dastool/Single_Sample/{sample}
    
    Args:
        sample: Sample ID
        cpus: Number of CPUs to use
        memory: Memory allocation
        time: Time allocation
        project_input: Path to project input directory
        project_output: Path to project output directory
        score_threshold: DAS Tool score threshold (default: 0.5)
    """
    print(f"Running DAS Tool (Single-Sample) for {sample}")

    base_output = f"{project_output}/Bin_Refinement/dastool/Single_Sample/{sample}"
    ensure_directory_exists(base_output)

    # Original bin directories for each method.
    metabat_bins = f"{project_output}/Binning/Single_Sample/Selected/{sample}/metabat2_bins"
    maxbin_bins  = f"{project_output}/Binning/Single_Sample/Selected/{sample}/maxbin2_bins"
    concoct_bins = f"{project_output}/Binning/Single_Sample/Selected/{sample}/concoct_bins"

    # Check which methods have valid bins
    valid_methods = []
    valid_bin_dirs = []
    method_labels = []
    
    # Define temporary filtered directories
    filtered_metabat = f"{base_output}/filtered_metabat"
    filtered_maxbin  = f"{base_output}/filtered_maxbin2"
    filtered_concoct = f"{base_output}/filtered_concoct"
    
    # Define TSV file names
    metabat_tsv = f"{base_output}/metabat_{sample}_contigs2bin.tsv"
    maxbin_tsv  = f"{base_output}/maxbin2_{sample}_contigs2bin.tsv"
    concoct_tsv = f"{base_output}/concoct_{sample}_contigs2bin.tsv"
    
    # Check metabat
    if os.path.exists(metabat_bins) and os.listdir(metabat_bins):
        valid_files = filter_valid_bins(metabat_bins, filtered_metabat)
        if valid_files:
            cmd_metabat = f"bash -c '{fasta_to_contig2bin} -i {filtered_metabat} -e fa > {metabat_tsv}'"
            run_command(cmd_metabat)
            valid_methods.append(metabat_tsv)
            method_labels.append("metabat")
            print(f"[INFO] Found {len(valid_files)} valid metabat bins for {sample}")
        else:
            print(f"[WARNING] No valid metabat bins found for {sample}")
    else:
        print(f"[WARNING] metabat bins directory is missing or empty for {sample}: {metabat_bins}")
    
    # Check maxbin
    if os.path.exists(maxbin_bins) and os.listdir(maxbin_bins):
        valid_files = filter_valid_bins(maxbin_bins, filtered_maxbin)
        if valid_files:
            cmd_maxbin = f"bash -c '{fasta_to_contig2bin} -i {filtered_maxbin} -e fa > {maxbin_tsv}'"
            run_command(cmd_maxbin)
            valid_methods.append(maxbin_tsv)
            method_labels.append("maxbin2")
            print(f"[INFO] Found {len(valid_files)} valid maxbin2 bins for {sample}")
        else:
            print(f"[WARNING] No valid maxbin2 bins found for {sample}")
    else:
        print(f"[WARNING] maxbin2 bins directory is missing or empty for {sample}: {maxbin_bins}")
    
    # Check concoct
    if os.path.exists(concoct_bins) and os.listdir(concoct_bins):
        valid_files = filter_valid_bins(concoct_bins, filtered_concoct)
        if valid_files:
            cmd_concoct = f"bash -c '{fasta_to_contig2bin} -i {filtered_concoct} -e fa > {concoct_tsv}'"
            run_command(cmd_concoct)
            valid_methods.append(concoct_tsv)
            method_labels.append("concoct")
            print(f"[INFO] Found {len(valid_files)} valid concoct bins for {sample}")
        else:
            print(f"[WARNING] No valid concoct bins found for {sample}")
    else:
        print(f"[WARNING] concoct bins directory is missing or empty for {sample}: {concoct_bins}")
    
    # Check if we have at least one valid method
    if not valid_methods:
        print(f"[ERROR] No valid bins found for any method for {sample}. Skipping DAS Tool.")
        return
    
    # Define contigs file and output for DAS Tool.
    contigs = f"{project_output}/Assembly/IDBA/{sample}_idba/scaffold.fa"
    das_tool_output = f"{base_output}/{sample}_DAS"
    
    # Run DAS Tool using only the valid methods
    das_tool_cmd = (
        f"bash -c 'source {das_activate} && "
        f"DAS_Tool -i {','.join(valid_methods)} -l {','.join(method_labels)} "
        f"-c {contigs} -o {das_tool_output} --write_bin_evals --write_bins --threads {cpus} "
        f"--score_threshold {score_threshold}'"
    )
    
    try:
        run_command(das_tool_cmd)
        print(f"[INFO] DAS Tool (Single-Sample) completed for {sample} using methods: {', '.join(method_labels)}")
        print(f"[INFO] Results saved in {das_tool_output}")
    except Exception as e:
        print(f"[WARNING] DAS Tool failed for {sample}: {str(e)}")
        print(f"[WARNING] This may be due to low-quality bins. Continuing pipeline...")

def run_das_tool_coassembly(cpus, memory, time, project_input, project_output, score_threshold=0.5):
    """
    Runs DAS Tool for co-assembly binning.
    For each binner, it checks if bins exist, then filters the bin folder so that only FASTA files 
    with numeric bin names are retained, generates a TSV file using Fasta_to_Contig2Bin.sh on the 
    filtered directory, then runs DAS Tool using only the methods that have valid bins.
    All outputs are saved under:
      {project_output}/Bin_Refinement/dastool/Coassembly
      
    Args:
        cpus: Number of CPUs to use
        memory: Memory allocation
        time: Time allocation
        project_input: Path to project input directory
        project_output: Path to project output directory
        score_threshold: DAS Tool score threshold (default: 0.5)
    """
    print("Running DAS Tool (Co-Assembly)")

    base_output = f"{project_output}/Bin_Refinement/dastool/Coassembly"
    ensure_directory_exists(base_output)

    metabat_bins = f"{project_output}/Binning/Coassembly/Selected/metabat2_bins"
    maxbin_bins  = f"{project_output}/Binning/Coassembly/Selected/maxbin2_bins"
    concoct_bins = f"{project_output}/Binning/Coassembly/Selected/concoct_bins"

    # Check which methods have valid bins
    valid_methods = []
    method_labels = []
    
    # Create temporary filtered directories.
    filtered_metabat = f"{base_output}/filtered_metabat"
    filtered_maxbin  = f"{base_output}/filtered_maxbin2"
    filtered_concoct = f"{base_output}/filtered_concoct"
    
    # Define TSV file names
    metabat_tsv = f"{base_output}/metabat_coassembly_contigs2bin.tsv"
    maxbin_tsv  = f"{base_output}/maxbin2_coassembly_contigs2bin.tsv"
    concoct_tsv = f"{base_output}/concoct_coassembly_contigs2bin.tsv"
    
    # Check metabat
    if os.path.exists(metabat_bins) and os.listdir(metabat_bins):
        valid_files = filter_valid_bins(metabat_bins, filtered_metabat)
        if valid_files:
            cmd_metabat = f"bash -c '{fasta_to_contig2bin} -i {filtered_metabat} -e fa > {metabat_tsv}'"
            run_command(cmd_metabat)
            valid_methods.append(metabat_tsv)
            method_labels.append("metabat")
            print(f"[INFO] Found {len(valid_files)} valid metabat bins for coassembly")
        else:
            print("[WARNING] No valid metabat bins found for coassembly")
    else:
        print(f"[WARNING] metabat bins directory is missing or empty: {metabat_bins}")
    
    # Check maxbin
    if os.path.exists(maxbin_bins) and os.listdir(maxbin_bins):
        valid_files = filter_valid_bins(maxbin_bins, filtered_maxbin)
        if valid_files:
            cmd_maxbin = f"bash -c '{fasta_to_contig2bin} -i {filtered_maxbin} -e fa > {maxbin_tsv}'"
            run_command(cmd_maxbin)
            valid_methods.append(maxbin_tsv)
            method_labels.append("maxbin2")
            print(f"[INFO] Found {len(valid_files)} valid maxbin2 bins for coassembly")
        else:
            print("[WARNING] No valid maxbin2 bins found for coassembly")
    else:
        print(f"[WARNING] maxbin2 bins directory is missing or empty: {maxbin_bins}")
    
    # Check concoct
    if os.path.exists(concoct_bins) and os.listdir(concoct_bins):
        valid_files = filter_valid_bins(concoct_bins, filtered_concoct)
        if valid_files:
            cmd_concoct = f"bash -c '{fasta_to_contig2bin} -i {filtered_concoct} -e fa > {concoct_tsv}'"
            run_command(cmd_concoct)
            valid_methods.append(concoct_tsv)
            method_labels.append("concoct")
            print(f"[INFO] Found {len(valid_files)} valid concoct bins for coassembly")
        else:
            print("[WARNING] No valid concoct bins found for coassembly")
    else:
        print(f"[WARNING] concoct bins directory is missing or empty: {concoct_bins}")
    
    # Check if we have at least one valid method
    if not valid_methods:
        print("[ERROR] No valid bins found for any method for coassembly. Skipping DAS Tool.")
        return
    
    contigs = f"{project_output}/Assembly/MEGAHIT/coassembly.contigs.fa"
    das_tool_output = f"{base_output}/Coassembly_DAS"
    
    # Run DAS Tool using only the valid methods
    das_tool_cmd = (
        f"bash -c 'source {das_activate} && "
        f"DAS_Tool -i {','.join(valid_methods)} -l {','.join(method_labels)} "
        f"-c {contigs} -o {das_tool_output} --write_bin_evals --write_bins --threads {cpus} "
        f"--score_threshold {score_threshold}'"
    )
    
    try:
        run_command(das_tool_cmd)
        print(f"[INFO] DAS Tool (Co-Assembly) completed using methods: {', '.join(method_labels)}")
        print(f"[INFO] Results saved in {das_tool_output}")
    except Exception as e:
        print(f"[WARNING] DAS Tool failed for co-assembly: {str(e)}")
        print(f"[WARNING] This may be due to low-quality bins. Continuing pipeline...")

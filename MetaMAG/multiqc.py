#!/usr/bin/env python3
"""
MultiQC module for MetaMAG pipeline.
Aggregates QC reports from various bioinformatics tools.
"""

import os
import glob
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run(sample, cpus, memory, time, project_input, project_output, **kwargs):
    """
    Run MultiQC to aggregate QC reports.
    
    Args:
        sample (str or list): Sample ID or list of sample IDs
        cpus (int): Number of CPU threads to use
        memory (str): Memory allocation
        time (str): Time limit for the step
        project_input (str): Input directory from project_config.yaml
        project_output (str): Output directory from project_config.yaml
        **kwargs: Additional keyword arguments
            - input_dir (str): Specific input directory to search for QC reports
    
    Returns:
        str: Path to the MultiQC report
    """
    # Validate and extract input directory
    input_dir = kwargs.get('input_dir')
    if not input_dir:
        print("[ERROR] No input directory provided for MultiQC")
        return None

    # Ensure sample is a string (extract if list)
    if isinstance(sample, list):
        sample = sample[0] if sample else None

    # Define output directory
    output_dir = os.path.join(project_output, "multiqc")
    ensure_directory_exists(output_dir)

    # MultiQC command
    multiqc = get_tool_path("multiqc")
    if not multiqc:
        print("[ERROR] MultiQC tool not found in configuration")
        return None

    # Construct command
    command = (
        f"{multiqc} {input_dir} "
        f"-o {output_dir} "
        f"-n multiqc_report "
        f"-f "
        f"--num-cpu-threads {cpus}"
    )

    print(f"[INFO] Running MultiQC command: {command}")
    
    # Run MultiQC
    run_command(command)

    # Find the generated report
    report_patterns = [
        os.path.join(output_dir, "multiqc_report.html"),
        os.path.join(output_dir, "multiqc_report_*.html")
    ]
    
    for pattern in report_patterns:
        matching_reports = glob.glob(pattern)
        if matching_reports:
            print(f"[INFO] MultiQC report generated: {matching_reports[0]}")
            return matching_reports[0]
    
    print("[WARNING] No MultiQC report was generated.")
    return None
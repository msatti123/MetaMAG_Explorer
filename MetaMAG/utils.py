import os
import subprocess
import logging
import sys
import shutil
from datetime import datetime
from MetaMAG.config import config

# Configure logging with real-time flush
log_file = "pipeline.log"
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_file, mode="a"),  # Log to file
        logging.StreamHandler(sys.stdout)         # Log to console in real-time
    ]
)

# Ensure Python prints logs immediately without buffering
sys.stdout.reconfigure(line_buffering=True)

def ensure_directory_exists(directory):
    """
    Ensures the specified directory exists. If it doesn't, it will be created.
    :param directory: Path to the directory to ensure.
    """
    os.makedirs(directory, exist_ok=True)
    logging.info(f"[DIRECTORY] Ensured directory exists: {directory}")
    sys.stdout.flush()

def get_tool_command(tool_name, default_path=None):
    """
    Get the command to run a tool, handling both direct paths and conda environments.
    
    Args:
        tool_name (str): Name of the tool
        default_path (str, optional): Default path to use if not found in config
    
    Returns:
        tuple: (command_prefix, tool_path)
            command_prefix: Any prefix needed (like conda activation)
            tool_path: Path to the tool executable
    """
    # First check if tool has a specific configuration in tool_configs
    if "tool_configs" in config and tool_name in config["tool_configs"] and config["tool_configs"][tool_name]:
        tool_config = config["tool_configs"][tool_name]
        
        if tool_config["requires_activation"]:
            # Tool requires conda environment activation
            logging.info(f"[TOOL] Using conda environment for {tool_name}: {tool_config['env_name']}")
            return (tool_config["activate_cmd"], tool_config["path"])
        else:
            # Tool can be run directly
            logging.info(f"[TOOL] Using direct path for {tool_name}: {tool_config['path']}")
            return ("", tool_config["path"])
    
    # Check if there's a direct path in tools section
    if tool_name in config["tools"]:
        tool_path = config["tools"][tool_name]
        
        # Check if this is an activation command (for backward compatibility)
        if tool_name.endswith("_activate"):
            # This is an environment activation command
            activate_cmd = tool_path
            # The actual tool name is without _activate
            actual_tool = tool_name.replace("_activate", "")
            tool_path = config["tools"].get(actual_tool, "")
            logging.info(f"[TOOL] Using activation command for {actual_tool}: {activate_cmd}")
            return (activate_cmd, tool_path)
        
        # Regular tool path
        logging.info(f"[TOOL] Using configured path for {tool_name}: {tool_path}")
        return ("", tool_path)
    
    # If tool not found in config, try to find it in PATH
    if default_path:
        logging.info(f"[TOOL] Using default path for {tool_name}: {default_path}")
        return ("", default_path)
    
    path_result = shutil.which(tool_name)
    if path_result:
        logging.info(f"[TOOL] Found {tool_name} in PATH: {path_result}")
        return ("", path_result)
    
    # Tool not found
    logging.error(f"[TOOL] Tool {tool_name} not found in config or PATH")
    return ("", tool_name)  # Return the tool name as a last resort

def run_command(command, log_prefix=""):
    """
    Executes a shell command and captures its output and errors in real-time.
    :param command: Command to execute.
    :param log_prefix: Prefix for log messages (e.g., step name).
    :return: The standard output of the command.
    :raises RuntimeError: If the command fails.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"[{log_prefix}] [COMMAND] Executing: {command}")
    sys.stdout.flush()

    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=1  # Line buffering
    )

    for line in process.stdout:
        logging.info(f"[{log_prefix}] {line.strip()}")
        sys.stdout.flush()  # Force immediate printing

    for line in process.stderr:
        logging.error(f"[{log_prefix}] [STDERR] {line.strip()}")
        sys.stderr.flush()  # Force immediate printing

    process.wait()

    if process.returncode != 0:
        raise RuntimeError(f"[{log_prefix}] Command failed: {command}")

def run_tool(tool_name, args, default_path=None, log_prefix=""):
    """
    Run a tool, handling both direct paths and conda environments.
    
    Args:
        tool_name (str): Name of the tool
        args (str): Arguments to pass to the tool
        default_path (str, optional): Default path to use if not found in config
        log_prefix (str): Prefix for log messages
    
    Returns:
        None - the command is executed directly
    """
    prefix, tool_path = get_tool_command(tool_name, default_path)
    
    if prefix:
        # Tool requires environment activation
        command = f"bash -c 'source {prefix} && {tool_path} {args}'"
    else:
        # Tool can be run directly
        command = f"{tool_path} {args}"
    
    # Use the existing run_command function to execute and log
    run_command(command, log_prefix=log_prefix or tool_name)

def create_slurm_script(step_name, command, output_dir, sample, cpus, memory, time, binning_methods=None):
    """
    Creates a SLURM job script for the pipeline.
    :param step_name: Name of the pipeline step.
    :param command: Command to run in the SLURM script.
    :param output_dir: Directory to store SLURM logs.
    :param sample: Sample name for the script.
    :param cpus: Number of CPUs to allocate.
    :param memory: Memory allocation.
    :param time: Time allocation.
    :param binning_methods: Optional binning methods (only used in binning steps).
    :return: Path to the created SLURM script.
    """
    script_path = os.path.join(output_dir, f"{step_name}_{sample}.sh")
    log_dir = os.path.join(output_dir, "logs")
    ensure_directory_exists(log_dir)

    # Append binning_methods only if they are provided
    binning_methods_str = f"--binning_methods {' '.join(binning_methods)}" if binning_methods else ""
    # Build the full command by appending binning_methods_str if needed
    full_command = f"{command} {binning_methods_str}".strip()

    script_content = f"""#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n {cpus}
#SBATCH --mem={memory}
#SBATCH -t {time}
#SBATCH --output={log_dir}/{step_name}_{sample}_%j.out
#SBATCH --error={log_dir}/{step_name}_{sample}_%j.err

echo "Starting step: {step_name} for sample: {sample} at $(date)"
{full_command}
echo "Completed step: {step_name} for sample: {sample} at $(date)"
"""

    with open(script_path, 'w') as script_file:
        script_file.write(script_content)

    logging.info(f"[SLURM] Created SLURM script: {script_path}")
    return script_path

def check_existing_output(output_files, log_prefix=""):
    """
    Checks if the specified output files already exist.
    :param output_files: List of output file paths to check.
    :param log_prefix: Prefix for log messages (e.g., step name).
    :return: True if all files exist, False otherwise.
    """
    existing_files = [file for file in output_files if os.path.exists(file)]
    missing_files = [file for file in output_files if not os.path.exists(file)]

    if existing_files:
        logging.info(f"[{log_prefix}] [CHECK] Existing files: {', '.join(existing_files)}")
    if missing_files:
        logging.warning(f"[{log_prefix}] [CHECK] Missing files: {', '.join(missing_files)}")

    sys.stdout.flush()
    return not missing_files

def create_slurm_tool_command(tool_name, args, default_path=None):
    """
    Creates a command string for a tool that can be used in a SLURM script,
    handling both direct paths and conda environments.
    
    Args:
        tool_name (str): Name of the tool
        args (str): Arguments to pass to the tool
        default_path (str, optional): Default path to use if not found in config
    
    Returns:
        str: Command string to be used in a SLURM script
    """
    prefix, tool_path = get_tool_command(tool_name, default_path)
    
    if prefix:
        # Tool requires environment activation
        command = f"source {prefix} && {tool_path} {args}"
    else:
        # Tool can be run directly
        command = f"{tool_path} {args}"
    
    return command
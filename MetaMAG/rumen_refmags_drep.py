import os
import sys
import argparse
import subprocess
import yaml

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    os.makedirs(directory, exist_ok=True)

def get_tool_path(tool_name):
    """Get tool path from config.py - handles multiple import scenarios"""
    try:
        # Method 1: Try to import config from MetaMAG package
        from MetaMAG.config import config
        return config["tools"].get(tool_name)
    except ImportError:
        try:
            # Method 2: Fallback to import from current directory
            from config import config
            return config["tools"].get(tool_name)
        except ImportError:
            try:
                # Method 3: Try importing from parent directory
                sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
                from config import config
                return config["tools"].get(tool_name)
            except ImportError:
                print(f"[WARNING] Could not import config module for tool: {tool_name}")
                return None
    except (KeyError, TypeError) as e:
        print(f"[WARNING] Error accessing config for {tool_name}: {e}")
        return None

def yaml_get_tool_path(yaml_config, tool_name):
    """Get tool path from YAML config - safe extraction"""
    try:
        if yaml_config and isinstance(yaml_config, dict):
            tools = yaml_config.get("tools", {})
            if isinstance(tools, dict):
                return tools.get(tool_name)
        return None
    except (KeyError, TypeError, AttributeError) as e:
        print(f"[WARNING] Error accessing YAML config for {tool_name}: {e}")
        return None

def validate_drep_availability(drep_path=None):
    """Validate that dRep is available either via direct path or system PATH"""
    if drep_path:
        # Check if dRep executable exists
        if os.path.exists(drep_path):
            print(f"[INFO] dRep executable found: {drep_path}")
            return drep_path, True
        else:
            print(f"[WARNING] dRep executable not found: {drep_path}")
    
    # Check if dRep is available in system PATH
    try:
        result = subprocess.run(
            ["dRep", "--version"], 
            check=True, 
            capture_output=True, 
            text=True,
            timeout=10
        )
        print(f"[INFO] System dRep found: {result.stdout.strip()}")
        return None, True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        print(f"[ERROR] dRep not found in system PATH")
        return None, False

def validate_input_directory(dataset_dir):
    """Validate input directory and count .fa files"""
    if not os.path.exists(dataset_dir):
        print(f"[ERROR] Dataset directory does not exist: {dataset_dir}")
        return False, []
    
    if not os.path.isdir(dataset_dir):
        print(f"[ERROR] Dataset path is not a directory: {dataset_dir}")
        return False, []
    
    # Find .fa files (case insensitive)
    fa_files = []
    for f in os.listdir(dataset_dir):
        if f.lower().endswith(('.fa', '.fasta', '.fna')):
            fa_files.append(f)
    
    if not fa_files:
        print(f"[ERROR] No FASTA files (.fa, .fasta, .fna) found in {dataset_dir}")
        return False, []
    
    print(f"[INFO] Found {len(fa_files)} FASTA files")
    return True, fa_files

def create_genome_list(dataset_dir, fa_files, output_dir):
    """Create genome list file for dRep"""
    genome_list_file = os.path.join(output_dir, "genome_list.txt")
    
    try:
        with open(genome_list_file, "w") as f:
            for file in fa_files:
                full_path = os.path.join(dataset_dir, file)
                # Verify each file exists
                if os.path.exists(full_path):
                    f.write(f"{full_path}\n")
                else:
                    print(f"[WARNING] File not found, skipping: {full_path}")
        
        print(f"[INFO] Created genome list file: {genome_list_file}")
        return genome_list_file
    
    except IOError as e:
        print(f"[ERROR] Failed to create genome list file: {e}")
        return None

def run_drep_command(drep_command, output_dir):
    """Execute dRep command with proper error handling"""
    print(f"[INFO] Executing dRep command: {drep_command}")
    
    try:
        result = subprocess.run(
            drep_command, 
            shell=True, 
            check=True, 
            stderr=subprocess.PIPE, 
            stdout=subprocess.PIPE,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        print("[INFO] dRep completed successfully")
        if result.stdout:
            print(f"[INFO] dRep stdout: {result.stdout}")
        
        # Check if expected output files were created
        expected_files = [
            os.path.join(output_dir, "data_tables", "genomeInfo.csv"),
            os.path.join(output_dir, "dereplicated_genomes")
        ]
        
        for expected_file in expected_files:
            if os.path.exists(expected_file):
                print(f"[INFO] Created expected output: {expected_file}")
            else:
                print(f"[WARNING] Expected output not found: {expected_file}")
        
        return 0
        
    except subprocess.TimeoutExpired:
        print("[ERROR] dRep execution timed out (2 hours)")
        return 1
    except subprocess.CalledProcessError as e:
        print("[ERROR] dRep execution failed")
        print(f"[ERROR] Return code: {e.returncode}")
        if e.stdout:
            print(f"[ERROR] Command output: {e.stdout}")
        if e.stderr:
            print(f"[ERROR] Command error: {e.stderr}")
        return 1
    except Exception as e:
        print(f"[ERROR] Unexpected error during dRep execution: {e}")
        return 1

def run(dataset_dir, output_dir, cpus, config_file=None):
    """
    Main function to run dRep on rumen reference MAGs.
    
    Args:
        dataset_dir (str): Directory containing input genome files
        output_dir (str): Directory to store dRep output
        cpus (int): Number of CPUs to use
        config_file (str, optional): Path to additional configuration file
    
    Returns:
        int: 0 for success, 1 for failure
    """
    print("=" * 60)
    print("Starting dRep Analysis on Rumen Reference MAGs")
    print("=" * 60)
    print(f"Dataset directory: {dataset_dir}")
    print(f"Output directory: {output_dir}")
    print(f"CPUs: {cpus}")
    print(f"Config file: {config_file if config_file else 'None'}")
    print("=" * 60)
    
    # Ensure output directory exists
    ensure_directory_exists(output_dir)
    
    # Validate input directory and find FASTA files
    valid_input, fa_files = validate_input_directory(dataset_dir)
    if not valid_input:
        return 1
    
    # Get dRep activation path from configuration
    drep_path = None
    
    # Check YAML config file if provided
    if config_file and os.path.exists(config_file):
        print(f"[INFO] Loading configuration from: {config_file}")
        try:
            with open(config_file, 'r') as f:
                yaml_config = yaml.safe_load(f)
                if yaml_config and 'tools' in yaml_config:
                    drep_path = yaml_get_tool_path(yaml_config, "drep")
                    if drep_path:
                        print(f"[INFO] Found drep_activate in YAML config: {drep_path}")
        except Exception as e:
            print(f"[WARNING] Error loading config file {config_file}: {e}")
    
    # Fallback to config.py if no activation path found
    if not drep_path:
        print("[INFO] Trying to get drep_path from config.py")
        drep_path = get_tool_path("drep")
        if drep_path:
            print(f"[INFO] Found drep_path in config.py: {drep_path}")
    
    # Validate dRep availability
    validated_activate, drep_available = validate_drep_availability(drep_path)
    if not drep_available:
        print("[ERROR] dRep is not available via activation script or system PATH")
        print("[ERROR] Please install dRep or configure the activation script path")
        return 1
    
    # Create genome list file
    genome_list_file = create_genome_list(dataset_dir, fa_files, output_dir)
    if not genome_list_file:
        return 1
    
    # Prepare dRep command
    if validated_activate:
        print(f"[INFO] Using dRep activation script: {validated_activate}")
        drep_command = f"{validated_activate} dereplicate {output_dir} -g {genome_list_file} -p {cpus}"        
    else:
        print("[INFO] Using system dRep installation")
        drep_command = f"dRep dereplicate {output_dir} -g {genome_list_file} -p {cpus}"
    
    # Execute dRep
    exit_code = run_drep_command(drep_command, output_dir)
    
    print("=" * 60)
    if exit_code == 0:
        print("? dRep analysis completed successfully")
    else:
        print("? dRep analysis failed")
    print("=" * 60)
    
    return exit_code

def main():
    """Command line interface for dRep analysis"""
    parser = argparse.ArgumentParser(
        description="Run dRep on rumen reference MAGs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rumen_refmags_drep.py --dataset_dir /path/to/mags --output_dir /path/to/output --cpus 8
  python rumen_refmags_drep.py --dataset_dir /path/to/mags --output_dir /path/to/output --cpus 8 --config config.yaml
        """
    )
    
    parser.add_argument(
        "--dataset_dir", 
        type=str, 
        required=True, 
        help="Directory containing input FASTA files (.fa, .fasta, .fna)"
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        required=True, 
        help="Directory to store dRep output"
    )
    parser.add_argument(
        "--cpus", 
        type=int, 
        default=1, 
        help="Number of CPUs to use (default: 1)"
    )
    parser.add_argument(
        "--config", 
        type=str, 
        help="Path to YAML configuration file containing tool paths"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.cpus < 1:
        print("[ERROR] Number of CPUs must be at least 1")
        sys.exit(1)
    
    # Run dRep analysis
    exit_code = run(
        dataset_dir=args.dataset_dir, 
        output_dir=args.output_dir, 
        cpus=args.cpus,
        config_file=args.config
    )
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
import os
import sys
import shutil
import subprocess
import pandas as pd
import argparse

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    os.makedirs(directory, exist_ok=True)

def get_tool_path(tool_name):
    """Get tool path from config - handles multiple import scenarios"""
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

def validate_drep_availability(drep_activate=None):
    """Validate that dRep is available either via activation script or system PATH"""
    if drep_activate:
        # Check if activation script exists
        if os.path.exists(drep_activate):
            print(f"[INFO] dRep activation script found: {drep_activate}")
            return drep_activate, True
        else:
            print(f"[WARNING] dRep activation script not found: {drep_activate}")
    
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

def run_drep_command(drep_command, drep_output_dir):
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
            timeout=14400  # 4 hour timeout
        )
        
        print(f"[INFO] dRep completed successfully. Results saved in {drep_output_dir}")
        if result.stdout:
            print(f"[INFO] dRep stdout: {result.stdout}")
        return True
        
    except subprocess.TimeoutExpired:
        print("[ERROR] dRep execution timed out (4 hours)")
        return False
    except subprocess.CalledProcessError as e:
        print("[ERROR] dRep execution failed")
        print(f"[ERROR] Return code: {e.returncode}")
        if e.stdout:
            print(f"[ERROR] Command output: {e.stdout}")
        if e.stderr:
            print(f"[ERROR] Command error: {e.stderr}")
        return False
    except Exception as e:
        print(f"[ERROR] Unexpected error during dRep execution: {e}")
        return False

def run(ref_mags_dir, cpus, memory, time, project_input, project_output):
    """
    Combines project MAGs with reference MAGs, runs dRep, and separates unique MAGs.
    
    Args:
        ref_mags_dir (str): Directory containing reference MAGs
        cpus (int): Number of CPUs to use
        memory (str): Memory allocation
        time (str): Time allocation
        project_input (str): Project input directory
        project_output (str): Project output directory
    
    Returns:
        bool: True for success, False for failure
    """
    print("=" * 70)
    print("Starting Combine and Dereplicate MAGs Process")
    print("=" * 70)
    print(f"Reference MAGs directory: {ref_mags_dir}")
    print(f"CPUs: {cpus}")
    print(f"Memory: {memory}")
    print(f"Time: {time}")
    print(f"Project input: {project_input}")
    print(f"Project output: {project_output}")
    print("=" * 70)
    
    # Base directory from config
    base_dir = project_output
    
    # Input directory with project dereplicated genomes
    proj_genomes_dir = os.path.join(base_dir, "Bin_Refinement", "drep", "dRep_output", "dereplicated_genomes")
    
    # Reference MAGs directory from argument
    ref_genomes_dir = ref_mags_dir
    
    # Create directory structure
    temp_dir = os.path.join(base_dir, "Novel_Mags", "temp_combined")
    drep_output_dir = os.path.join(base_dir, "Novel_Mags", "rumen_drep_work")
    unique_proj_dir = os.path.join(base_dir, "Novel_Mags", "UniqueProjMags")
    unique_ref_dir = os.path.join(base_dir, "Novel_Mags", "UniqueRefMags")
    
    for directory in [temp_dir, drep_output_dir, unique_proj_dir, unique_ref_dir]:
        ensure_directory_exists(directory)
    
    # Step 1: Verify directories exist
    if not os.path.exists(proj_genomes_dir):
        print(f"[ERROR] Project genomes directory not found: {proj_genomes_dir}")
        print("[INFO] Checking alternative locations...")
        
        # Try alternative paths
        alt_paths = [
            os.path.join(base_dir, "Bin_Refinement", "dereplicate", "dRep_output", "dereplicated_genomes"),
            os.path.join(base_dir, "Bin_Refinement", "dRep", "dereplicated_genomes"),
            os.path.join(base_dir, "dereplicate", "dereplicated_genomes"),
            os.path.join(base_dir, "dRep_output", "dereplicated_genomes")
        ]
        
        found_alternative = False
        for alt_path in alt_paths:
            if os.path.exists(alt_path):
                print(f"[INFO] Found alternative project genomes directory: {alt_path}")
                proj_genomes_dir = alt_path
                found_alternative = True
                break
        
        if not found_alternative:
            print(f"[ERROR] Could not find project genomes directory in any expected location")
            return False
    
    if not os.path.exists(ref_genomes_dir):
        print(f"[ERROR] Reference MAGs directory not found: {ref_genomes_dir}")
        return False
    
    # Validate input directories have files
    proj_files = [f for f in os.listdir(proj_genomes_dir) if f.endswith(('.fa', '.fasta'))]
    ref_files = [f for f in os.listdir(ref_genomes_dir) if f.endswith(('.fa', '.fasta'))]
    
    if not proj_files:
        print(f"[ERROR] No FASTA files found in project genomes directory: {proj_genomes_dir}")
        return False
    
    if not ref_files:
        print(f"[ERROR] No FASTA files found in reference MAGs directory: {ref_genomes_dir}")
        return False
    
    print(f"[INFO] Found {len(proj_files)} project genomes and {len(ref_files)} reference genomes")
    
    # Check if MAGs have already been copied to the unique directories
    proj_files_in_unique = len([f for f in os.listdir(unique_proj_dir) if f.endswith(('.fa', '.fasta'))])
    ref_files_in_unique = len([f for f in os.listdir(unique_ref_dir) if f.endswith(('.fa', '.fasta'))])
    
    # Define the possible locations for dereplicated genomes
    derep_genomes_dir = os.path.join(drep_output_dir, "dereplicated_genomes")
    alt_derep_genomes_dir = os.path.join(drep_output_dir, "dRep_output", "dereplicated_genomes")
    
    # Check if the work is already done
    if proj_files_in_unique > 0 and ref_files_in_unique > 0:
        print(f"[INFO] Found {proj_files_in_unique} project MAGs and {ref_files_in_unique} reference MAGs already in output directories.")
        print("[INFO] Skipping processing as files already exist in output directories.")
        return True
    elif os.path.exists(derep_genomes_dir) and len([f for f in os.listdir(derep_genomes_dir) if f.endswith('.fa')]) > 0:
        print(f"[INFO] dRep has already been run. Found dereplicated genomes in {derep_genomes_dir}")
        print("[INFO] Will skip to separating MAGs step.")
        # Continue to the separating step
    elif os.path.exists(alt_derep_genomes_dir) and len([f for f in os.listdir(alt_derep_genomes_dir) if f.endswith('.fa')]) > 0:
        print(f"[INFO] dRep has already been run. Found dereplicated genomes in {alt_derep_genomes_dir}")
        derep_genomes_dir = alt_derep_genomes_dir
        print("[INFO] Will skip to separating MAGs step.")
        # Continue to the separating step
    else:
        # Need to run the full process
        
        # Step 2: Copy all genomes to temporary directory (if not already done)
        temp_files = [f for f in os.listdir(temp_dir) if f.endswith(('.fa', '.fasta'))]
        if len(temp_files) > 0:
            print(f"[INFO] Found {len(temp_files)} files already in temporary directory. Skipping copy step.")
            # Count the existing files
            proj_files_copied = len([f for f in temp_files if f.startswith('proj_')])
            ref_files_copied = len([f for f in temp_files if f.startswith('ref_')])
        else:
            print(f"[INFO] Copying project genomes to temporary directory: {temp_dir}")
            proj_files_copied = 0
            ref_files_copied = 0
            
            # Copy project genomes
            for filename in proj_files:
                src_path = os.path.join(proj_genomes_dir, filename)
                dest_path = os.path.join(temp_dir, f"proj_{filename}")
                try:
                    shutil.copy2(src_path, dest_path)
                    proj_files_copied += 1
                except Exception as e:
                    print(f"[WARNING] Failed to copy {filename}: {e}")
            
            # Copy reference genomes
            for filename in ref_files:
                src_path = os.path.join(ref_genomes_dir, filename)
                # Use a prefix to distinguish reference genomes
                dest_path = os.path.join(temp_dir, f"ref_{filename}")
                try:
                    shutil.copy2(src_path, dest_path)
                    ref_files_copied += 1
                except Exception as e:
                    print(f"[WARNING] Failed to copy {filename}: {e}")
            
            print(f"[INFO] Copied {proj_files_copied} project genomes and {ref_files_copied} reference genomes")
        
        # Step 3: Run dRep on the combined set
        print(f"[INFO] Running dRep on combined genomes")
        
        # Create genome list file
        genome_list_file = os.path.join(base_dir, "Novel_Mags", "combined_genome_list.txt")
        with open(genome_list_file, "w") as f:
            for filename in os.listdir(temp_dir):
                if filename.endswith(('.fa', '.fasta')):
                    file_path = os.path.join(temp_dir, filename)
                    f.write(f"{file_path}\n")
        
        print(f"[INFO] Created genome list file with {sum(1 for line in open(genome_list_file))} genomes")
        
        # Get dRep activation path from config
        drep_activate = get_tool_path('drep_activate')
        
        # Validate dRep availability
        validated_activate, drep_available = validate_drep_availability(drep_activate)
        if not drep_available:
            print("[ERROR] dRep is not available via activation script or system PATH")
            print("[ERROR] Please install dRep or configure the activation script path")
            return False
        
        # Prepare dRep command following the same pattern as the example
        if validated_activate:
            print(f"[INFO] Using dRep activation script: {validated_activate}")
            drep_command = (
                f"bash -c 'source {validated_activate} && "
                f"dRep dereplicate {drep_output_dir} -g {genome_list_file} -p {cpus} "
                f"--S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95'"
            )
        else:
            print("[INFO] Using system dRep installation")
            drep_command = (
                f"dRep dereplicate {drep_output_dir} -g {genome_list_file} -p {cpus} "
                f"--S_algorithm fastANI -comp 80 -con 10 -str 100 -strW 0 -pa 0.95"
            )
        
        # Run dRep command
        if not run_drep_command(drep_command, drep_output_dir):
            return False
    
    # Step 4: Find the dereplicated_genomes directory
    if not os.path.exists(derep_genomes_dir):
        if os.path.exists(alt_derep_genomes_dir):
            derep_genomes_dir = alt_derep_genomes_dir
        else:
            print(f"[ERROR] Could not find dereplicated genomes directory at:")
            print(f"  - {derep_genomes_dir}")
            print(f"  - {alt_derep_genomes_dir}")
            
            # List contents of drep_output_dir to help diagnose
            print(f"[DEBUG] Contents of {drep_output_dir}:")
            try:
                for item in os.listdir(drep_output_dir):
                    item_path = os.path.join(drep_output_dir, item)
                    if os.path.isdir(item_path):
                        print(f"  - {item}/ (directory)")
                        try:
                            subcontents = os.listdir(item_path)
                            if len(subcontents) > 5:  # If there are many files, just show a sample
                                print(f"    - Contains {len(subcontents)} items, including: {', '.join(subcontents[:5])}...")
                            else:
                                print(f"    - Contains: {', '.join(subcontents)}")
                        except Exception as e:
                            print(f"    - Error listing contents: {e}")
                    else:
                        print(f"  - {item} (file)")
            except Exception as e:
                print(f"[ERROR] Could not list contents of {drep_output_dir}: {e}")
            return False
    
    # Step 5: Separate unique MAGs by source if not already done
    if proj_files_in_unique == 0 and ref_files_in_unique == 0:
        print(f"[INFO] Separating MAGs from {derep_genomes_dir}")
        proj_winners = 0
        ref_winners = 0
        
        derep_files = [f for f in os.listdir(derep_genomes_dir) if f.endswith(('.fa', '.fasta'))]
        print(f"[INFO] Found {len(derep_files)} dereplicated genomes to separate")
        
        for filename in derep_files:
            src_path = os.path.join(derep_genomes_dir, filename)
            
            if filename.startswith('proj_'):
                # This is a project genome
                dest_filename = filename.replace('proj_', '', 1)  # Remove only the first occurrence
                dest_path = os.path.join(unique_proj_dir, dest_filename)
                try:
                    shutil.copy2(src_path, dest_path)
                    proj_winners += 1
                except Exception as e:
                    print(f"[WARNING] Failed to copy project MAG {filename}: {e}")
            elif filename.startswith('ref_'):
                # This is a reference genome
                dest_filename = filename.replace('ref_', '', 1)  # Remove only the first occurrence
                dest_path = os.path.join(unique_ref_dir, dest_filename)
                try:
                    shutil.copy2(src_path, dest_path)
                    ref_winners += 1
                except Exception as e:
                    print(f"[WARNING] Failed to copy reference MAG {filename}: {e}")
            else:
                print(f"[WARNING] Unknown MAG prefix for file: {filename}")
        
        print(f"[INFO] Separated unique MAGs: {proj_winners} project-specific and {ref_winners} reference-specific")
    else:
        print(f"[INFO] MAGs already separated: {proj_files_in_unique} project-specific and {ref_files_in_unique} reference-specific")
        proj_winners = proj_files_in_unique
        ref_winners = ref_files_in_unique
    
    # Get the counts if we skipped the copy step but need to include in summary
    if 'proj_files_copied' not in locals():
        if os.path.exists(temp_dir):
            temp_files = [f for f in os.listdir(temp_dir) if f.endswith(('.fa', '.fasta'))]
            proj_files_copied = len([f for f in temp_files if f.startswith('proj_')])
            ref_files_copied = len([f for f in temp_files if f.startswith('ref_')])
        else:
            # If we don't have exact counts, estimate based on the source directories
            proj_files_copied = len(proj_files)
            ref_files_copied = len(ref_files)
    
    # Create a simple summary report if it doesn't exist
    summary_report_path = os.path.join(base_dir, "Novel_Mags", "uniqueness_summary.txt")
    if not os.path.exists(summary_report_path):
        try:
            with open(summary_report_path, "w") as f:
                f.write("=== MAG Uniqueness Analysis Summary ===\n\n")
                f.write(f"Total input project MAGs: {proj_files_copied}\n")
                f.write(f"Total input reference MAGs: {ref_files_copied}\n")
                f.write(f"Unique project-specific MAGs: {proj_winners}\n")
                f.write(f"Unique reference-specific MAGs: {ref_winners}\n")
                f.write(f"\nProject uniqueness rate: {(proj_winners/proj_files_copied)*100:.1f}%\n")
                f.write(f"Reference uniqueness rate: {(ref_winners/ref_files_copied)*100:.1f}%\n")
                f.write(f"\nTotal dereplicated MAGs: {proj_winners + ref_winners}\n")
                f.write(f"Dereplication rate: {((proj_winners + ref_winners)/(proj_files_copied + ref_files_copied))*100:.1f}%\n")
            
            print(f"[INFO] Summary report created: {summary_report_path}")
        except Exception as e:
            print(f"[WARNING] Could not create summary report: {e}")
    else:
        print(f"[INFO] Summary report already exists: {summary_report_path}")
    
    print("=" * 70)
    print("? Combine and dereplicate MAGs process completed successfully")
    print(f"? Project-specific MAGs: {proj_winners}")
    print(f"? Reference-specific MAGs: {ref_winners}")
    print(f"? Results saved in: {base_dir}/Novel_Mags/")
    print("=" * 70)
    
    return True

def main():
    """Command line interface for MAG combination and dereplication"""
    parser = argparse.ArgumentParser(
        description="Combine project MAGs with reference MAGs, run dRep, and separate unique MAGs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rumen_drep.py --ref_mags_dir /path/to/ref_mags --cpus 16 --memory 100G --time 12:00:00 \\
                       --project_input /path/to/input --project_output /path/to/output
        """
    )
    
    parser.add_argument("--ref_mags_dir", type=str, required=True, 
                        help="Directory containing reference MAGs")
    parser.add_argument("--cpus", type=int, required=True, 
                        help="Number of CPUs to use")
    parser.add_argument("--memory", required=True, 
                        help="Memory allocation (e.g., 100G)")
    parser.add_argument("--time", required=True, 
                        help="Time allocation (e.g., 12:00:00)")
    parser.add_argument("--project_input", type=str, required=True, 
                        help="Project input directory")
    parser.add_argument("--project_output", type=str, required=True, 
                        help="Project output directory")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.cpus < 1:
        print("[ERROR] Number of CPUs must be at least 1")
        sys.exit(1)
    
    if not os.path.exists(args.project_output):
        print(f"[ERROR] Project output directory does not exist: {args.project_output}")
        sys.exit(1)
    
    # Run the main function
    success = run(
        args.ref_mags_dir,
        args.cpus,
        args.memory,
        args.time,
        args.project_input,
        args.project_output
    )
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
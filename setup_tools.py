#!/usr/bin/env python3
import os
import subprocess
import sys
import argparse
import shutil
import json

# Define tool groups and their components with tools from the config
TOOL_GROUPS = {
    "qc": ["fastqc", "multiqc"],
    "preprocessing": ["fastp", "bwa", "samtools"],
    "assembly": ["idba", "megahit", "metaquast"],
    "binning": ["jgi_summarize", "das_tool", "metawrap", "metabat2"],
    "evaluation": ["checkm2", "drep"],
    "taxonomy": ["gtdbtk", "kraken2", "bracken"],
    "annotation": ["eggnog-mapper","eggnog_db_dir", "dbcan", "dbcan_db_dir", "prodigal"],
    }

# Dictionary mapping tools to their executable names
TOOL_EXECUTABLE_MAP = {
    "eggnog-mapper": "emapper.py",
    "eggnog_mapper": "emapper.py",
    "dbcan": "run_dbcan.py",
    "das_tool": "DAS_Tool",
    "maxbin2": "run_MaxBin.pl",
    "drep": "dRep",
    "metaquast": "metaquast.py",
    "jgi_summarize": "jgi_summarize_bam_contig_depths",
    "metawrap": "python -c 'import metawrap'"
}
# Dictionary mapping tool names to their conda package names
TOOL_CONDA_PACKAGE_MAP = {
    "metaquast": "quast",  # metaquast.py comes with quast package
    "idba_fq2fa": None,    # comes with idba, skip during installation
    "DAS_Tool_Fasta2Contig": None,  # comes with das_tool, skip
    "eggnog_db_dir": None,  # database directory, not a package
    "dbcan_db_dir": None,   # database directory, not a package
    "jgi_summarize": "metabat2"  # comes with metabat2 package
}


# Define dependencies between steps
STEP_DEPENDENCIES = {
    "qc": [],
    "trimming": ["qc"],
    "host_removal": ["trimming"],
    "single_assembly": ["host_removal"],
    "co_assembly": ["host_removal"],
    "metaquast": ["single_assembly"],
    "metaquast_coassembly": ["co_assembly"],
    "single_binning": ["single_assembly"],
    "co_binning": ["co_assembly"],
    "single_bin_refinement": ["single_binning"],
    "co_bin_refinement": ["co_binning"],
    "dRep": ["single_bin_refinement", "co_bin_refinement"],
    "evaluation": ["dRep"],
    "gtdbtk": ["dRep"],
    "identify_novel_mags": ["gtdbtk"],
    "process_novel_mags": ["identify_novel_mags"],
    "add_mags_to_repo": ["process_novel_mags"],
    "build_kraken_db": ["add_mags_to_repo"],
    "eggnog_annotation": ["dRep"],
    "dbcan_annotation": ["dRep"],
    "functional_analysis": ["eggnog_annotation", "dbcan_annotation"],
    "mags_tree": ["gtdbtk"],
    "tree_visualization": ["mags_tree"],
    "kraken_abundance": ["build_kraken_db"],
    "abundance_estimation": ["kraken_abundance"]
}

# Map pipeline steps to tool groups
STEP_TO_TOOLS = {
    "qc": ["qc"],
    "trimming": ["preprocessing"],
    "host_removal": ["preprocessing"],
    "single_assembly": ["assembly"],
    "co_assembly": ["assembly"],
    "metaquast": ["assembly"],
    "metaquast_coassembly": ["assembly"],
    "single_binning": ["binning"],
    "co_binning": ["binning"],
    "single_bin_refinement": ["binning"],
    "co_bin_refinement": ["binning"],
    "dRep": ["evaluation"],
    "evaluation": ["evaluation"],
    "gtdbtk": ["taxonomy"],
    "identify_novel_mags": ["taxonomy"],
    "process_novel_mags": ["taxonomy"],
    "add_mags_to_repo": ["taxonomy"],
    "build_kraken_db": ["taxonomy"],
    "eggnog_annotation": ["annotation"],
    "dbcan_annotation": ["annotation"],
    "functional_analysis": ["annotation", "visualization"],
    "mags_tree": ["visualization"],
    "tree_visualization": ["visualization"],
    "kraken_abundance": ["taxonomy"],
    "abundance_estimation": ["taxonomy", "visualization"]
}

# Tool version check commands
TOOL_VERSION_COMMANDS = {
    "fastqc": "fastqc --version",
    "multiqc": "multiqc --version",
    "fastp": "fastp --version",
    "bwa": "bwa",
    "samtools": "samtools --version",
    "idba": "idba --version",
    "idba_fq2fa": "idba_fq2fa --help",
    "megahit": "megahit --version",
    "metaquast": "metaquast.py --version",
    "jgi_summarize": "jgi_summarize_bam_contig_depths --help",
    "metabat2": "metabat2 --help",
    "maxbin2": "run_MaxBin.pl -v",
    "concoct": "concoct --version",
    "das_tool": "DAS_Tool -v",
    "metawrap": "metawrap -v",
    "checkm2": "checkm2 -h",
    "drep": "dRep -h",
    "gtdbtk": "gtdbtk --version",
    "kraken2": "kraken2 --version",
    "bracken": "bracken -v",
    "eggnog-mapper": "emapper.py -v",
    "prodigal": "prodigal -h",
    "dbcan": "run_dbcan -v"
}

# Dictionary mapping tools to their executable names

# Conda environments for tools
TOOL_CONDA_ENVS = {
    "fastqc": "metamag",
    "multiqc": "metamag",
    "fastp": "metamag",
    "bwa": "metamag",
    "samtools": "metamag",
    "megahit": "metamag",
    "metaquast": "quast",
    "jgi_summarize": "metabat_env",
    "metabat2": "metabat_env",
    "maxbin2": "metamag",
    "concoct": "metamag",
    "das_tool": "das_env",
    "metawrap": "metawrap_env",
    "checkm2": "checkm_env",
    "drep": "metamag",
    "gtdbtk": "gtdbtk",
    "kraken2": "metamag",
    "bracken": "metamag",
    "eggnog-mapper": "metamag",
    "prodigal": "metamag",
    "dbcan": "dbcan"
}

def is_conda_available():
    """Check if conda is available in PATH."""
    return shutil.which("conda") is not None

def detect_perl_lib_path():
    """Detect the appropriate PERL5LIB path for the current system."""
    # Try to find the current Perl version
    try:
        result = subprocess.run(['perl', '-e', 'print $^V'], 
                               capture_output=True, text=True, check=True)
        perl_version = result.stdout.strip()
        # Convert from format like "v5.36.0" to "5.36"
        if perl_version.startswith('v'):
            perl_version = perl_version[1:]  # Remove the 'v' prefix
        # Extract major and minor version
        version_parts = perl_version.split('.')
        if len(version_parts) >= 2:
            perl_version = f"{version_parts[0]}.{version_parts[1]}"
    except Exception:
        # Default to a common version if detection fails
        perl_version = "5.36"
    
    # Check common Perl library paths
    possible_paths = [
        f"/opt/ghpc/lib64/perl-{perl_version}",
        f"/usr/lib/perl5/{perl_version}",
        f"/usr/local/lib/perl5/{perl_version}",
        f"/usr/local/lib/perl5/site_perl/{perl_version}",
        f"/usr/share/perl5",
        os.path.expanduser(f"~/perl5/lib/perl5")
    ]
    
    # Also check the current PERL5LIB environment variable
    current_perl5lib = os.environ.get("PERL5LIB", "")
    if current_perl5lib:
        possible_paths.insert(0, current_perl5lib)
    
    # Return the first path that exists
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # If no path is found, return a default or empty string
    return "/opt/ghpc/lib64/perl-5.36"  # Default fallback to the original value

def find_config_file():
    """
    Search for config.py in multiple locations
    """
    potential_paths = [
        # Current directory
        "config.py",
        "MetaMAG/config.py",
        
        # Script's directory
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.py"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "MetaMAG", "config.py"),
        
        # Absolute path in the current working directory
        os.path.join(os.getcwd(), "config.py"),
        os.path.join(os.getcwd(), "MetaMAG", "config.py"),
        
        # Home directory
        os.path.expanduser("~/config.py"),
        os.path.expanduser("~/MetaMAG/config.py")
    ]
    
    # Print search paths for debugging
    print("Searching for config.py in:")
    for path in potential_paths:
        print(f"  - {path}")
        if os.path.exists(path):
            print(f"    ? Found config file at {path}")
            return path
    
    return None

def check_env_exists(env_name):
    """Check if a conda environment exists."""
    try:
        # Run conda env list and get JSON output
        result = subprocess.run(
            ['conda', 'env', 'list', '--json'],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse JSON output
        env_list = json.loads(result.stdout)
        
        # Check if environment exists in the list
        env_names = [os.path.basename(env) for env in env_list["envs"]]
        return env_name in env_names
    except Exception as e:
        print("Error checking environments: {}".format(e))
        return False

def update_existing_environment(env_name, tools):
    """Update an existing conda environment with new tools.
    
    Supports tool version specifications like 'gtdbtk=2.3.2'
    """
    print("Updating existing environment '{}' with new tools...".format(env_name))
    
    # Map tool names to their conda package names
    conda_packages = []
    for tool in tools:
        # Check if tool has a special conda package name
        if tool in TOOL_CONDA_PACKAGE_MAP:
            conda_pkg = TOOL_CONDA_PACKAGE_MAP[tool]
            if conda_pkg:  # If not None, add it
                conda_packages.append(conda_pkg)
            # If None, skip this tool (it comes bundled with another)
        else:
            # Use the tool name as-is
            conda_packages.append(tool)
    
    # Remove duplicates
    conda_packages = list(set(conda_packages))
    
    if not conda_packages:
        print("No tools to install after filtering.")
        return True
    
    tools = conda_packages
   
    # For better dependency resolution, try to use mamba if available
    use_mamba = shutil.which("mamba") is not None
    
    if use_mamba:
        print("Using mamba for better dependency resolution")
        installer = "mamba"
    else:
        installer = "conda"
    
    # Prepare command to update existing environment
    cmd = [installer, 'install', '--name', env_name, '--yes',
           '-c', 'bioconda', '-c', 'conda-forge']
    
    # Add tools with version specifications if present
    cmd.extend(tools)
    
    try:
        # Run the update command
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print("Successfully updated environment '{}'".format(env_name))
        return True
    except subprocess.CalledProcessError as e:
        print("Error updating environment: {}".format(e))
        
        # If conda failed but we weren't using mamba, try mamba as fallback
        if not use_mamba and shutil.which("mamba") is not None:
            print("Trying again with mamba for better dependency resolution...")
            try:
                cmd[0] = "mamba"
                print(f"Running: {' '.join(cmd)}")
                subprocess.run(cmd, check=True)
                print("Successfully updated environment with mamba '{}'".format(env_name))
                return True
            except subprocess.CalledProcessError as e2:
                print("Error updating environment with mamba: {}".format(e2))
                return False
        
        return False

def smart_install_tools(tools_list):
    """
    Smart installation strategy:
    1. Try to install in metamag first (for tools with Python 3.9)
    2. Create separate environments only for known conflicts
    """
    # Tools that MUST have separate environments (known conflicts)
    MUST_SEPARATE = {
        "checkm2": ("checkm_env", "3.8"),   # Python version conflict
        "gtdbtk": ("gtdbtk", "3.9"),        # Complex dependencies
        "dbcan": ("dbcan", "3.9")           # Specific dependencies
    }
    
    # Separate tools into groups
    metamag_tools = []
    separate_env_tools = {}
    
    for tool in tools_list:
        tool_name = tool.replace('-', '_')
        
        # Check if this tool MUST be in separate environment
        if tool_name in MUST_SEPARATE:
            env_name, python_ver = MUST_SEPARATE[tool_name]
            if env_name not in separate_env_tools:
                separate_env_tools[env_name] = {'python': python_ver, 'tools': []}
            
            # Map to conda package name
            if tool in TOOL_CONDA_PACKAGE_MAP:
                conda_pkg = TOOL_CONDA_PACKAGE_MAP[tool]
                if conda_pkg:
                    separate_env_tools[env_name]['tools'].append(conda_pkg)
            else:
                separate_env_tools[env_name]['tools'].append(tool)
        else:
            # Try in metamag first
            if tool in TOOL_CONDA_PACKAGE_MAP:
                conda_pkg = TOOL_CONDA_PACKAGE_MAP[tool]
                if conda_pkg:
                    metamag_tools.append(conda_pkg)
            else:
                metamag_tools.append(tool)
    
    # Remove duplicates
    metamag_tools = list(set(metamag_tools))
    
    print("\n" + "="*60)
    print("INSTALLATION PLAN:")
    print("="*60)
    if metamag_tools:
        print(f"metamag (Python 3.9): {', '.join(metamag_tools)}")
    for env_name, env_info in separate_env_tools.items():
        tools = list(set(env_info['tools']))
        python_ver = env_info['python']
        print(f"{env_name} (Python {python_ver}): {', '.join(tools)}")
    print("="*60 + "\n")
    
    success = True
    
    # Step 1: Install most tools in metamag
    if metamag_tools:
        print(f"Installing {len(metamag_tools)} tools in 'metamag' environment...")
        if not update_existing_environment("metamag", metamag_tools):
            print("Warning: Some tools may have failed in metamag")
            success = False
    
    # Step 2: Install problematic tools in separate environments
    for env_name, env_info in separate_env_tools.items():
        tools = list(set(env_info['tools']))
        python_ver = env_info['python']
        
        env_exists = check_env_exists(env_name)
        
        if env_exists:
            print(f"\nUpdating '{env_name}' environment...")
            if not update_existing_environment(env_name, tools):
                success = False
        else:
            print(f"\nCreating '{env_name}' environment (Python {python_ver})...")
            if not create_environment_with_tools(env_name, python_ver, tools):
                success = False
    
    return success

def create_environment_with_tools(env_name, python_ver, tools):
    """Create a new conda environment with specified Python version and tools."""
    use_mamba = shutil.which("mamba") is not None
    installer = "mamba" if use_mamba else "conda"
    
    cmd = [installer, 'create', '--name', env_name, '--yes',
           '-c', 'bioconda', '-c', 'conda-forge',
           f'python={python_ver}'] + tools
    
    try:
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"Successfully created environment '{env_name}'")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating environment '{env_name}': {e}")
        return False


def verify_conda_environment(env_name, tool_paths):
    """Verify all tools in the conda environment are properly installed and check for missing pipeline tools."""
    print("\nVerifying tools in '{}' environment:".format(env_name))
    print("=" * 60)
    print("{:<20} {:<10} {:<30}".format("Tool", "Status", "Version/Notes"))
    print("=" * 60)
    
    all_success = True
    working_tools = []
    failed_tools = []
    
    # Get conda environment path
    conda_prefix = None
    if env_name != "system":  # Only try to get conda path for conda environments
        try:
            # Direct method to get conda env path
            result = subprocess.run(
                ['conda', 'info', '--envs', '--json'],
                capture_output=True,
                text=True,
                check=True
            )
            envs_info = json.loads(result.stdout)
            
            # Find the named environment in the list
            for env_path in envs_info.get("envs", []):
                if os.path.basename(env_path) == env_name:
                    conda_prefix = env_path
                    break
        except Exception as e:
            print(f"Warning: Could not determine conda prefix via conda info: {e}")
            
            # Fallback to environment variable if the direct method fails
            conda_prefix = os.environ.get('CONDA_PREFIX')
    
    # Collect all tools across all tool groups for pipeline requirements
    all_pipeline_tools = set()
    for group_tools in TOOL_GROUPS.values():
        for tool in group_tools:
            all_pipeline_tools.add(tool.replace('-', '_'))
    
    if conda_prefix:
        print(f"Using conda environment at: {conda_prefix}")
        bin_path = os.path.join(conda_prefix, "bin")
        
        # Map all executables
        found_tools = {}
        
        # Scan bin directory and map all executables
        if os.path.exists(bin_path):
            print(f"Scanning bin directory: {bin_path}")
            for executable in os.listdir(bin_path):
                exec_path = os.path.join(bin_path, executable)
                if os.access(exec_path, os.X_OK):
                    found_tools[executable] = exec_path
        else:
            print(f"Bin directory does not exist: {bin_path}")
            
        # Reverse mapping to also check tool names
        reverse_map = {executable: tool for tool, executable in TOOL_EXECUTABLE_MAP.items()}
            
        # Now check each tool with improved name resolution
        for tool_name, tool_path in sorted(tool_paths.items()):
            # Skip activation commands
            if tool_name.endswith("_activate"):
                continue
                
            original_tool_name = tool_name.replace('_', '-')
            found = False
            error_msg = "Not found in environment"
            
            # First, check if the tool path exists as specified
            if tool_path and os.path.exists(tool_path):
                found = True
            if tool_path:
                # Handle Python module case
                if tool_path.startswith("python -c"):
                    found = True
                elif os.path.exists(tool_path):
                    found = True
            # Check if the tool is directly in found_tools
            elif tool_name in found_tools:
                tool_path = found_tools[tool_name]
                found = True
            # Check original (hyphenated) name
            elif original_tool_name in found_tools:
                tool_path = found_tools[original_tool_name]
                found = True
            # Check mapped executable name
            elif tool_name in TOOL_EXECUTABLE_MAP and TOOL_EXECUTABLE_MAP[tool_name] in found_tools:
                tool_path = found_tools[TOOL_EXECUTABLE_MAP[tool_name]]
                found = True
            # Check original name with mapped executable
            elif original_tool_name in TOOL_EXECUTABLE_MAP and TOOL_EXECUTABLE_MAP[original_tool_name] in found_tools:
                tool_path = found_tools[TOOL_EXECUTABLE_MAP[original_tool_name]]
                found = True
            # Check reverse mapping (executable to tool)
            else:
                for exec_name, found_path in found_tools.items():
                    if exec_name in reverse_map and (reverse_map[exec_name] == tool_name or 
                                                  reverse_map[exec_name] == original_tool_name):
                        tool_path = found_path
                        found = True
                        break
                        
            # If we found the tool, verify it works
            if found and tool_path:
                success, message = verify_tool_installation(original_tool_name, tool_path)
                
                if success:
                    print("{:<20} {:<10} {:<30}".format(tool_name, "OK", message[:30]))
                    working_tools.append(tool_name)
                else:
                    print("{:<20} {:<10} {:<30}".format(tool_name, "FAILED", message[:30]))
                    all_success = False
                    failed_tools.append(tool_name)
            else:
                print("{:<20} {:<10} {:<30}".format(tool_name, "MISSING", error_msg))
                all_success = False
                failed_tools.append(tool_name)
    else:
        print("WARNING: Could not determine conda environment path.")
        for tool_name, tool_path in sorted(tool_paths.items()):
            # Skip activation commands
            if tool_name.endswith("_activate"):
                continue
                
            if not tool_path:
                print("{:<20} {:<10} {:<30}".format(tool_name, "MISSING", "Not found in PATH"))
                all_success = False
                failed_tools.append(tool_name)
                continue
                
            if not os.path.exists(tool_path):
                print("{:<20} {:<10} {:<30}".format(tool_name, "MISSING", "Path does not exist"))
                all_success = False
                failed_tools.append(tool_name)
                continue
                
            original_tool_name = tool_name.replace('_', '-')
            success, message = verify_tool_installation(original_tool_name, tool_path)
            
            if success:
                print("{:<20} {:<10} {:<30}".format(tool_name, "OK", message[:30]))
                working_tools.append(tool_name)
            else:
                print("{:<20} {:<10} {:<30}".format(tool_name, "FAILED", message[:30]))
                all_success = False
                failed_tools.append(tool_name)
    
    # Now check which pipeline tools are missing
    missing_pipeline_tools = []
    for tool in all_pipeline_tools:
        if tool not in working_tools and not tool.endswith("_db_dir"):
            missing_pipeline_tools.append(tool)
    
    print("=" * 60)
    print("Summary: {} tools working, {} tools failed".format(len(working_tools), len(failed_tools)))
    
    if failed_tools:
        print("\nFailed tools:")
        for tool in failed_tools:
            print("  - {}".format(tool))
            
        if conda_prefix:
            print("\nTroubleshooting tips:")
            print("  1. Activate the environment: conda activate {}".format(env_name))
            print("  2. Try reinstalling a tool: conda install -c bioconda {}".format(failed_tools[0]))
            print("  3. Check for conda environment conflicts")
    
    # Report missing pipeline tools
    if missing_pipeline_tools:
        print("\nMissing tools required for the pipeline:")
        # Group by tool category for better organization
        missing_by_group = {}
        for tool in missing_pipeline_tools:
            for group, group_tools in TOOL_GROUPS.items():
                if tool.replace('_', '-') in group_tools or tool in group_tools:
                    if group not in missing_by_group:
                        missing_by_group[group] = []
                    missing_by_group[group].append(tool)
                    break
        
        # Display missing tools by group
        for group, tools in sorted(missing_by_group.items()):
            print(f"\n{group.upper()} tools:")
            for tool in sorted(tools):
                print(f"  - {tool}")
        
        print("\nTo install missing tools, run:")
        print("  python setup_tools.py --update --tools " + " ".join(missing_pipeline_tools[:5]) + " ...")
        print("Or install all tools for specific steps:")
        print("  python setup_tools.py --update --steps qc preprocessing ...")
    
    return all_success, working_tools, failed_tools

def verify_tool_installation(tool_name, tool_path):
    """Verify a tool is properly installed and can be executed."""
   
    # Handle Python module case FIRST (before checking path existence)
    if tool_path and tool_path.startswith("python -c"):
        return True, f"{tool_name} module available"
    
    if not tool_path or not os.path.exists(tool_path):
        return False, "Tool not found at specified path"
    
    # Check if the file is executable
    if not os.access(tool_path, os.X_OK):
        return False, "Tool exists but is not executable"
    
    # Special handling for tools that show usage when run without arguments
    try:
        # Just check if the tool runs without arguments (for 1 second max)
        result = subprocess.run(
            [tool_path], 
            capture_output=True, 
            text=True,
            timeout=1
        )
        # If we get a usage message, that's a good sign the tool is working
        output = result.stdout + result.stderr
        if "usage:" in output.lower() or "usage :" in output.lower():
            return True, f"{os.path.basename(tool_path)} available"
    except Exception:
        # If this fails, continue with standard checks
        pass
    
    # Try to run the version check command if available
    if tool_name in TOOL_VERSION_COMMANDS:
        try:
            # Parse the version command
            cmd_parts = TOOL_VERSION_COMMANDS[tool_name].split()
            
            # Handle relative vs absolute paths
            if os.path.basename(tool_path) != cmd_parts[0]:
                # If tool_path doesn't match the command name, use tool_path directly
                cmd = [tool_path] + cmd_parts[1:]
            else:
                # Otherwise use the command as specified
                cmd = cmd_parts
                
            # Run the command with a short timeout
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True,
                timeout=5
            )
            
            # Check if the command was successful
            if result.returncode == 0:
                # Extract version information
                output = result.stdout if result.stdout else result.stderr
                version_info = output.strip().split('\n')[0]
                return True, version_info
            
            # Try with common version flags if the standard command failed
            for flag in ['--version', '-v', '--v', '-version', '-h', '--help']:
                try:
                    result = subprocess.run(
                        [tool_path, flag], 
                        capture_output=True, 
                        text=True,
                        timeout=3
                    )
                    if result.returncode == 0 or "version" in (result.stdout + result.stderr).lower():
                        output = result.stdout if result.stdout else result.stderr
                        return True, output.strip().split('\n')[0]
                except:
                    continue
                    
            # Tool exists but version command failed
            return True, f"{os.path.basename(tool_path)} available (no version)"
                
        except subprocess.TimeoutExpired:
            # Tool exists but took too long
            return True, f"{os.path.basename(tool_path)} available (timeout)"
        except Exception as e:
            # Something went wrong but the tool exists
            return True, f"{os.path.basename(tool_path)} available"
    
    # For tools without version commands, just return that they're available
    return True, f"{os.path.basename(tool_path)} available"

def main():
    parser = argparse.ArgumentParser(description="MetaMAG Setup Utility")
    parser.add_argument("--steps", nargs="+", help="Specific pipeline steps to install tools for")
    parser.add_argument("--tools", nargs="+", help="Specific tools to install")
    parser.add_argument("--all", action="store_true", help="Install all tools")
    parser.add_argument("--use-existing", action="store_true", help="Use existing tools instead of installing")
    parser.add_argument("--update", action="store_true", help="Update existing environment instead of creating a new one")
    parser.add_argument("--force-recreate", action="store_true", help="Force recreation of environment even if it exists")
    parser.add_argument("--verify", action="store_true", help="Verify tool installation and functionality")
    parser.add_argument("--interactive", action="store_true", help="Run in interactive mode (with user prompts)")
    parser.add_argument("--verify-after", action="store_true", help="Verify tool installation after setup/installation")
    parser.add_argument("--generate-config-on-error", action="store_true", help="Generate a config from PATH if conda installation fails")
    
    args = parser.parse_args()

    print("MetaMAG Setup Utility")
    print("=======================")
    
    # Default behavior is now non-interactive unless specifically requested
    interactive = args.interactive
    verify_after = args.verify_after
    
    # If just verifying tools, do that and exit
    if args.verify:
        # Try to find the config file
        config_path = find_config_file()
        
        if config_path:
            try:
                # Add the directory containing the config file to Python path
                sys.path.insert(0, os.path.dirname(os.path.abspath(config_path)))
                
                # Dynamically import the config
                import importlib.util
                spec = importlib.util.spec_from_file_location("config", config_path)
                config_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(config_module)
                

                # Use the imported config
                config = config_module.config
                
                # Get tools from both locations
                direct_tools = config.get("tools", {})
                detected_tools = config.get("tool_configs", {})
                
                # Merge them - prefer detected_tools (from all environments)
                tool_paths = {}
                
                # First add direct tools
                for tool_name, tool_path in direct_tools.items():
                    if tool_path:  # Only add if not empty
                        tool_paths[tool_name] = tool_path
                
                # Then add/override with detected tools from all environments
                for tool_name, tool_info in detected_tools.items():
                    if tool_info:
                        if isinstance(tool_info, dict) and "path" in tool_info:
                            tool_paths[tool_name] = tool_info["path"]
                        elif isinstance(tool_info, str):
                            tool_paths[tool_name] = tool_info


                
                # Determine environment name
                env_name = "metamag"  # Default environment name
                conda_env_path = config.get("environment", {}).get("conda_env", "")
                if conda_env_path and "activate" in conda_env_path:
                    env_name = conda_env_path.split()[-1]
                
                verify_conda_environment(env_name, tool_paths)
            
            except Exception as e:
                print(f"Error loading config file: {e}")
                print("Detailed debugging information:")
                print(f"Config file path: {config_path}")
                print("Current working directory:", os.getcwd())
                print("Python path:", sys.path)
        
        else:
            # If no config file found, provide detailed guidance
            print("ERROR: No config.py found.")
            print("\nPossible solutions:")
            print("1. Ensure you're in the correct directory")
            print("2. Run setup_tools.py with --use-existing to generate config")
            print("3. Manually create config.py in the current directory")
            
            # Print current directory contents for debugging
            print("\nCurrent directory contents:")
            try:
                print(os.listdir('.'))
            except Exception as e:
                print(f"Could not list directory contents: {e}")
        
        return   

    # Check if conda is available if not using existing tools
    if not is_conda_available() and not args.use_existing:
        print("\nERROR: Conda is not available in your PATH.")
        print("Conda is required for automatic tool installation.")
        print("\nYou have the following options:")
        print("  1. Install Conda first (https://docs.conda.io/en/latest/miniconda.html)")
        print("  2. Use --use-existing if you already have the tools installed elsewhere")
        print("  3. Manually create and edit config.py with your tool paths")
        
        # If interactive, ask if user wants to proceed with --use-existing
        if interactive:
            choice = input("\nDo you want to proceed with the --use-existing option to scan for tools in your PATH? (yes/no): ")
            if choice.lower() == 'yes':
                args.use_existing = True
            else:
                return
        else:
            # In non-interactive mode, default to using existing tools
            print("\nAutomatically proceeding with --use-existing in non-interactive mode")
            args.use_existing = True
    
    if args.use_existing:
        # Use existing tools
        print("\nUsing existing tools in PATH...")
        generate_path_config()
        print("\nSetup complete! A config file has been generated based on tools found in your PATH.")
        print("If some tools weren't found, edit the config.py file to provide correct paths.")
        
        # If interactive and not explicitly set, ask if user wants to verify tool installation
        if interactive and not verify_after:
            choice = input("\nDo you want to verify tool installation now? (yes/no): ")
            verify_after = choice.lower() == 'yes'
        
        if verify_after:
            # Load the newly generated config
            try:
                sys.path.append(os.path.dirname(os.path.abspath(__file__)))
                import importlib
                if 'config' in sys.modules:
                    importlib.reload(sys.modules['config'])
                else:
                    from config import config
                tool_paths = config.get("tools", {})
                verify_conda_environment("system", tool_paths)
            except Exception as e:
                print("Error during verification: {}".format(e))
        
        return

    # Check for conflicting arguments
    if args.update and args.force_recreate:
        print("Error: --update and --force-recreate cannot be used together")
        return

    # Check if metamag environment already exists
    env_exists = False
    if is_conda_available():
        env_exists = check_env_exists("metamag")
    
    # Handle existing environment based on arguments
    if env_exists:
        if args.force_recreate:
            print("Environment 'metamag' exists. Removing and recreating...")
            subprocess.run(['conda', 'env', 'remove', '--name', 'metamag', '--yes'], check=True)
            env_exists = False
        elif not args.update:
            print("Environment 'metamag' already exists.")
            print("Use --update to add tools to the existing environment")
            print("Or use --force-recreate to remove and recreate the environment")
            return
    elif args.update:
        print("Environment 'metamag' does not exist. Cannot update.")
        print("Remove the --update flag to create a new environment.")
        return

    if args.all:
        # Get all tools
        all_tools = set()
        for group in TOOL_GROUPS.values():
            all_tools.update(group)
        
        if env_exists and args.update:
            # Update existing environment with all tools
            print("Updating with all tools...")
            update_existing_environment("metamag", sorted(all_tools))
            generate_conda_config()
        else:
            # Install all tools in new environment
            print("Installing all tools via Conda...")
            install_tools(list(TOOL_GROUPS.keys()), args.generate_config_on_error)
      
    elif args.steps:
        # Determine required tools based on requested steps
        tool_groups = set()
        for step in args.steps:
            # First check if it's a pipeline step
            if step in STEP_TO_TOOLS:
                tool_groups.update(STEP_TO_TOOLS[step])
            # Then check if it's a direct tool group name
            elif step in TOOL_GROUPS:
                tool_groups.add(step)
            else:
                print("Warning: Unknown step '{}'. Skipping.".format(step))
                continue
        
        print("Installing tools for steps: {}".format(', '.join(args.steps)))
        print("Required tool groups: {}".format(', '.join(tool_groups)))
        
        # Collect all tools from the specified groups
        tools_to_install = set()
        for group in tool_groups:
            if group in TOOL_GROUPS:
                tools_to_install.update(TOOL_GROUPS[group])
        
        if args.update:
            # Use smart installation strategy
            smart_install_tools(sorted(tools_to_install))
            generate_conda_config()
        else:
            # Install tools in new environment
            install_tools(tool_groups, args.generate_config_on_error)


    elif args.tools:
        # Install specific tools
        if args.update:
            # Use smart installation strategy
            smart_install_tools(args.tools)
            generate_conda_config()
        else:
            # Install tools in new environment
            install_specific_tools(args.tools, args.generate_config_on_error)
    
    # If no specific action was specified and we're in interactive mode, show the menu
    elif interactive:
        show_interactive_menu(env_exists, args.update)
    else:
        print("No specific action specified. Use --all, --steps, --tools, or --use-existing.")
        return
    
    # Check if verification is requested after installation
    if verify_after:
        # Load the newly generated config
        try:
            sys.path.append(os.path.dirname(os.path.abspath(__file__)))
            import importlib
            if 'config' in sys.modules:
                importlib.reload(sys.modules['config'])
            else:
                from config import config
            tool_paths = config.get("tools", {})
            verify_conda_environment("metamag", tool_paths)
        except Exception as e:
            print("Error during verification: {}".format(e))

def show_interactive_menu(env_exists=False, update=False):
    """Show an interactive menu for tool selection."""
    print("\nPlease select an installation option:")
    
    if env_exists:
        print("1. Update existing environment with all tools")
        print("2. Update existing environment with tools for specific pipeline steps")
        print("3. Update existing environment with specific tools")
        print("4. Use existing tools (no installation)")
        print("5. Force recreate environment with all tools (warning: will delete existing environment)")
        print("6. Verify tool installation")
    else:
        print("1. Install all tools")
        print("2. Install tools for specific pipeline steps")
        print("3. Install specific tools")
        print("4. Use existing tools (no installation)")
        print("5. Verify tool installation")
    
    choice = input("\nEnter your choice: ")
    
    if choice == "1":
        if env_exists and update:
            # Update with all tools
            all_tools = set()
            for group in TOOL_GROUPS.values():
                all_tools.update(group)
            update_existing_environment("metamag", sorted(all_tools))
            generate_conda_config()
        else:
            # Install all tools
            install_tools(list(TOOL_GROUPS.keys()), False)
    
    elif choice == "2":
        available_steps = sorted(STEP_TO_TOOLS.keys())
        print("\nAvailable pipeline steps:")
        for i, step in enumerate(available_steps, 1):
            print("{}. {}".format(i, step))
        
        step_indices = input("\nEnter the numbers of steps you want to install tools for (comma-separated): ")
        try:
            selected_indices = [int(idx.strip()) - 1 for idx in step_indices.split(",")]
            selected_steps = [available_steps[idx] for idx in selected_indices if 0 <= idx < len(available_steps)]
            
            tool_groups = set()
            for step in selected_steps:
                tool_groups.update(STEP_TO_TOOLS[step])
            
            print("Selected steps: {}".format(', '.join(selected_steps)))
            print("Required tool groups: {}".format(', '.join(tool_groups)))
            
            # Collect all tools from the specified groups
            tools_to_install = set()
            for group in tool_groups:
                if group in TOOL_GROUPS:
                    tools_to_install.update(TOOL_GROUPS[group])
            
            if env_exists and update:
                # Update existing environment with selected tools
                update_existing_environment("metamag", sorted(tools_to_install))
                generate_conda_config()
            else:
                # Install tools in new environment
                install_tools(tool_groups, False)
                
        except (ValueError, IndexError):
            print("Invalid selection. Please run the setup script again.")
    
    elif choice == "3":
        all_tools = []
        for group, tools in TOOL_GROUPS.items():
            all_tools.extend(tools)
        
        print("\nAvailable tools:")
        for i, tool in enumerate(sorted(all_tools), 1):
            print("{}. {}".format(i, tool))
        
        tool_indices = input("\nEnter the numbers of tools you want to install (comma-separated): ")
        try:
            selected_indices = [int(idx.strip()) - 1 for idx in tool_indices.split(",")]
            selected_tools = [sorted(all_tools)[idx] for idx in selected_indices if 0 <= idx < len(all_tools)]
            
            if env_exists and update:
                # Update existing environment with selected tools
                update_existing_environment("metamag", selected_tools)
                generate_conda_config()
            else:
                # Install tools in new environment
                install_specific_tools(selected_tools, False)
                
        except (ValueError, IndexError):
            print("Invalid selection. Please run the setup script again.")
    
    elif choice == "4":
        generate_path_config()
        print("\nSetup complete! A config file has been generated based on tools found in your PATH.")
        print("If some tools weren't found, edit the config.py file to provide correct paths.")
    
    elif (choice == "5" and env_exists) or (choice == "6" and env_exists):
        if choice == "5":
            print("Warning: This will delete your existing 'metamag' environment.")
            confirm = input("Type 'yes' to confirm: ")
            if confirm.lower() == 'yes':
                subprocess.run(['conda', 'env', 'remove', '--name', 'metamag', '--yes'], check=True)
                install_tools(list(TOOL_GROUPS.keys()), False)
            else:
                print("Operation cancelled.")
        else:  # choice == "6"
            # Load the config file
            try:
                sys.path.append(os.path.dirname(os.path.abspath(__file__)))
                import importlib
                if 'config' in sys.modules:
                    importlib.reload(sys.modules['config'])
                else:
                    from config import config
                tool_paths = config.get("tools", {})
                verify_conda_environment("metamag", tool_paths)
            except Exception as e:
                print("Error loading config.py: {}".format(e))
                print("Please run setup_tools.py first to generate a config file.")
    
    elif choice == "5" and not env_exists:
        # Verify tool installation
        if os.path.exists("config.py"):
            try:
                sys.path.append(os.path.dirname(os.path.abspath(__file__)))
                from config import config
                tool_paths = config.get("tools", {})
                verify_conda_environment("system", tool_paths)
            except Exception as e:
                print("Error loading config.py: {}".format(e))
                print("Please run setup_tools.py first to generate a config file.")
        else:
            print("Error: config.py not found. Please run setup_tools.py first.")
    
    else:
        print("Invalid choice. Please run setup.py again.")

def install_tools(tool_groups, generate_config_on_error=False):
    """Install tools from specified tool groups."""
    # Create a custom environment file
    create_custom_environment_file(tool_groups)
    
    # Install via Conda
    print("Installing selected tools via Conda...")
    try:
        subprocess.run(["conda", "env", "create", "-f", "custom_environment.yml"], check=True)
        
        # Generate config with Conda paths
        generate_conda_config()
        
        print("\nSetup complete! Activate the environment with: conda activate metamag")
        print("Then run the pipeline as described in the README.md")
    except subprocess.CalledProcessError as e:
        print("\nError during environment creation: {}".format(e))
        print("You may need to check your Conda installation or network connection.")
        
        # Generate a config file anyway if requested
        if generate_config_on_error:
            print("\nGenerating a config file based on existing tools...")
            generate_path_config()
            print("\nConfig file generated based on tools found in PATH.")

def install_specific_tools(tools, generate_config_on_error=False):
    """Install specific tools.
    
    Supports tool version specifications like 'gtdbtk=2.3.2'
    """
    print("Installing specific tools: {}".format(', '.join(tools)))
    
    # Create a custom environment file with only the specified tools
    with open("custom_environment.yml", "w") as f:
        f.write("name: metamag\n")
        f.write("channels:\n")
        f.write("  - conda-forge\n")
        f.write("  - bioconda\n")
        f.write("  - defaults\n")
        f.write("dependencies:\n")
        f.write("  - python=3.9\n")
        
        # Add requested tools (preserving version specifications)
        for tool in tools:
            if tool not in ["idba_fq2fa", "DAS_Tool_Fasta2Contig", "eggnog_db_dir", "dbcan_db_dir"]:
                f.write("  - {}\n".format(tool))
        
        # Always include essential Python packages
        f.write("  - pandas\n")
        f.write("  - numpy\n")
        f.write("  - pyyaml\n")
        f.write("  - biopython\n")
    
    # Try to install with conda first
    try:
        subprocess.run(["conda", "env", "create", "-f", "custom_environment.yml"], check=True)
        
        # Generate config with Conda paths
        generate_conda_config()
        
        print("\nSetup complete! Activate the environment with: conda activate metamag")
        print("Then run the pipeline as described in the README.md")
    except subprocess.CalledProcessError as e:
        print("\nError during environment creation with conda: {}".format(e))
        
        # Try mamba if available
        if shutil.which("mamba") is not None:
            print("Trying with mamba for better dependency resolution...")
            try:
                subprocess.run(["mamba", "env", "create", "-f", "custom_environment.yml"], check=True)
                
                # Generate config with Conda paths
                generate_conda_config()
                
                print("\nSetup complete! Activate the environment with: conda activate metamag")
                print("Then run the pipeline as described in the README.md")
                return
            except subprocess.CalledProcessError as e2:
                print("Error during environment creation with mamba: {}".format(e2))
        
        print("You may need to check your Conda installation or network connection.")
        
        # Generate a config file anyway if requested
        if generate_config_on_error:
            print("\nGenerating a config file based on existing tools...")
            generate_path_config()
            print("\nConfig file generated based on tools found in PATH.")

def generate_conda_config():
    """Generate a config file for Conda-installed tools that only includes tools that actually exist."""
    # Get the Conda environment path
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        # Try to get it from conda info
        try:
            result = subprocess.run(
                ['conda', 'info', '--json'],
                capture_output=True,
                text=True,
                check=True
            )
            conda_info = json.loads(result.stdout)
            
            # Find metamag environment path
            for env_path in conda_info.get("envs", []):
                if os.path.basename(env_path) == "metamag":
                    conda_prefix = env_path
                    break
                    
            # If not found, use base conda path
            if not conda_prefix:
                conda_prefix = conda_info.get("conda_prefix", "")
                
        except Exception as e:
            print("Error getting Conda info: {}".format(e))
            conda_prefix = ""
    
    if not conda_prefix:
        print("Error: Could not determine conda prefix. Cannot generate config.")
        return
        
    # First, scan the environment and find all available tools
    bin_path = os.path.join(conda_prefix, "bin")
    available_tools = {}
    
    if os.path.exists(bin_path):
        # Map of expected executable names to their tool names
        exec_to_tool = {executable: tool for tool, executable in TOOL_EXECUTABLE_MAP.items()}
        
        # Scan the bin directory
        for executable in os.listdir(bin_path):
            exec_path = os.path.join(bin_path, executable)
            
            # Skip if not an executable file
            if not os.path.isfile(exec_path) or not os.access(exec_path, os.X_OK):
                continue
                
            # Check if this is a known tool
            tool_name = exec_to_tool.get(executable, executable)
            
            # Store the path to the tool
            available_tools[tool_name.replace('-', '_')] = exec_path
    
    # Also check for Python modules that aren't in bin/
   
    for tool_name, exec_name in TOOL_EXECUTABLE_MAP.items():
        if exec_name and exec_name.startswith("python -c"):
            # This is a Python module, try to import it
            print(f"Checking Python module: {tool_name} with command: {exec_name}")
            try:
                # Extract module name from the import statement
                import_cmd = exec_name.replace("python -c '", "").replace("'", "")
                print(f"  Import command: {import_cmd}")
                result = subprocess.run(
                    ['python', '-c', import_cmd],
                    capture_output=True,
                    text=True,
                    timeout=2
                )
                print(f"  Return code: {result.returncode}")
                if result.returncode == 0:
                    # Module is available
                    available_tools[tool_name.replace('-', '_')] = exec_name
                    print(f"  SUCCESS: Added {tool_name} to available_tools")
                else:
                    print(f"  FAILED: stderr: {result.stderr}")
            except Exception as e:
                print(f"  EXCEPTION: {e}")

    config_path = "MetaMAG/config.py"
    if not os.path.exists("MetaMAG"):
        # If MetaMAG directory doesn't exist, create it
        os.makedirs("MetaMAG", exist_ok=True)
        # Also create __init__.py if it doesn't exist
        if not os.path.exists("MetaMAG/__init__.py"):
            with open("MetaMAG/__init__.py", "w") as f:
                f.write('__version__ = "1.0.0"\n')
 
    # Write the config file with only the tools that exist
    with open(config_path, "w") as f:
        f.write("import os\n")
        f.write("import shutil\n")
        f.write("import yaml\n")
        f.write("import subprocess\n\n")
        
        f.write("def find_tool(tool_name, default_path=None):\n")
        f.write("    \"\"\"Try to find a tool in PATH, otherwise return the default path.\"\"\"\n")
        f.write("    path_result = shutil.which(tool_name)\n")
        f.write("    if path_result:\n")
        f.write("        return path_result\n")
        f.write("    return default_path\n\n")
        
        f.write("def is_conda_available():\n")
        f.write("    \"\"\"Check if conda is available in the system.\"\"\"\n")
        f.write("    return shutil.which('conda') is not None\n\n")
        
        f.write("def detect_tool_locations():\n")
        f.write("    \"\"\"Detect tool locations based on system paths and conda environments.\"\"\"\n")
        f.write("    tools = {}\n")
        f.write("    \n")
        f.write("    # Dictionary of tool names and their possible conda environments\n")
        f.write("    conda_environments = {\n")
        
        # Add tool-to-env mappings
        for tool, env in TOOL_CONDA_ENVS.items():
            f.write(f"        \"{tool.replace('-', '_')}\": \"{env}\",\n")
        
        f.write("    }\n")
        f.write("    \n")
        f.write("    # Check if conda is available\n")
        f.write("    conda_available = is_conda_available()\n")
        f.write("    \n")
        f.write("    # Try to find conda root path\n")
        f.write("    conda_root = None\n")
        f.write("    if conda_available:\n")
        f.write("        try:\n")
        f.write("            result = subprocess.run(['conda', 'info', '--json'], \n")
        f.write("                                   capture_output=True, text=True, check=True)\n")
        f.write("            import json\n")
        f.write("            conda_info = json.loads(result.stdout)\n")
        f.write("            conda_root = conda_info.get('root_prefix')\n")
        f.write("        except Exception as e:\n")
        f.write("            print(f\"Warning: Conda available but couldn't get info: {e}\")\n")
        f.write("    \n")
        f.write("    # Map of expected tool names to their executable names\n")
        f.write("    tool_to_exec = {\n")
        
        # Add tool-to-executable mappings
        for tool, executable in TOOL_EXECUTABLE_MAP.items():
            f.write(f"        \"{tool.replace('-', '_')}\": \"{executable}\",\n")
        
        f.write("    }\n")
        f.write("    \n")
        f.write("    # Process each tool\n")
        f.write("    for tool_name, env_name in conda_environments.items():\n")
        f.write("        # Get the executable name if it differs\n")
        f.write("        exec_name = tool_to_exec.get(tool_name, tool_name)\n")
        f.write("        \n")
        f.write("        # First try to find the tool in PATH\n")
        f.write("        direct_path = shutil.which(exec_name)\n")
        f.write("        if direct_path:\n")
        f.write("            tools[tool_name] = {\n")
        f.write("                \"path\": direct_path,\n")
        f.write("                \"requires_activation\": False\n")
        f.write("            }\n")
        f.write("            continue\n")
        f.write("        \n")
        f.write("        # If not in PATH and conda is available, check in conda environment\n")
        f.write("        if conda_available and conda_root:\n")
        f.write("            # Check if environment exists\n")
        f.write("            env_path = os.path.join(conda_root, \"envs\", env_name)\n")
        f.write("            if os.path.exists(env_path):\n")
        f.write("                # Check if tool exists in this environment\n")
        f.write("                tool_in_env = os.path.join(env_path, \"bin\", exec_name)\n")
        f.write("                if os.path.exists(tool_in_env):\n")
        f.write("                    tools[tool_name] = {\n")
        f.write("                        \"path\": tool_in_env,\n")
        f.write("                        \"requires_activation\": True,\n")
        f.write("                        \"env_name\": env_name,\n")
        f.write("                        \"activate_cmd\": f\"conda activate {env_name}\"\n")
        f.write("                    }\n")
        f.write("                    continue\n")
        f.write("        \n")
        f.write("        # Tool not found in PATH or conda\n")
        f.write("        tools[tool_name] = None\n")
        f.write("    \n")
        f.write("    return tools\n\n")
        
        # Start config dictionary
        f.write("# Base config - only includes tools that actually exist\n")
        f.write("config = {\n")
        f.write("    \"environment\": {\n")
        perl_lib_path = detect_perl_lib_path()
        f.write(f"        \"PERL5LIB\": \"{perl_lib_path}\",\n")
        
        conda_path = shutil.which("conda")
        if conda_path:
            conda_base = os.path.dirname(os.path.dirname(conda_path))  # Go up two levels
            conda_activate_path = os.path.join(conda_base, "bin", "activate")
            f.write(f"        \"conda_env\": \"source {conda_activate_path} metamag\",\n")
        else:
            f.write("        \"conda_env\": \"source activate metamag\",  # Update this path\n")
      
        f.write("    },\n")
        f.write("    \"tools\": {\n")
        
        # Group tools by category for better organization in the config
        tools_by_group = {}
        for group, tools in TOOL_GROUPS.items():
            tools_by_group[group] = []
            for tool in tools:
                tool_name = tool.replace('-', '_')
                if tool_name in available_tools:
                    tools_by_group[group].append((tool_name, available_tools[tool_name]))
        
        # Add each tool with its path, but ONLY if it actually exists
        for group, tools in tools_by_group.items():
            if tools:  # Only add group comment if we have tools in this group
                f.write(f"        # {group} tools\n")
                for tool_name, tool_path in tools:
                    f.write("        \"{}\": \"{}\",\n".format(tool_name, tool_path))
                
                # Add activation commands for known environments that have tools
                if group in ["binning", "evaluation", "taxonomy", "annotation"]:
                    env_names = set()
                    for tool in TOOL_GROUPS[group]:
                        tool_name = tool.replace('-', '_')
                        env_name = TOOL_CONDA_ENVS.get(tool)
                        if tool_name in available_tools and env_name and env_name != "metamag":
                            env_names.add(env_name)
                    
                    for env_name in env_names:
                        activate_path = os.path.join(os.path.dirname(conda_prefix), "bin", "activate")
                        if os.path.exists(activate_path):
                            f.write(f"        \"{env_name}_activate\": \"{activate_path} {env_name}\",\n")
                
                f.write("\n")
        
        # Check for available database directories
        if "dbcan" in available_tools:
            dbcan_env = TOOL_CONDA_ENVS.get("dbcan")
            if dbcan_env:
                dbcan_db_dir = os.path.join(os.path.dirname(conda_prefix), "envs", dbcan_env, "databases")
                if os.path.exists(dbcan_db_dir):
                    f.write("        # Database directories\n")
                    f.write(f"        \"dbcan_db_dir\": \"{dbcan_db_dir}\",\n\n")
        
        f.write("    },\n")
        f.write("    \"tool_configs\": {}\n")
        f.write("}\n\n")
        
        # Add auto-detection code
        f.write("# Try to auto-detect tools\n")
        f.write("try:\n")
        f.write("    detected_tools = detect_tool_locations()\n")
        f.write("    config[\"tool_configs\"] = detected_tools\n")
        f.write("except Exception as e:\n")
        f.write("    print(f\"Warning: Could not auto-detect tools: {e}\")\n")
        f.write("    # Continue with default paths\n\n")
        
        # Add user config support
        f.write("# Allow for user config to override default settings\n")
        f.write("user_config_path = os.environ.get('METAMAG_CONFIG', os.path.expanduser('~/.metamag/config.yaml'))\n")
        f.write("if os.path.exists(user_config_path):\n")
        f.write("    try:\n")
        f.write("        with open(user_config_path, 'r') as f:\n")
        f.write("            user_config = yaml.safe_load(f)\n")
        f.write("        \n")
        f.write("        # Update config with user settings\n")
        f.write("        if user_config.get('environment'):\n")
        f.write("            config['environment'].update(user_config['environment'])\n")
        f.write("        if user_config.get('tools'):\n")
        f.write("            config['tools'].update(user_config['tools'])\n")
        f.write("        if user_config.get('tool_configs'):\n")
        f.write("            config['tool_configs'].update(user_config['tool_configs'])\n")
        f.write("            \n")
        f.write("    except Exception as e:\n")
        f.write("        print(f\"Warning: Could not load user config from {user_config_path}: {e}\")\n")
    
    # Add helper functions for backward compatibility
    with open(config_path, "a") as f:  # Use append mode
        f.write("\n# Helper functions for backward compatibility and conda support\n")
        f.write("def get_tool_path(tool_name):\n")
        f.write("    \"\"\"\n")
        f.write("    Get the correct tool path from either direct tools config or detected tools.\n")
        f.write("    This function provides backward compatibility with old config structure.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        tool_name (str): Name of the tool to get path for\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        str: Path to the tool executable, or None if not found\n")
        f.write("    \"\"\"\n")
        f.write("    # First check if tool is in the direct tools config\n")
        f.write("    if tool_name in config[\"tools\"]:\n")
        f.write("        return config[\"tools\"][tool_name]\n")
        f.write("    \n")
        f.write("    # Then check if tool is in detected tool configs\n")
        f.write("    if tool_name in config.get(\"tool_configs\", {}):\n")
        f.write("        tool_config = config[\"tool_configs\"][tool_name]\n")
        f.write("        if tool_config and isinstance(tool_config, dict):\n")
        f.write("            return tool_config.get(\"path\")\n")
        f.write("        return tool_config\n")
        f.write("    \n")
        f.write("    # Tool not found\n")
        f.write("    print(f\"Warning: Tool '{tool_name}' not found in configuration\")\n")
        f.write("    return None\n\n")
        
        f.write("def get_tool_command(tool_name):\n")
        f.write("    \"\"\"\n")
        f.write("    Get the complete command to run a tool, including conda activation if needed.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        tool_name (str): Name of the tool\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        tuple: (command_prefix, tool_path) where command_prefix includes conda activation if needed\n")
        f.write("    \"\"\"\n")
        f.write("    # First check direct tools config\n")
        f.write("    if tool_name in config[\"tools\"]:\n")
        f.write("        return \"\", config[\"tools\"][tool_name]\n")
        f.write("    \n")
        f.write("    # Check detected tool configs\n")
        f.write("    if tool_name in config.get(\"tool_configs\", {}):\n")
        f.write("        tool_config = config[\"tool_configs\"][tool_name]\n")
        f.write("        if tool_config and isinstance(tool_config, dict):\n")
        f.write("            if tool_config.get(\"requires_activation\", False):\n")
        f.write("                activate_cmd = tool_config.get(\"activate_cmd\", \"\")\n")
        f.write("                return f\"source {activate_cmd} && \", tool_config.get(\"path\")\n")
        f.write("            else:\n")
        f.write("                return \"\", tool_config.get(\"path\")\n")
        f.write("        return \"\", tool_config\n")
        f.write("    \n")
        f.write("    print(f\"Warning: Tool '{tool_name}' not found in configuration\")\n")
        f.write("    return \"\", None\n")
    
    # Print summary (keep your existing print statements)
    print("\nGenerated config.py with {} available tools.".format(len(available_tools)))
    print("Only tools that actually exist have been added to the configuration.")

def generate_path_config():
    """Generate a config file with all pipeline tools (with empty paths for unavailable tools)."""
    # Collect all available tools from PATH
    available_tools = {}
    
    for group in TOOL_GROUPS.values():
        for tool in group:
            # Handle tool name differences
            tool_name = tool.replace('-', '_')
            search_name = TOOL_EXECUTABLE_MAP.get(tool, tool)
            
            path = shutil.which(search_name)
            if path:
                available_tools[tool_name] = path
    config_path = "MetaMAG/config.py"
    if not os.path.exists("MetaMAG"):
        # If MetaMAG directory doesn't exist, create it
        os.makedirs("MetaMAG", exist_ok=True)
        # Also create __init__.py if it doesn't exist
        if not os.path.exists("MetaMAG/__init__.py"):
            with open("MetaMAG/__init__.py", "w") as f:
                f.write('__version__ = "1.0.0"\n')

    with open(config_path, "w") as f:
        f.write("import os\n")
        f.write("import shutil\n")
        f.write("import yaml\n")
        f.write("import subprocess\n\n")
        
        f.write("def find_tool(tool_name, default_path=None):\n")
        f.write("    \"\"\"Try to find a tool in PATH, otherwise return the default path.\"\"\"\n")
        f.write("    path_result = shutil.which(tool_name)\n")
        f.write("    if path_result:\n")
        f.write("        return path_result\n")
        f.write("    return default_path\n\n")
        
        f.write("def is_conda_available():\n")
        f.write("    \"\"\"Check if conda is available in the system.\"\"\"\n")
        f.write("    return shutil.which('conda') is not None\n\n")
        
        f.write("def detect_tool_locations():\n")
        f.write("    \"\"\"Detect tool locations based on system paths and conda environments.\"\"\"\n")
        f.write("    tools = {}\n")
        f.write("    \n")
        f.write("    # Dictionary of tool names and their possible conda environments\n")
        f.write("    conda_environments = {\n")
        
        # Add tool-to-env mappings
        for tool, env in TOOL_CONDA_ENVS.items():
            f.write(f"        \"{tool.replace('-', '_')}\": \"{env}\",\n")
        
        f.write("    }\n")
        f.write("    \n")
        f.write("    # Map of expected tool names to their executable names\n")
        f.write("    tool_to_exec = {\n")
        
        # Add tool-to-executable mappings
        for tool, executable in TOOL_EXECUTABLE_MAP.items():
            f.write(f"        \"{tool.replace('-', '_')}\": \"{executable}\",\n")
        
        f.write("    }\n")
        f.write("    \n")
        f.write("    # Check if conda is available\n")
        f.write("    conda_available = is_conda_available()\n")
        f.write("    \n")
        f.write("    # Try to find conda root path\n")
        f.write("    conda_root = None\n")
        f.write("    if conda_available:\n")
        f.write("        try:\n")
        f.write("            result = subprocess.run(['conda', 'info', '--json'], \n")
        f.write("                                   capture_output=True, text=True, check=True)\n")
        f.write("            import json\n")
        f.write("            conda_info = json.loads(result.stdout)\n")
        f.write("            conda_root = conda_info.get('root_prefix')\n")
        f.write("        except Exception as e:\n")
        f.write("            print(f\"Warning: Conda available but couldn't get info: {e}\")\n")
        f.write("    \n")
        f.write("    # Process each tool\n")
        f.write("    for tool_name, env_name in conda_environments.items():\n")
        f.write("        # Get the executable name if it differs\n")
        f.write("        exec_name = tool_to_exec.get(tool_name, tool_name)\n")
        f.write("        \n")
        f.write("        # First try to find the tool in PATH\n")
        f.write("        direct_path = shutil.which(exec_name)\n")
        f.write("        if direct_path:\n")
        f.write("            tools[tool_name] = {\n")
        f.write("                \"path\": direct_path,\n")
        f.write("                \"requires_activation\": False\n")
        f.write("            }\n")
        f.write("            continue\n")
        f.write("        \n")
        f.write("        # If not in PATH and conda is available, check in conda environment\n")
        f.write("        if conda_available and conda_root:\n")
        f.write("            # Check if environment exists\n")
        f.write("            env_path = os.path.join(conda_root, \"envs\", env_name)\n")
        f.write("            if os.path.exists(env_path):\n")
        f.write("                # Check if tool exists in this environment\n")
        f.write("                tool_in_env = os.path.join(env_path, \"bin\", exec_name)\n")
        f.write("                if os.path.exists(tool_in_env):\n")
        f.write("                    tools[tool_name] = {\n")
        f.write("                        \"path\": tool_in_env,\n")
        f.write("                        \"requires_activation\": True,\n")
        f.write("                        \"env_name\": env_name,\n")
        f.write("                        \"activate_cmd\": f\"conda activate {env_name}\"\n")
        f.write("                    }\n")
        f.write("                    continue\n")
        f.write("        \n")
        f.write("        # Tool not found in PATH or conda\n")
        f.write("        tools[tool_name] = None\n")
        f.write("    \n")
        f.write("    return tools\n\n")
        
        # Start config dictionary
        f.write("# Base config - includes all tools needed for the pipeline\n")
        f.write("config = {\n")
        f.write("    \"environment\": {\n")
        perl_lib_path = detect_perl_lib_path()
        f.write(f"        \"PERL5LIB\": \"{perl_lib_path}\",\n")
       
        conda_path = shutil.which("conda")
        if conda_path:
            conda_base = os.path.dirname(os.path.dirname(conda_path))  # Go up two levels
            conda_activate_path = os.path.join(conda_base, "bin", "activate")
            f.write(f"        \"conda_env\": \"source {conda_activate_path} metamag\",\n")
        else:
            f.write("        \"conda_env\": \"source activate metamag\",  # Update this path\n")

        
        f.write("    },\n")
        f.write("    \"tools\": {\n")
        
        # List all tools by group, populated with paths when available
        for group, tools in TOOL_GROUPS.items():
            f.write(f"        # {group} tools\n")
            for tool in tools:
                tool_name = tool.replace('-', '_')
                tool_path = available_tools.get(tool_name, "")
                f.write("        \"{}\": \"{}\",\n".format(tool_name, tool_path))
            f.write("\n")
        
        f.write("    },\n")
        f.write("    \"tool_configs\": {}\n")
        f.write("}\n\n")
        
        # Add auto-detection code
        f.write("# Try to auto-detect tools\n")
        f.write("try:\n")
        f.write("    detected_tools = detect_tool_locations()\n")
        f.write("    config[\"tool_configs\"] = detected_tools\n")
        f.write("except Exception as e:\n")
        f.write("    print(f\"Warning: Could not auto-detect tools: {e}\")\n")
        f.write("    # Continue with default paths\n\n")
        
        # Add user config support
        f.write("# Allow for user config to override default settings\n")
        f.write("user_config_path = os.environ.get('METAMAG_CONFIG', os.path.expanduser('~/.metamag/config.yaml'))\n")
        f.write("if os.path.exists(user_config_path):\n")
        f.write("    try:\n")
        f.write("        with open(user_config_path, 'r') as f:\n")
        f.write("            user_config = yaml.safe_load(f)\n")
        f.write("        \n")
        f.write("        # Update config with user settings\n")
        f.write("        if user_config.get('environment'):\n")
        f.write("            config['environment'].update(user_config['environment'])\n")
        f.write("        if user_config.get('tools'):\n")
        f.write("            config['tools'].update(user_config['tools'])\n")
        f.write("        if user_config.get('tool_configs'):\n")
        f.write("            config['tool_configs'].update(user_config['tool_configs'])\n")
        f.write("            \n")
        f.write("    except Exception as e:\n")
        f.write("        print(f\"Warning: Could not load user config from {user_config_path}: {e}\")\n")
    
    with open(config_path, "a") as f:  # Use append mode
        f.write("\n# Helper functions for backward compatibility and conda support\n")
        f.write("def get_tool_path(tool_name):\n")
        f.write("    \"\"\"\n")
        f.write("    Get the correct tool path from either direct tools config or detected tools.\n")
        f.write("    This function provides backward compatibility with old config structure.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        tool_name (str): Name of the tool to get path for\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        str: Path to the tool executable, or None if not found\n")
        f.write("    \"\"\"\n")
        f.write("    # First check if tool is in the direct tools config\n")
        f.write("    if tool_name in config[\"tools\"]:\n")
        f.write("        return config[\"tools\"][tool_name]\n")
        f.write("    \n")
        f.write("    # Then check if tool is in detected tool configs\n")
        f.write("    if tool_name in config.get(\"tool_configs\", {}):\n")
        f.write("        tool_config = config[\"tool_configs\"][tool_name]\n")
        f.write("        if tool_config and isinstance(tool_config, dict):\n")
        f.write("            return tool_config.get(\"path\")\n")
        f.write("        return tool_config\n")
        f.write("    \n")
        f.write("    # Tool not found\n")
        f.write("    print(f\"Warning: Tool '{tool_name}' not found in configuration\")\n")
        f.write("    return None\n\n")
        
        f.write("def get_tool_command(tool_name):\n")
        f.write("    \"\"\"\n")
        f.write("    Get the complete command to run a tool, including conda activation if needed.\n")
        f.write("    \n")
        f.write("    Args:\n")
        f.write("        tool_name (str): Name of the tool\n")
        f.write("        \n")
        f.write("    Returns:\n")
        f.write("        tuple: (command_prefix, tool_path) where command_prefix includes conda activation if needed\n")
        f.write("    \"\"\"\n")
        f.write("    # First check direct tools config\n")
        f.write("    if tool_name in config[\"tools\"]:\n")
        f.write("        return \"\", config[\"tools\"][tool_name]\n")
        f.write("    \n")
        f.write("    # Check detected tool configs\n")
        f.write("    if tool_name in config.get(\"tool_configs\", {}):\n")
        f.write("        tool_config = config[\"tool_configs\"][tool_name]\n")
        f.write("        if tool_config and isinstance(tool_config, dict):\n")
        f.write("            if tool_config.get(\"requires_activation\", False):\n")
        f.write("                activate_cmd = tool_config.get(\"activate_cmd\", \"\")\n")
        f.write("                return f\"source {activate_cmd} && \", tool_config.get(\"path\")\n")
        f.write("            else:\n")
        f.write("                return \"\", tool_config.get(\"path\")\n")
        f.write("        return \"\", tool_config\n")
        f.write("    \n")
        f.write("    print(f\"Warning: Tool '{tool_name}' not found in configuration\")\n")
        f.write("    return \"\", None\n")
    
    # Print summary (keep your existing print statements)
    found_tools = len(available_tools)
    total_tools = sum(len(tools) for tools in TOOL_GROUPS.values())
    missing_tools = total_tools - found_tools
    
    print("\nFound {} out of {} tools in PATH.".format(found_tools, total_tools))
    print("{} tools were not found and have empty paths in the config.".format(missing_tools))
    
    if missing_tools > 0:
        print("\nYou have several options for missing tools:")
        print("  1. Install missing tools with Conda: python setup_tools.py --tools [tool1] [tool2] ...")
        print("  2. Install missing tools manually and update config.py with their paths")
        print("  3. Skip missing tools if they're not needed for your workflow")

def create_custom_environment_file(tool_groups):
    """Create a custom environment file with tools from specified groups."""
    with open("custom_environment.yml", "w") as f:
        f.write("name: metamag\n")
        f.write("channels:\n")
        f.write("  - conda-forge\n")
        f.write("  - bioconda\n")
        f.write("  - defaults\n")
        f.write("dependencies:\n")
        f.write("  - python=3.9\n")
        
        # Add tools from specified groups
        for group in tool_groups:
            if group in TOOL_GROUPS:
                f.write("  # {} tools\n".format(group))
                for tool in TOOL_GROUPS[group]:
                    # Skip tools that might not be in conda
                    if tool in ["idba_fq2fa"]:
                        continue
                    f.write("  - {}\n".format(tool))
        
        # Always include essential Python packages
        f.write("  # Essential packages\n")
        f.write("  - pandas\n")
        f.write("  - numpy\n")
        f.write("  - pyyaml\n")
        f.write("  - biopython\n")
if __name__ == "__main__":
    main()

import os
import shutil
import yaml
import subprocess
import json

# Helper functions for backward compatibility and conda support
def get_tool_path(tool_name):
    """
    Get the correct tool path from either direct tools config or detected tools.
    This function provides backward compatibility with old config structure.
    
    Args:
        tool_name (str): Name of the tool to get path for
        
    Returns:
        str: Path to the tool executable, or None if not found
    """
    # First check if tool is in the direct tools config
    if tool_name in config["tools"]:
        return config["tools"][tool_name]
    
    # Then check if tool is in detected tool configs
    if tool_name in config.get("tool_configs", {}):
        tool_config = config["tool_configs"][tool_name]
        if tool_config and isinstance(tool_config, dict):
            return tool_config.get("path")
        return tool_config
    
    # Tool not found
    print(f"Warning: Tool '{tool_name}' not found in configuration")
    return None

def get_tool_command(tool_name):
    """
    Get the complete command to run a tool, including conda activation if needed.
    
    Args:
        tool_name (str): Name of the tool
        
    Returns:
        tuple: (command_prefix, tool_path) where command_prefix includes conda activation if needed
    """
    # First check direct tools config
    if tool_name in config["tools"]:
        return "", config["tools"][tool_name]
    
    # Check detected tool configs
    if tool_name in config.get("tool_configs", {}):
        tool_config = config["tool_configs"][tool_name]
        if tool_config and isinstance(tool_config, dict):
            if tool_config.get("requires_activation", False):
                activate_cmd = tool_config.get("activate_cmd", "")
                return f"source {activate_cmd} && ", tool_config.get("path")
            else:
                return "", tool_config.get("path")
        return "", tool_config
    
    print(f"Warning: Tool '{tool_name}' not found in configuration")
    return "", None

def get_conda_base_path():
    """
    Dynamically detect conda installation path.
    Returns the conda base path or None if not found.
    """
    # Method 1: Try conda info command
    try:
        result = subprocess.run(['conda', 'info', '--json'], 
                               capture_output=True, text=True, check=True)
        conda_info = json.loads(result.stdout)
        conda_root = conda_info.get('root_prefix')
        if conda_root and os.path.exists(conda_root):
            return conda_root
    except Exception:
        pass
    
    # Method 2: Try to find from CONDA_PREFIX environment variable
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        # If we're in a conda environment, get the base
        base_path = conda_prefix
        while base_path and base_path != '/' and 'envs' in base_path:
            base_path = os.path.dirname(base_path)
        if base_path and (base_path.endswith('miniconda3') or base_path.endswith('anaconda3')):
            return base_path
    
    # Method 3: Check common conda installation locations
    import getpass
    username = getpass.getuser()
    common_paths = [
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/anaconda3"),
        f"/usr/home/{username}/.local/bin/miniconda3",
        f"/usr/home/{username}/.local/bin/anaconda3",
        "/opt/miniconda3",
        "/opt/anaconda3",
        "/usr/local/miniconda3",
        "/usr/local/anaconda3"
    ]
    
    for path in common_paths:
        if os.path.exists(os.path.join(path, "etc", "profile.d", "conda.sh")):
            return path
    
    print("[WARNING] Could not automatically detect conda installation path")
    return None

def detect_r_environment():
    """
    Detect R environment configuration.
    Returns R environment configuration or None if not found.
    """
    conda_base = get_conda_base_path()
    if not conda_base:
        return None
    
    # Check for common R environment names
    r_env_names = ['R', 'r-base', 'rstats', 'r-env']
    
    for env_name in r_env_names:
        env_path = os.path.join(conda_base, "envs", env_name)
        if os.path.exists(env_path):
            r_bin = os.path.join(env_path, "bin", "R")
            rscript_bin = os.path.join(env_path, "bin", "Rscript")
            
            if os.path.exists(r_bin) and os.path.exists(rscript_bin):
                return {
                    "env_name": env_name,
                    "env_path": env_path,
                    "conda_base": conda_base,
                    "r_bin": r_bin,
                    "rscript_bin": rscript_bin,
                    "r_home": os.path.join(env_path, "lib", "R"),
                    "r_libs": os.path.join(env_path, "lib", "R", "library"),
                    "activate_cmd": f"conda activate {env_name}"
                }
    
    # Check if R is available in system PATH
    r_system = shutil.which('R')
    rscript_system = shutil.which('Rscript')
    
    if r_system and rscript_system:
        try:
            # Try to get R home
            result = subprocess.run(['R', '--slave', '-e', 'cat(R.home())'], 
                                   capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                r_home = result.stdout.strip()
                return {
                    "env_name": "system",
                    "env_path": None,
                    "conda_base": conda_base,
                    "r_bin": r_system,
                    "rscript_bin": rscript_system,
                    "r_home": r_home,
                    "r_libs": os.path.join(r_home, "library"),
                    "activate_cmd": None
                }
        except Exception:
            pass
    
    return None

def get_r_config():
    """
    Get R environment configuration.
    
    Returns:
        dict: R configuration dictionary or None if not found
    """
    return config.get("r_config")

def setup_r_environment():
    """
    Set up R environment variables based on detected configuration.
    
    Returns:
        bool: True if R environment was set up successfully, False otherwise
    """
    r_config = get_r_config()
    if not r_config:
        print("[WARNING] No R configuration found")
        return False
    
    # Set environment variables
    if r_config.get("r_home"):
        os.environ["R_HOME"] = r_config["r_home"]
    
    if r_config.get("r_libs"):
        os.environ["R_LIBS_USER"] = r_config["r_libs"]
    
    if r_config.get("env_path"):
        # Add R bin to PATH
        r_bin_dir = os.path.join(r_config["env_path"], "bin")
        os.environ["PATH"] = f"{r_bin_dir}:" + os.environ.get("PATH", "")
        
        # Add R lib to LD_LIBRARY_PATH
        r_lib_dir = os.path.join(r_config["env_path"], "lib")
        if os.path.exists(r_lib_dir):
            os.environ["LD_LIBRARY_PATH"] = f"{r_lib_dir}:" + os.environ.get("LD_LIBRARY_PATH", "")
    
    return True




def find_tool(tool_name, default_path=None):
    """Try to find a tool in PATH, otherwise return the default path."""
    path_result = shutil.which(tool_name)
    if path_result:
        return path_result
    return default_path

def is_conda_available():
    """Check if conda is available in the system."""
    return shutil.which('conda') is not None

def detect_tool_locations():
    """Detect tool locations based on system paths and conda environments."""
    tools = {}
    
    # Dictionary of tool names and their possible conda environments
    conda_environments = {
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
        "checkm": "checkm_env",
        "drep": "metamag",
        "gtdbtk": "gtdbtk",
        "kraken2": "metamag",
        "bracken": "metamag",
        "eggnog_mapper": "metamag",
        "prodigal": "metamag",
        "dbcan": "dbcan",
    }
    
    # Map of expected tool names to their executable names
    tool_to_exec = {
        "eggnog_mapper": "emapper.py",
        "eggnog_mapper": "emapper.py",
        "dbcan": "run_dbcan.py",
        "das_tool": "DAS_Tool",
        "maxbin2": "run_MaxBin.pl",
        "drep": "dRep",
        "metaquast": "metaquast.py",
        "jgi_summarize": "jgi_summarize_bam_contig_depths",
    }
    
    # Check if conda is available
    conda_available = is_conda_available()
    
    # Try to find conda root path
    conda_root = None
    if conda_available:
        try:
            result = subprocess.run(['conda', 'info', '--json'], 
                                   capture_output=True, text=True, check=True)
            import json
            conda_info = json.loads(result.stdout)
            conda_root = conda_info.get('root_prefix')
        except Exception as e:
            print(f"Warning: Conda available but couldn't get info: {e}")
    
    # Process each tool
    for tool_name, env_name in conda_environments.items():
        # Get the executable name if it differs
        exec_name = tool_to_exec.get(tool_name, tool_name)
        
        # First try to find the tool in PATH
        direct_path = shutil.which(exec_name)
        if direct_path:
            tools[tool_name] = {
                "path": direct_path,
                "requires_activation": False
            }
            continue
        
        # If not in PATH and conda is available, check in conda environment
        if conda_available and conda_root:
            # Check if environment exists
            env_path = os.path.join(conda_root, "envs", env_name)
            if os.path.exists(env_path):
                # Check if tool exists in this environment
                tool_in_env = os.path.join(env_path, "bin", exec_name)
                if os.path.exists(tool_in_env):
                    tools[tool_name] = {
                        "path": tool_in_env,
                        "requires_activation": True,
                        "env_name": env_name,
                        "activate_cmd": f"conda activate {env_name}"
                    }
                    continue
        
        # Tool not found in PATH or conda
        tools[tool_name] = None
    
    return tools

# Base config - includes all tools needed for the pipeline
config = {
    "environment": {
         "PERL5LIB": "/opt/ghpc/lib64/perl-5.36",
         "conda_env": "source /usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate metamag",
    },
    "tools": {
        # qc tools
        "fastqc": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/metaforge/bin/fastqc",
        "multiqc": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/metaforge/bin/multiqc",

        # preprocessing tools
        "fastp": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/2.1_Trim_short/fastp-0.23.2/fastp",
        "bwa": "/opt/ghpc/bin/bwa",
        "samtools": "/opt/ghpc/bin/samtools",
        
        # assembly tools
        "idba": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/3.1_assemble_short/idba-1.1.3/bin/idba_ud",
        "idba_fq2fa": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/3.1_assemble_short/idba-1.1.3/bin/fq2fa",
        "megahit": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/3.1_assemble_short/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit",
        "metaquast": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/quast/bin/metaquast.py",


        # binning tools
        "jgi_summarize": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/metabat_env/bin/jgi_summarize_bam_contig_depths",
        "das_tool": "/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate das_env",
        "DAS_Tool_Fasta2Contig": "/usr/home/workspace/tmp.0/zexi/app/DAS_Tool-1.1.4/src/Fasta_to_Contig2Bin.sh",
        "metawrap": "/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate metawrap_env",
        

        # evaluation tools
        "checkm": "/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate checkm_env",
        "drep": "/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate drep",

        # taxonomy tools
        "gtdbtk": "/usr/home/qgg/maralta/.local/bin/miniconda3/bin/activate gtdbtk-2.4.1",
        "kraken2": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/6_Taxonomy/kraken2-2.1.2/kraken2",
        "kraken2_build": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/6_Taxonomy/kraken2-2.1.2/kraken2-build",
        "bracken": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/9_Abundance/Bracken-2.7/src/est_abundance.py",
        "bracken_build": "/usr/home/workspace/tmp.0/zexi/app/Bracken-2.7/bracken-build",


        # annotation tools
        "eggnog_mapper": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/4.2_gene_annoation/eggnog-mapper-2.1.11/emapper.py",
        "eggnog_db_dir": "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/db/eggNOGdb",
        "dbcan": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/dbcan/bin/run_dbcan.py",
        "dbcan_db_dir": "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/dbcan/databases/",
        "dbcan_activate": "/usr/home/qgg/maralta/.local/bin/miniconda3/etc/profile.d/conda.sh && conda activate dbcan",
        "prodigal": "/usr/home/workspace/tmp.0/zexi/app/Prodigal_v2.6.3/prodigal.linux",
    },
    "tool_configs": {},
    "r_config": None
}

# Try to auto-detect tools
try:
    detected_tools = detect_tool_locations()
    config["tool_configs"] = detected_tools
except Exception as e:
    print(f"Warning: Could not auto-detect tools: {e}")
    # Continue with default paths

# Try to auto-detect R environment
try:
    r_config = detect_r_environment()
    config["r_config"] = r_config
    if r_config:
        print(f"[INFO] Detected R environment: {r_config['env_name']} at {r_config['env_path'] or 'system'}")
    else:
        print("[WARNING] No R environment detected")
except Exception as e:
    print(f"Warning: Could not auto-detect R environment: {e}")

# Allow for user config to override default settings
user_config_path = os.environ.get('METAMAG_CONFIG', os.path.expanduser('~/.metamag/config.yaml'))
if os.path.exists(user_config_path):
    try:
        with open(user_config_path, 'r') as f:
            user_config = yaml.safe_load(f)
        
        # Update config with user settings
        if user_config.get('environment'):
            config['environment'].update(user_config['environment'])
        if user_config.get('tools'):
            config['tools'].update(user_config['tools'])
        if user_config.get('tool_configs'):
            config['tool_configs'].update(user_config['tool_configs'])
            
    except Exception as e:
        print(f"Warning: Could not load user config from {user_config_path}: {e}")


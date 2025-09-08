import os
import sys
import shutil
import subprocess
import pandas as pd
import tempfile
from Bio import SeqIO
from collections import defaultdict

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist."""
    os.makedirs(directory, exist_ok=True)

def run_command(command):
    """Run a shell command with proper error handling."""
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True
        )
        print(f"[INFO] Command executed successfully: {command}")
        if result.stdout:
            print(f"[INFO] Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {command}")
        print(f"[ERROR] Return code: {e.returncode}")
        if e.stdout:
            print(f"[ERROR] stdout: {e.stdout}")
        if e.stderr:
            print(f"[ERROR] stderr: {e.stderr}")
        return False

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

def get_tool_command(tool_name):
    """Get tool command - wrapper around get_tool_path for compatibility"""
    return get_tool_path(tool_name)

def find_script_in_main_directory(script_name):
    """Find a script in the main script directory of the pipeline"""
    # Get the directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Look for the script in the main script directory
    # Assuming the main script directory is where main.py and other scripts are located
    main_script_dir = os.path.dirname(current_dir)  # Go up one level from MetaMAG package
    
    # Check multiple possible locations
    possible_paths = [
        os.path.join(main_script_dir, script_name),
        os.path.join(main_script_dir, "scripts", script_name),
        os.path.join(current_dir, script_name),  # Same directory as this file
        os.path.join(current_dir, "..", script_name),  # Parent directory
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            print(f"[INFO] Found script: {path}")
            return path
    
    print(f"[ERROR] Could not find script {script_name} in any of these locations:")
    for path in possible_paths:
        print(f"  - {path}")
    return None

def add_nmags_to_kraken(project_output, cpus, kraken_db_path=None, rumen_ref_mags_dir=None, 
                       read_length=150, merge_mode=True, custom_gtdbtk_dir=None,
                       advanced_script_path=None):
    """
    Add Novel MAGs to GTDB-Kraken2 database using advanced taxonomy resolver and build Bracken database.
    
    Steps:
    1. Check existing database and get next taxonomy ID
    2. Placeholder taxa assignment
    3. GTDB to Taxdump Conversion using advanced resolver (with merge support)
    4. Kraken-Compatible Header Assignment
    5. Add Kraken Library
    6. Build Kraken Database
    7. Build Bracken Database
    8. Build Database Distribution
    
    Args:
        project_output (str): Project output directory
        cpus (int): Number of CPUs to use
        kraken_db_path (str): Optional user-specified path to build Kraken database
                             If not provided, defaults to project_output/Kraken_Database
        rumen_ref_mags_dir (str): Path to rumen reference MAGs directory (optional)
        read_length (int): Read length for Bracken database building (default: 150)
        merge_mode (bool): Whether to merge with existing database (default: True)
        custom_gtdbtk_dir (str): Custom path to GTDB-Tk results directory (optional)
        advanced_script_path (str): Path to gtdb_to_taxdump_latest_resolve.py script (optional)
    
    Returns:
        bool: Success or failure
    """
    print("=" * 70)
    print("Starting Kraken Database Building Process with Advanced Resolver")
    print("=" * 70)
    print(f"Project output: {project_output}")
    print(f"CPUs: {cpus}")
    print(f"Kraken DB path: {kraken_db_path}")
    print(f"Rumen ref MAGs dir: {rumen_ref_mags_dir}")
    print(f"Read length: {read_length}")
    print(f"Merge mode: {merge_mode}")
    print(f"Custom GTDB-Tk dir: {custom_gtdbtk_dir}")
    print("=" * 70)
    
    # Define directories
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    true_novel_mags_dir = os.path.join(project_output, "Novel_Mags", "true_novel_MAGs")
    
    # Set Kraken database path - use user-provided path if available
    if not kraken_db_path:
        kraken_db_path = os.path.join(project_output, "Kraken_Database")
    
    print(f"[INFO] Using Kraken database path: {kraken_db_path}")
    
    # Create subdirectories in the Kraken database directory
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    processed_data_dir = os.path.join(kraken_db_path, "processed_data")
    taxonomy_dir = os.path.join(kraken_db_path, "taxonomy")
    library_dir = os.path.join(kraken_db_path, "library")
    
    # Create backup directory for existing taxonomy files
    backup_dir = os.path.join(kraken_db_path, "taxonomy_backup")
    
    # Create required directories
    for directory in [kraken_db_path, scripts_dir, processed_data_dir, taxonomy_dir, library_dir, backup_dir]:
        ensure_directory_exists(directory)
    
    # Check if existing database exists and get taxonomy file paths
    existing_nodes = os.path.join(taxonomy_dir, "nodes.dmp")
    existing_names = os.path.join(taxonomy_dir, "names.dmp")
    
    # Find the resolve script in the main script directory
    if not advanced_script_path:
        advanced_script_path = find_script_in_main_directory("gtdb_to_taxdump_latest_resolve.py")
        if not advanced_script_path:
            print("[ERROR] Could not find gtdb_to_taxdump_latest_resolve.py script")
            return False
    
    print(f"[INFO] Using resolve script: {advanced_script_path}")
    
    # Check if MAGs repository exists
    if not os.path.exists(mags_repo_dir):
        print(f"[ERROR] MAGs Repository directory not found: {mags_repo_dir}")
        return False
    
    # Count MAGs in repository and novel MAGs
    repo_mags = [f for f in os.listdir(mags_repo_dir) if f.endswith(('.fa', '.fasta'))]
    
    if len(repo_mags) == 0:
        print(f"[WARNING] No MAGs found in repository: {mags_repo_dir}")
        return False
    
    print(f"[INFO] Found {len(repo_mags)} MAGs in repository")
    
    # Check for reference MAGs in rumen-specific directory
    ref_mags = []
    ref_taxonomy_dir = None
    mags_subdir = None
    
    if rumen_ref_mags_dir and os.path.exists(rumen_ref_mags_dir):
        # Check if it has mags and taxonomy subdirectories
        mags_subdir = os.path.join(rumen_ref_mags_dir, "mags")
        tax_subdir = os.path.join(rumen_ref_mags_dir, "taxonomy")
        
        if os.path.exists(mags_subdir):
            ref_mags = [f for f in os.listdir(mags_subdir) if f.endswith(('.fa', '.fasta'))]
            print(f"[INFO] Found {len(ref_mags)} reference MAGs in {mags_subdir}")
        else:
            # Check if the main directory has FASTA files
            ref_mags = [f for f in os.listdir(rumen_ref_mags_dir) if f.endswith(('.fa', '.fasta'))]
            print(f"[INFO] Found {len(ref_mags)} reference MAGs in {rumen_ref_mags_dir}")
            mags_subdir = rumen_ref_mags_dir
        
        if os.path.exists(tax_subdir):
            ref_taxonomy_dir = tax_subdir
            print(f"[INFO] Found taxonomy information in {tax_subdir}")
    
    print(f"[INFO] Adding {len(repo_mags)} MAGs to Kraken2 database at {kraken_db_path}")
    
    # Step 1: Placeholder taxa assignment with your custom script
    print(f"[INFO] Step 1: Running custom placeholder taxa assignment")
    
    # Get the list of all MAGs to process (both repo MAGs and novel MAGs)
    all_mags = []
    for mag_file in repo_mags:
        mag_path = os.path.join(mags_repo_dir, mag_file)
        all_mags.append((os.path.splitext(mag_file)[0], mag_path))
    
    # Also include reference MAGs if they exist
    if ref_mags and mags_subdir:
        for mag_file in ref_mags:
            mag_path = os.path.join(mags_subdir, mag_file)
            all_mags.append((os.path.splitext(mag_file)[0], mag_path))
    
    # Save MAG list to a temporary file
    tax_input_file = os.path.join(processed_data_dir, "mags_list.xlsx")
    df = pd.DataFrame({
        'user_genome': [mag[0] for mag in all_mags],
        'classification': ["d__Bacteria;p__;c__;o__;f__;g__;s__"] * len(all_mags)  # Default classification
    })
    df.to_excel(tax_input_file, index=False)
    
    # Get GTDB-Tk classifications if available - USE CUSTOM DIRECTORY IF PROVIDED
    if custom_gtdbtk_dir:
        gtdbtk_dir = custom_gtdbtk_dir
        print(f"[INFO] Using custom GTDB-Tk directory: {gtdbtk_dir}")
    else:
        gtdbtk_dir = os.path.join(project_output, "Novel_Mags", "gtdbtk")
        print(f"[INFO] Using default GTDB-Tk directory: {gtdbtk_dir}")
    
    bac_summary = os.path.join(gtdbtk_dir, "classify", "gtdbtk.bac120.summary.tsv")
    ar_summary = os.path.join(gtdbtk_dir, "classify", "gtdbtk.ar53.summary.tsv")
    
    # Create the placeholder taxa assignment script
    placeholder_script = os.path.join(scripts_dir, "1_placeholder_taxa.py")
    tax_fixed_file = os.path.join(processed_data_dir, "tax_fixed.xlsx")
    
    with open(placeholder_script, 'w') as f:
        f.write(f"""
import pandas as pd
import os

# Load the input data
input_file = '{tax_input_file}'
data = pd.read_excel(input_file)
print(f"Loaded {{len(data)}} rows from {{input_file}}")

# Function to clean data and replace all empty levels without underscore
def clean_and_replace_all_levels(row):
    \"\"\"
    1. Remove duplicated taxonomy info
    2. Make sure it starts with domain (d__)
    3. Replace ALL empty taxonomic levels with appropriate prefix + genome name + asterisk
       without underscore between level prefix and genome name
    \"\"\"
    classification = row['classification']
    genome = row['user_genome']
    
    # Check if taxonomy appears to be duplicated (contains d__Bacteria after the start)
    if ';d__Bacteria;' in classification:
        # Split and take the part that starts with domain
        parts = classification.split(';d__Bacteria;')
        if len(parts) >= 2:
            classification = 'd__Bacteria;' + parts[1]
    
    # If it still doesn't start with domain, add it
    if not classification.startswith('d__'):
        if classification.startswith('p__'):
            classification = 'd__Bacteria;' + classification
    
    # Split the classification string into parts
    levels = classification.split(';')
    
    # Define the prefixes for each taxonomic level
    prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    level_names = ['d', 'p', 'c', 'o', 'f', 'g', 'sp']  # 'sp' for species level
    
    # Process each level and replace if empty
    for i in range(min(len(levels), len(prefixes))):
        level = levels[i].strip()
        prefix = prefixes[i]
        level_name = level_names[i]
        
        # Check if this level is empty (just the prefix)
        if level == prefix:
            # For species level, use "sp" prefix without underscore
            if prefix == 's__':
                levels[i] = f"{{prefix}}{{level_name}}{{genome}}*"
            else:
                # For other levels, use the level's letter without underscore
                levels[i] = f"{{prefix}}{{level_name}}{{genome}}*"
    
    # Rejoin the levels into a classification string
    return ';'.join(levels)

# Update with GTDB-Tk classifications if available
bac_path = '{bac_summary}'
ar_path = '{ar_summary}'

if os.path.exists(bac_path):
    bac_df = pd.read_csv(bac_path, sep='\\t')
    for _, row in bac_df.iterrows():
        idx = data[data['user_genome'] == row['user_genome']].index
        if len(idx) > 0:
            data.loc[idx, 'classification'] = row['classification']

if os.path.exists(ar_path):
    ar_df = pd.read_csv(ar_path, sep='\\t')
    for _, row in ar_df.iterrows():
        idx = data[data['user_genome'] == row['user_genome']].index
        if len(idx) > 0:
            data.loc[idx, 'classification'] = row['classification']

# If reference taxonomy directory exists, use that information
ref_taxonomy_dir = '{ref_taxonomy_dir}'
if ref_taxonomy_dir and os.path.exists(ref_taxonomy_dir):
    for file in os.listdir(ref_taxonomy_dir):
        if file.endswith('.tsv') and not file.startswith('.'):
            try:
                # Try to load the taxonomy file (assuming it's tab-separated)
                tax_file = os.path.join(ref_taxonomy_dir, file)
                print(f"Loading reference taxonomy from: {{tax_file}}")
                ref_tax_df = pd.read_csv(tax_file, sep='\\t')
                
                # Check if it has the expected columns
                if 'user_genome' in ref_tax_df.columns and 'classification' in ref_tax_df.columns:
                    for _, row in ref_tax_df.iterrows():
                        # Find genome in our data
                        idx = data[data['user_genome'] == row['user_genome']].index
                        if len(idx) > 0:
                            data.loc[idx, 'classification'] = row['classification']
                            print(f"Updated taxonomy for {{row['user_genome']}}")
            except Exception as e:
                print(f"Error loading reference taxonomy file {{file}}: {{e}}")

# Create a new clean dataframe
clean_data = pd.DataFrame()
clean_data['user_genome'] = data['user_genome']
clean_data['classification'] = data.apply(clean_and_replace_all_levels, axis=1)

# Save the results
output_file = '{tax_fixed_file}'
clean_data.to_excel(output_file, index=False)
print(f"Saved cleaned data to {{output_file}}")
""")
    
    # Run the placeholder taxa assignment script
    if not run_command(f"python {placeholder_script}"):
        print("[ERROR] Failed to run placeholder taxa assignment")
        return False
    
    # Step 2: Convert to TSV for GTDB to taxdump conversion
    print(f"[INFO] Step 2: Converting taxonomy to TSV format")
    
    # Create the TSV conversion script
    tsv_conversion_script = os.path.join(scripts_dir, "2_convert_to_tsv.py")
    tax_tsv_file = os.path.join(processed_data_dir, "taxonomy.tsv")
    tax_tsv_noheader_file = os.path.join(processed_data_dir, "taxonomy_noheader.tsv")
    
    with open(tsv_conversion_script, 'w') as f:
        f.write(f"""
import pandas as pd

# Input Excel file
input_file = '{tax_fixed_file}'

# Output TSV file
output_file = '{tax_tsv_file}'
noheader_file = '{tax_tsv_noheader_file}'

# Read the Excel file
print(f"Reading Excel file from: {{input_file}}")
data = pd.read_excel(input_file)

# Save as TSV file (tab-separated values)
print(f"\\nSaving data to TSV file: {{output_file}}")
data.to_csv(output_file, sep='\\t', index=False)

# Save as TSV file without header
print(f"\\nSaving data to TSV file without header: {{noheader_file}}")
data.to_csv(noheader_file, sep='\\t', index=False, header=False)

print(f"Conversion complete! TSV files saved.")
""")
    
    # Run the TSV conversion script
    if not run_command(f"python {tsv_conversion_script}"):
        print("[ERROR] Failed to run TSV conversion")
        return False
    
    # Step 3: GTDB to taxdump conversion using advanced resolver
    print(f"[INFO] Step 3: Converting GTDB taxonomy to NCBI taxdump format using advanced resolver")
    
    # Copy the advanced resolver script to scripts directory
    advanced_resolver_script = os.path.join(scripts_dir, "3_gtdb_to_taxdump_latest_resolve.py")
    
    try:
        shutil.copy2(advanced_script_path, advanced_resolver_script)
        os.chmod(advanced_resolver_script, 0o755)
        print(f"[INFO] Copied advanced resolver script from {advanced_script_path}")
    except Exception as e:
        print(f"[ERROR] Failed to copy advanced resolver script: {e}")
        return False
    
    # Prepare command arguments based on merge mode
    if merge_mode and os.path.exists(existing_nodes):
        # Backup existing files before merging
        try:
            shutil.copy2(existing_nodes, os.path.join(backup_dir, "nodes.dmp.backup"))
            shutil.copy2(existing_names, os.path.join(backup_dir, "names.dmp.backup"))
            print(f"[INFO] Backed up existing taxonomy files to {backup_dir}")
        except Exception as e:
            print(f"[WARNING] Could not backup existing files: {e}")
        
        # Use existing taxonomy for merging
        resolver_command = f"""python3 {advanced_resolver_script} \\
    --existing-nodes {existing_nodes} \\
    --existing-names {existing_names} \\
    --new-taxonomy {tax_tsv_noheader_file} \\
    --output-dir {taxonomy_dir}"""
        print(f"[INFO] Running in merge mode with existing taxonomy")
    else:
        # For non-merge mode, we need to create minimal existing taxonomy first
        temp_existing_dir = os.path.join(processed_data_dir, "temp_existing")
        ensure_directory_exists(temp_existing_dir)
        
        # Create minimal root taxonomy
        temp_nodes = os.path.join(temp_existing_dir, "nodes.dmp")
        temp_names = os.path.join(temp_existing_dir, "names.dmp")
        
        with open(temp_nodes, 'w') as f:
            f.write("1\t|\t1\t|\tno rank\t|\t\t|\n")
        
        with open(temp_names, 'w') as f:
            f.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        
        resolver_command = f"""python3 {advanced_resolver_script} \\
    --existing-nodes {temp_nodes} \\
    --existing-names {temp_names} \\
    --new-taxonomy {tax_tsv_noheader_file} \\
    --output-dir {taxonomy_dir}"""
        print(f"[INFO] Running in new database mode")
    
    # Run the advanced GTDB to taxdump conversion
    if not run_command(resolver_command):
        print("[ERROR] Failed to run advanced GTDB to taxdump conversion")
        return False
    
    # The resolver script should have created new nodes.dmp, names.dmp, and taxid.map files
    # in the taxonomy_dir, replacing the old ones
    print(f"[INFO] Taxonomy files updated in {taxonomy_dir}")
    
    # Step 4: Rename headers for Kraken compatibility using the header script from main directory
    print(f"[INFO] Step 4: Renaming headers for Kraken compatibility")
    
    # Find the header renaming script in the main script directory
    header_script_path = find_script_in_main_directory("header_kraken.py")
    if not header_script_path:
        print("[ERROR] Could not find header_kraken.py script")
        return False
    
    print(f"[INFO] Using header script: {header_script_path}")
    
    # Create the corrected MAGs directory
    corrected_mags_dir = os.path.join(processed_data_dir, "mags_corrected_header")
    ensure_directory_exists(corrected_mags_dir)
    
    # Copy the header script to the local scripts directory and modify it for our paths
    local_header_script = os.path.join(scripts_dir, "4_rename_headers.py")
    
    # Read the original header script and modify the paths
    with open(header_script_path, 'r') as f:
        header_script_content = f.read()
    
    # Replace the paths in the header script
    taxid_map_file = os.path.join(taxonomy_dir, "taxid.map")
    
    # Create modified header script with correct paths
    modified_header_content = f"""import os
import multiprocessing as mp
from Bio import SeqIO
from functools import partial

# Path to your MAG files folder
mags_folder_path = '{mags_repo_dir}'
# Path to your taxid.map file  
taxid_map_path = '{taxid_map_file}'
# Output folder where renamed files will be saved
output_folder = '{corrected_mags_dir}'

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Load the taxid.map file into a dictionary
def load_taxid_map(map_path):
    genome_to_taxid = {{}}
    
    with open(map_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) == 2:
                genome_id, taxid = parts
                genome_to_taxid[genome_id] = taxid
    
    return genome_to_taxid

# Function to rename headers based on the filename and taxid map
def rename_fasta_headers(fasta_file, mags_folder_path, output_folder, genome_to_taxid):
    input_file_path = os.path.join(mags_folder_path, fasta_file)
    output_file_path = os.path.join(output_folder, fasta_file)
    
    # Extract the genome name from the filename (remove extension)
    genome_name = fasta_file.replace('.fa', '')
    
    # Look up the taxid for this genome
    taxid = genome_to_taxid.get(genome_name, 'UnknownTaxID')
    
    try:
        with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
            # Initialize counter for contig numbers
            contig_counter = 1
            
            for record in SeqIO.parse(input_handle, 'fasta'):
                # Create the new standardized header with sequential contig numbers
                new_id = f"{{genome_name}}_{{contig_counter}}|kraken:taxid|{{taxid}}"
                contig_counter += 1
                
                # Update the record
                record.id = new_id
                record.description = ''
                
                # Write the renamed record
                SeqIO.write(record, output_handle, 'fasta')
        
        return f"Successfully processed {{fasta_file}} - renamed {{contig_counter-1}} contigs"
    except Exception as e:
        return f"Error processing {{fasta_file}}: {{str(e)}}"

def main():
    # Load taxid map
    genome_to_taxid = load_taxid_map(taxid_map_path)
    
    # Get list of FASTA files
    fasta_files = [f for f in os.listdir(mags_folder_path) if f.endswith('.fa')]
    
    # Set up the process pool
    num_cpus = {cpus}  # Use specified CPUs
    print(f"Using {{num_cpus}} CPUs for processing")
    
    # Create a partial function with the fixed arguments
    process_file = partial(
        rename_fasta_headers,
        mags_folder_path=mags_folder_path,
        output_folder=output_folder,
        genome_to_taxid=genome_to_taxid
    )
    
    # Process files in parallel
    with mp.Pool(processes=num_cpus) as pool:
        results = pool.map(process_file, fasta_files)
    
    # Print results
    for result in results:
        print(result)
    
    print("Renaming completed. The renamed files are saved in:", output_folder)

if __name__ == "__main__":
    main()
"""
    
    # Write the modified header script
    with open(local_header_script, 'w') as f:
        f.write(modified_header_content)
    
    # Run the header renaming script
    if not run_command(f"python {local_header_script}"):
        print("[ERROR] Failed to run header renaming")
        return False
    
    # Step 5: Combine all MAGs into a single file
    print(f"[INFO] Step 5: Combining all MAGs into a single file")
    
    combined_mags_file = os.path.join(processed_data_dir, "combined_mags.fa")
    try:
        with open(combined_mags_file, 'w') as outfile:
            for filename in os.listdir(corrected_mags_dir):
                if filename.endswith('.fa') or filename.endswith('.fasta'):
                    filepath = os.path.join(corrected_mags_dir, filename)
                    with open(filepath, 'r') as infile:
                        outfile.write(infile.read())
        print(f"[INFO] Successfully combined MAGs into: {combined_mags_file}")
    except Exception as e:
        print(f"[ERROR] Failed to combine MAGs: {e}")
        return False
    
    # Step 6: Add to Kraken library and build database
    print(f"[INFO] Step 6: Adding to Kraken library and building database")
    
    # Get Kraken path from config
    kraken_build_path = get_tool_path("kraken2_build")
    if not kraken_build_path:
        print(f"[ERROR] kraken2_build tool not found in configuration")
        return False
    
    # Create add-to-library script
    add_lib_script = os.path.join(scripts_dir, "5_add_to_library.sh")
    
    with open(add_lib_script, 'w') as f:
        f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {kraken_build_path} \\
  --add-to-library {combined_mags_file} \\
  --threads {cpus} \\
  --db {kraken_db_path}
""")
    
    # Create build database script
    build_db_script = os.path.join(scripts_dir, "6_build_database.sh")
    
    with open(build_db_script, 'w') as f:
        f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {kraken_build_path} \\
  --build \\
  --threads {cpus} \\
  --db {kraken_db_path}
""")
    
    # Make scripts executable
    os.chmod(add_lib_script, 0o755)
    os.chmod(build_db_script, 0o755)
    
    # Run add to library
    try:
        print(f"[INFO] Adding MAGs to Kraken library")
        if not run_command(add_lib_script):
            print("[ERROR] Failed to add MAGs to Kraken library")
            return False
        
        # Run build database
        print(f"[INFO] Building Kraken database")
        if not run_command(build_db_script):
            print("[ERROR] Failed to build Kraken database")
            return False
        
        print(f"[INFO] Kraken database built successfully in {kraken_db_path}")
        
        # Step 7: Build Bracken database
        print(f"[INFO] Step 7: Building Bracken database")
        
        # Get Bracken path from config or use default path
        bracken_build_path = get_tool_path("bracken_build")
        if not bracken_build_path:
            bracken_build_path = "/usr/home/workspace/tmp.0/zexi/app/Bracken-2.7/bracken-build"
            print(f"[WARNING] Using default bracken_build path: {bracken_build_path}")
        
        kraken2_path = get_tool_path("kraken2")
        if not kraken2_path:
            kraken2_path = "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/6_Taxonomy/kraken2-2.1.2"
            print(f"[WARNING] Using default kraken2 path: {kraken2_path}")
        
        # Create build Bracken database script
        build_bracken_script = os.path.join(scripts_dir, "7_build_bracken.sh")
        
        with open(build_bracken_script, 'w') as f:
            f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {bracken_build_path} \\
  -d {kraken_db_path} \\
  -t {cpus} \\
  -l {read_length} \\
  -x {kraken2_path}
""")
        
        # Make Bracken script executable
        os.chmod(build_bracken_script, 0o755)
        
        # Run build Bracken database
        print(f"[INFO] Building Bracken database with read length {read_length}")
        if not run_command(build_bracken_script):
            print("[WARNING] Bracken database build failed, but Kraken database is still usable")
        else:
            print(f"[INFO] Bracken database built successfully in {kraken_db_path}")
        
        # Step 8: Build Database Distribution
        print(f"[INFO] Step 8: Building database distribution")
        
        # Create database distribution script
        dist_script = os.path.join(scripts_dir, "8_build_distribution.sh")
        
        with open(dist_script, 'w') as f:
            f.write(f"""#!/bin/bash
#SBATCH -p ghpc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {cpus}
#SBATCH --mem 720G
#SBATCH -t 30-12:00:00
#SBATCH --output=combine_%db_build.out
#SBATCH --error=combine_%db_build.err

cd {kraken_db_path}

echo "[INFO] Building database distribution with Bracken"
echo "Database: {kraken_db_path}"
echo "CPUs: {cpus}"
echo "Read length: {read_length}"

PERL5LIB=/opt/ghpc/lib64/perl-5.36 {bracken_build_path} \\
  -d {kraken_db_path} \\
  -t {cpus} \\
  -l {read_length} \\
  -x {kraken2_path}

echo "[INFO] Database distribution build completed"
""")
        
        # Make script executable
        os.chmod(dist_script, 0o755)
        
        # Run database distribution build
        print(f"[INFO] Building database distribution with enhanced Bracken")
        if not run_command(f"bash {dist_script}"):
            print("[WARNING] Database distribution build failed, but core database is still usable")
        else:
            print(f"[INFO] Database distribution built successfully")
        
        print("=" * 70)
        print("? Complete Kraken + Bracken database with distribution built successfully")
        print(f"? Database location: {kraken_db_path}")
        print(f"? Advanced taxonomy resolver used for conflict resolution")
        print(f"? Conflicts log available at: {os.path.join(taxonomy_dir, 'conflicts.log')}")
        print("=" * 70)
        return True
        
    except Exception as e:
        print(f"[ERROR] Database build failed: {e}")
        return False


# Add the other functions with similar modifications...
def create_taxonomy_files(project_output, kraken_db_path, merge_mode=True, custom_gtdbtk_dir=None):
    """Create taxonomy files (nodes.dmp, names.dmp) only"""
    
    # Setup directories
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    processed_data_dir = os.path.join(kraken_db_path, "processed_data")
    taxonomy_dir = os.path.join(kraken_db_path, "taxonomy")
    
    for directory in [kraken_db_path, scripts_dir, processed_data_dir, taxonomy_dir]:
        ensure_directory_exists(directory)
    
    # Get existing taxonomy info from the specified kraken_db_path
    existing_nodes = os.path.join(taxonomy_dir, "nodes.dmp")
    existing_names = os.path.join(taxonomy_dir, "names.dmp")
    
    # Prepare taxonomy input file
    tax_tsv_noheader_file = os.path.join(processed_data_dir, "taxonomy_noheader.tsv")
    
    # Check if taxonomy input exists
    if not os.path.exists(tax_tsv_noheader_file):
        print("[ERROR] Taxonomy input file not found. Run previous steps first.")
        return False
    
    # Find and copy advanced resolver script
    advanced_script_path = find_script_in_main_directory("gtdb_to_taxdump_latest_resolve.py")
    if not advanced_script_path:
        print("[ERROR] Could not find gtdb_to_taxdump_latest_resolve.py script")
        return False
    
    advanced_resolver_script = os.path.join(scripts_dir, "3_gtdb_to_taxdump_latest_resolve.py")
    
    try:
        shutil.copy2(advanced_script_path, advanced_resolver_script)
        os.chmod(advanced_resolver_script, 0o755)
    except Exception as e:
        print(f"[ERROR] Failed to copy advanced resolver script: {e}")
        return False
    
    # Prepare resolver command
    if merge_mode and os.path.exists(existing_nodes):
        resolver_command = f"""python3 {advanced_resolver_script} \\
    --existing-nodes {existing_nodes} \\
    --existing-names {existing_names} \\
    --new-taxonomy {tax_tsv_noheader_file} \\
    --output-dir {taxonomy_dir}"""
    else:
        # Create minimal existing taxonomy
        temp_existing_dir = os.path.join(processed_data_dir, "temp_existing")
        ensure_directory_exists(temp_existing_dir)
        
        temp_nodes = os.path.join(temp_existing_dir, "nodes.dmp")
        temp_names = os.path.join(temp_existing_dir, "names.dmp")
        
        with open(temp_nodes, 'w') as f:
            f.write("1\t|\t1\t|\tno rank\t|\t\t|\n")
        
        with open(temp_names, 'w') as f:
            f.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        
        resolver_command = f"""python3 {advanced_resolver_script} \\
    --existing-nodes {temp_nodes} \\
    --existing-names {temp_names} \\
    --new-taxonomy {tax_tsv_noheader_file} \\
    --output-dir {taxonomy_dir}"""
    
    # Run taxonomy creation
    return run_command(resolver_command)

def rename_fasta_headers(project_output, kraken_db_path, cpus):
    """Rename FASTA headers for Kraken compatibility only"""
    
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    processed_data_dir = os.path.join(kraken_db_path, "processed_data")
    taxonomy_dir = os.path.join(kraken_db_path, "taxonomy")
    corrected_mags_dir = os.path.join(processed_data_dir, "mags_corrected_header")
    
    ensure_directory_exists(corrected_mags_dir)
    
    # Check required files
    taxid_map_file = os.path.join(taxonomy_dir, "taxid.map")
    if not os.path.exists(taxid_map_file):
        print("[ERROR] taxid.map not found. Run create_taxonomy step first.")
        return False
    
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    if not os.path.exists(mags_repo_dir):
        print(f"[ERROR] MAGs Repository directory not found: {mags_repo_dir}")
        return False
    
    # Find the header renaming script in the main script directory
    header_script_path = find_script_in_main_directory("header_kraken.py")
    if not header_script_path:
        print("[ERROR] Could not find header_kraken.py script")
        return False
    
    # Create modified header script with correct paths
    local_header_script = os.path.join(scripts_dir, "4_rename_headers.py")
    
    modified_header_content = f"""import os
import multiprocessing as mp
from Bio import SeqIO
from functools import partial

# Path to your MAG files folder
mags_folder_path = '{mags_repo_dir}'
# Path to your taxid.map file  
taxid_map_path = '{taxid_map_file}'
# Output folder where renamed files will be saved
output_folder = '{corrected_mags_dir}'

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Load the taxid.map file into a dictionary
def load_taxid_map(map_path):
    genome_to_taxid = {{}}
    
    with open(map_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) == 2:
                genome_id, taxid = parts
                genome_to_taxid[genome_id] = taxid
    
    return genome_to_taxid

# Function to rename headers based on the filename and taxid map
def rename_fasta_headers(fasta_file, mags_folder_path, output_folder, genome_to_taxid):
    input_file_path = os.path.join(mags_folder_path, fasta_file)
    output_file_path = os.path.join(output_folder, fasta_file)
    
    # Extract the genome name from the filename (remove extension)
    genome_name = fasta_file.replace('.fa', '')
    
    # Look up the taxid for this genome
    taxid = genome_to_taxid.get(genome_name, 'UnknownTaxID')
    
    try:
        with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
            # Initialize counter for contig numbers
            contig_counter = 1
            
            for record in SeqIO.parse(input_handle, 'fasta'):
                # Create the new standardized header with sequential contig numbers
                new_id = f"{{genome_name}}_{{contig_counter}}|kraken:taxid|{{taxid}}"
                contig_counter += 1
                
                # Update the record
                record.id = new_id
                record.description = ''
                
                # Write the renamed record
                SeqIO.write(record, output_handle, 'fasta')
        
        return f"Successfully processed {{fasta_file}} - renamed {{contig_counter-1}} contigs"
    except Exception as e:
        return f"Error processing {{fasta_file}}: {{str(e)}}"

def main():
    # Load taxid map
    genome_to_taxid = load_taxid_map(taxid_map_path)
    
    # Get list of FASTA files
    fasta_files = [f for f in os.listdir(mags_folder_path) if f.endswith('.fa')]
    
    # Set up the process pool
    num_cpus = {cpus}  # Use specified CPUs
    print(f"Using {{num_cpus}} CPUs for processing")
    
    # Create a partial function with the fixed arguments
    process_file = partial(
        rename_fasta_headers,
        mags_folder_path=mags_folder_path,
        output_folder=output_folder,
        genome_to_taxid=genome_to_taxid
    )
    
    # Process files in parallel
    with mp.Pool(processes=num_cpus) as pool:
        results = pool.map(process_file, fasta_files)
    
    # Print results
    for result in results:
        print(result)
    
    print("Renaming completed. The renamed files are saved in:", output_folder)

if __name__ == "__main__":
    main()
"""
    
    # Write the modified header script
    with open(local_header_script, 'w') as f:
        f.write(modified_header_content)
    
    return run_command(f"python {local_header_script}")

def add_to_kraken_library(kraken_db_path, cpus):
    """Add MAGs to Kraken library only"""
    
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    processed_data_dir = os.path.join(kraken_db_path, "processed_data")
    corrected_mags_dir = os.path.join(processed_data_dir, "mags_corrected_header")
    
    # Check if corrected MAGs exist
    if not os.path.exists(corrected_mags_dir):
        print("[ERROR] Corrected MAGs directory not found. Run rename_headers step first.")
        return False
    
    # Combine all MAGs into single file
    combined_mags_file = os.path.join(processed_data_dir, "combined_mags.fa")
    try:
        with open(combined_mags_file, 'w') as outfile:
            for filename in os.listdir(corrected_mags_dir):
                if filename.endswith('.fa') or filename.endswith('.fasta'):
                    filepath = os.path.join(corrected_mags_dir, filename)
                    with open(filepath, 'r') as infile:
                        outfile.write(infile.read())
    except Exception as e:
        print(f"[ERROR] Failed to combine MAGs: {e}")
        return False
    
    # Get Kraken build path
    kraken_build_path = get_tool_path("kraken2_build")
    if not kraken_build_path:
        print("[ERROR] kraken2_build tool not found in configuration")
        return False
    
    # Create and run add-to-library script
    add_lib_script = os.path.join(scripts_dir, "5_add_to_library.sh")
    
    with open(add_lib_script, 'w') as f:
        f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {kraken_build_path} \\
  --add-to-library {combined_mags_file} \\
  --threads {cpus} \\
  --db {kraken_db_path}
""")
    
    os.chmod(add_lib_script, 0o755)
    return run_command(add_lib_script)

def build_kraken_db(kraken_db_path, cpus):
    """Build main Kraken database only"""
    
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    
    # Get Kraken build path
    kraken_build_path = get_tool_path("kraken2_build")
    if not kraken_build_path:
        print("[ERROR] kraken2_build tool not found in configuration")
        return False
    
    # Create and run build database script
    build_db_script = os.path.join(scripts_dir, "6_build_database.sh")
    
    with open(build_db_script, 'w') as f:
        f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {kraken_build_path} \\
  --build \\
  --threads {cpus} \\
  --db {kraken_db_path}
""")
    
    os.chmod(build_db_script, 0o755)
    return run_command(build_db_script)

def build_bracken_db(kraken_db_path, cpus, read_length=150):
    """Build Bracken database only"""
    
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    
    # Get tool paths
    bracken_build_path = get_tool_path("bracken_build")
    if not bracken_build_path:
        bracken_build_path = "/usr/home/workspace/tmp.0/zexi/app/Bracken-2.7/bracken-build"
    
    kraken2_path = get_tool_path("kraken2")
    if not kraken2_path:
        kraken2_path = "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/6_Taxonomy/kraken2-2.1.2"
    
    # Create and run Bracken build script
    build_bracken_script = os.path.join(scripts_dir, "7_build_bracken.sh")
    
    with open(build_bracken_script, 'w') as f:
        f.write(f"""#!/bin/bash
PERL5LIB=/opt/ghpc/lib64/perl-5.36 {bracken_build_path} \\
  -d {kraken_db_path} \\
  -t {cpus} \\
  -l {read_length} \\
  -x {kraken2_path}
""")
    
    os.chmod(build_bracken_script, 0o755)
    return run_command(build_bracken_script)

def build_db_distribution(kraken_db_path, cpus, read_length=150):
    """Build database distribution only"""
    
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    
    # Get tool paths
    bracken_build_path = get_tool_path("bracken_build")
    if not bracken_build_path:
        bracken_build_path = "/usr/home/workspace/tmp.0/zexi/app/Bracken-2.7/bracken-build"
    
    kraken2_path = get_tool_path("kraken2")
    if not kraken2_path:
        kraken2_path = "/usr/home/ironbank/researchdata/metagenome/meta_pipeline/6_Taxonomy/kraken2-2.1.2"
    
    # Create database distribution script
    dist_script = os.path.join(scripts_dir, "8_build_distribution.sh")
    
    with open(dist_script, 'w') as f:
        f.write(f"""#!/bin/bash
cd {kraken_db_path}

echo "[INFO] Building database distribution with Bracken"
echo "Database: {kraken_db_path}"

PERL5LIB=/opt/ghpc/lib64/perl-5.36 {bracken_build_path} \\
  -d {kraken_db_path} \\
  -t {cpus} \\
  -l {read_length} \\
  -x {kraken2_path}

echo "[INFO] Database distribution build completed"
""")
    
    os.chmod(dist_script, 0o755)
    return run_command(f"bash {dist_script}")

def prepare_taxonomy_input(project_output, kraken_db_path, custom_gtdbtk_dir=None, rumen_ref_mags_dir=None):
    """Prepare taxonomy input files (Steps 1-2 combined)"""
    
    print("[INFO] Preparing taxonomy input files")
    
    # Setup directories
    scripts_dir = os.path.join(kraken_db_path, "scripts")
    processed_data_dir = os.path.join(kraken_db_path, "processed_data")
    
    for directory in [scripts_dir, processed_data_dir]:
        ensure_directory_exists(directory)
    
    # Define MAG sources
    mags_repo_dir = os.path.join(project_output, "MAGs_Repository")
    
    if not os.path.exists(mags_repo_dir):
        print(f"[ERROR] MAGs Repository directory not found: {mags_repo_dir}")
        return False
    
    # Get MAGs list
    repo_mags = [f for f in os.listdir(mags_repo_dir) if f.endswith(('.fa', '.fasta'))]
    
    if len(repo_mags) == 0:
        print(f"[WARNING] No MAGs found in repository: {mags_repo_dir}")
        return False
    
    all_mags = []
    for mag_file in repo_mags:
        all_mags.append(os.path.splitext(mag_file)[0])
    
    # Add rumen reference MAGs if they exist
    if rumen_ref_mags_dir and os.path.exists(rumen_ref_mags_dir):
        mags_subdir = os.path.join(rumen_ref_mags_dir, "mags")
        if os.path.exists(mags_subdir):
            ref_mags = [f for f in os.listdir(mags_subdir) if f.endswith(('.fa', '.fasta'))]
            for mag_file in ref_mags:
                all_mags.append(os.path.splitext(mag_file)[0])
    
    # Create input file
    tax_input_file = os.path.join(processed_data_dir, "mags_list.xlsx")
    df = pd.DataFrame({
        'user_genome': all_mags,
        'classification': ["d__Bacteria;p__;c__;o__;f__;g__;s__"] * len(all_mags)
    })
    df.to_excel(tax_input_file, index=False)
    
    # Run steps 1 and 2 from main function
    # Step 1: Placeholder taxa assignment
    if custom_gtdbtk_dir:
        gtdbtk_dir = custom_gtdbtk_dir
    else:
        gtdbtk_dir = os.path.join(project_output, "Novel_Mags", "gtdbtk")
    
    bac_summary = os.path.join(gtdbtk_dir, "classify", "gtdbtk.bac120.summary.tsv")
    ar_summary = os.path.join(gtdbtk_dir, "classify", "gtdbtk.ar53.summary.tsv")
    
    ref_taxonomy_dir = None
    if rumen_ref_mags_dir:
        ref_taxonomy_dir = os.path.join(rumen_ref_mags_dir, "taxonomy")
    
    placeholder_script = os.path.join(scripts_dir, "1_placeholder_taxa.py")
    tax_fixed_file = os.path.join(processed_data_dir, "tax_fixed.xlsx")
    
    # Create placeholder taxa script (same as main function)
    with open(placeholder_script, 'w') as f:
        f.write(f"""
import pandas as pd
import os

# Load the input data
input_file = '{tax_input_file}'
data = pd.read_excel(input_file)
print(f"Loaded {{len(data)}} rows from {{input_file}}")

# Function to clean data and replace all empty levels
def clean_and_replace_all_levels(row):
    classification = row['classification']
    genome = row['user_genome']
    
    if ';d__Bacteria;' in classification:
        parts = classification.split(';d__Bacteria;')
        if len(parts) >= 2:
            classification = 'd__Bacteria;' + parts[1]
    
    if not classification.startswith('d__'):
        if classification.startswith('p__'):
            classification = 'd__Bacteria;' + classification
    
    levels = classification.split(';')
    prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    level_names = ['d', 'p', 'c', 'o', 'f', 'g', 'sp']
    
    for i in range(min(len(levels), len(prefixes))):
        level = levels[i].strip()
        prefix = prefixes[i]
        level_name = level_names[i]
        
        if level == prefix:
            if prefix == 's__':
                levels[i] = f"{{prefix}}{{level_name}}{{genome}}*"
            else:
                levels[i] = f"{{prefix}}{{level_name}}{{genome}}*"
    
    return ';'.join(levels)

# Update with GTDB-Tk classifications
bac_path = '{bac_summary}'
ar_path = '{ar_summary}'

if os.path.exists(bac_path):
    bac_df = pd.read_csv(bac_path, sep='\\t')
    for _, row in bac_df.iterrows():
        idx = data[data['user_genome'] == row['user_genome']].index
        if len(idx) > 0:
            data.loc[idx, 'classification'] = row['classification']

if os.path.exists(ar_path):
    ar_df = pd.read_csv(ar_path, sep='\\t')
    for _, row in ar_df.iterrows():
        idx = data[data['user_genome'] == row['user_genome']].index
        if len(idx) > 0:
            data.loc[idx, 'classification'] = row['classification']

# Process reference taxonomy if available
ref_taxonomy_dir = '{ref_taxonomy_dir}'
if ref_taxonomy_dir and os.path.exists(ref_taxonomy_dir):
    for file in os.listdir(ref_taxonomy_dir):
        if file.endswith('.tsv') and not file.startswith('.'):
            try:
                tax_file = os.path.join(ref_taxonomy_dir, file)
                print(f"Loading reference taxonomy from: {{tax_file}}")
                ref_tax_df = pd.read_csv(tax_file, sep='\\t')
                
                if 'user_genome' in ref_tax_df.columns and 'classification' in ref_tax_df.columns:
                    for _, row in ref_tax_df.iterrows():
                        idx = data[data['user_genome'] == row['user_genome']].index
                        if len(idx) > 0:
                            data.loc[idx, 'classification'] = row['classification']
                            print(f"Updated taxonomy for {{row['user_genome']}}")
            except Exception as e:
                print(f"Error loading reference taxonomy file {{file}}: {{e}}")

# Clean and save
clean_data = pd.DataFrame()
clean_data['user_genome'] = data['user_genome']
clean_data['classification'] = data.apply(clean_and_replace_all_levels, axis=1)

output_file = '{tax_fixed_file}'
clean_data.to_excel(output_file, index=False)
print(f"Saved cleaned data to {{output_file}}")
""")
    
    # Run placeholder taxa script
    if not run_command(f"python {placeholder_script}"):
        print("[ERROR] Failed to run placeholder taxa assignment")
        return False
    
    # Step 2: Convert to TSV
    tsv_conversion_script = os.path.join(scripts_dir, "2_convert_to_tsv.py")
    tax_tsv_noheader_file = os.path.join(processed_data_dir, "taxonomy_noheader.tsv")
    
    with open(tsv_conversion_script, 'w') as f:
        f.write(f"""
import pandas as pd

input_file = '{tax_fixed_file}'
noheader_file = '{tax_tsv_noheader_file}'

data = pd.read_excel(input_file)
data.to_csv(noheader_file, sep='\\t', index=False, header=False)
print(f"Conversion complete! TSV file saved to {{noheader_file}}")
""")
    
    # Run TSV conversion
    if not run_command(f"python {tsv_conversion_script}"):
        print("[ERROR] Failed to run TSV conversion")
        return False
    
    print("[INFO] Taxonomy input preparation completed successfully")
    return True
#!/usr/bin/env python3
"""
Module for running dbCAN on MAGs.
This module handles the annotation of MAGs using dbCAN for CAZyme identification.
"""

import os
import subprocess
import yaml
import shutil
import time
import glob
from MetaMAG.config import config, get_tool_path, get_tool_command

def run(cpus, memory, time_limit, project_input, project_output, mags_dir=None):
    """
    Run dbCAN on the provided MAGs directory.
    
    Parameters:
    -----------
    cpus : int
        Number of CPUs to use for annotation
    memory : str
        Memory allocation (e.g., "200G")
    time_limit : str
        Time allocation (e.g., "30-12:00:00")
    project_input : str
        Path to project input directory
    project_output : str
        Path to project output directory
    mags_dir : str, optional
        Path to MAGs directory. If None, uses default from project configuration.
    """
    print(f"[INFO] Starting dbCAN annotation for CAZyme identification")
    print(f"[INFO] Using {cpus} CPUs")
    
    # Set up input and output directories
    if not mags_dir:
        mags_dir = os.path.join(project_output, "final_mags")
    
    if not os.path.exists(mags_dir):
        print(f"[ERROR] MAGs directory not found: {mags_dir}")
        return

    # Create output directory
    dbcan_output_dir = os.path.join(project_output, "dbcan_output")
    os.makedirs(dbcan_output_dir, exist_ok=True)
    
    # Get dbCAN and Prodigal paths from config
    dbcan_path = get_tool_path("dbcan")
    if not dbcan_path:
        print(f"[ERROR] Dbcan tool not found in configuration")
        return
    dbcan_db_dir = get_tool_path("dbcan_db_dir")
    if not dbcan_db_dir:
        print(f"[ERROR] Dbcan_Db_Dir tool not found in configuration")
        return
    prodigal_path = get_tool_path("prodigal")
    if not prodigal_path:
        print(f"[ERROR] Prodigal tool not found in configuration")
        return  # Default to "prodigal" if not specified
    
    if not dbcan_path:
        print(f"[ERROR] dbCAN path not found in config. Make sure 'dbcan' is defined in the tools section.")
        return
    
    if not os.path.exists(dbcan_path):
        print(f"[ERROR] dbCAN path does not exist: {dbcan_path}")
        return
    
    if not dbcan_db_dir:
        print(f"[ERROR] dbCAN database directory not found in config. Make sure 'dbcan_db_dir' is defined in the tools section.")
        return
    
    if not os.path.exists(dbcan_db_dir):
        print(f"[ERROR] dbCAN database directory does not exist: {dbcan_db_dir}")
        return

    # Find all MAG files
    mag_files = []
    for root, dirs, files in os.walk(mags_dir):
        for file in files:
            if file.endswith((".fa", ".fasta", ".fna")):
                mag_files.append(os.path.join(root, file))
    
    # Check if all MAGs have already been processed
    need_processing = False
    mags_to_process = []
    
    for mag_file in mag_files:
        mag_base_name = os.path.basename(mag_file).replace(".", "_")
        
        # Check if dbCAN result directory exists
        dbcan_result_dir = os.path.join(dbcan_output_dir, f"{mag_base_name}_dbcan_results")
        
        if not os.path.exists(dbcan_result_dir) or not os.listdir(dbcan_result_dir):
            need_processing = True
            mags_to_process.append(mag_file)
    
    if not need_processing:
        print(f"[INFO] All {len(mag_files)} MAGs already have dbCAN results. Skipping annotation step.")
        print(f"[INFO] Summarizing CAZyme annotations from existing files.")
        cazyme_dir = summarize_cazyme_annotations(dbcan_output_dir, project_output)
        
        # Create visualizations
        try:
            print(f"[INFO] Creating CAZyme visualization plots")
            plots_dir = create_visualization_plots(cazyme_dir)
            if plots_dir:
                print(f"[INFO] Visualization plots created and available in {plots_dir}")
        except Exception as e:
            print(f"[WARNING] Failed to create visualization plots: {e}")
            print(f"[WARNING] Make sure matplotlib, seaborn, and pandas are installed")
            import traceback
            traceback.print_exc()
            
        return
    else:
        print(f"[INFO] Found {len(mags_to_process)} MAGs that need processing out of {len(mag_files)} total MAGs.")
    
    # Step 1: First look for protein files for MAGs that need processing
    protein_files_map = {}  # Map MAG file to its protein file

    # Check for existing protein files in MAGs directory
    for mag_file in mags_to_process:
        mag_base = os.path.basename(mag_file)
        mag_dir = os.path.dirname(mag_file)
        
        # Common naming patterns for protein files
        potential_protein_names = [
            f"{mag_base}.faa",
            f"{mag_base}.proteins.fa",
            f"{mag_base}.pep",
            f"{os.path.splitext(mag_base)[0]}.faa",
            f"{os.path.splitext(mag_base)[0]}.proteins.fa"
        ]
        
        # Check MAG directory for protein files
        protein_found = False
        for protein_name in potential_protein_names:
            protein_path = os.path.join(mag_dir, protein_name)
            if os.path.exists(protein_path):
                protein_files_map[mag_file] = protein_path
                protein_found = True
                break
        
        # Check output directory for existing protein files
        if not protein_found:
            for protein_name in potential_protein_names:
                protein_path = os.path.join(dbcan_output_dir, protein_name)
                if os.path.exists(protein_path):
                    protein_files_map[mag_file] = protein_path
                    protein_found = True
                    break
    
    # Step 2: For MAGs without protein files, predict them using Prodigal
    for mag_file in mags_to_process:
        if mag_file not in protein_files_map:
            protein_file = os.path.join(dbcan_output_dir, os.path.basename(mag_file) + ".faa")
            
            # Skip if protein prediction already exists
            if os.path.exists(protein_file):
                print(f"[INFO] Protein file already exists: {protein_file}")
                protein_files_map[mag_file] = protein_file
                continue
                
            print(f"[INFO] Predicting proteins from {mag_file}")
            prodigal_cmd = f"{prodigal_path} -i {mag_file} -a {protein_file} -p meta"
            subprocess.run(prodigal_cmd, shell=True, check=True)
            protein_files_map[mag_file] = protein_file
    
    # Step 3: Run dbCAN on each protein file
    for mag_file, protein_file in protein_files_map.items():
        # Get mag base name for output directory
        mag_base_name = os.path.basename(mag_file).replace(".", "_")
        
        # Define output directory for this MAG
        mag_dbcan_out_dir = os.path.join(dbcan_output_dir, f"{mag_base_name}_dbcan_results")
        
        # Skip if output already exists
        if os.path.exists(mag_dbcan_out_dir) and os.listdir(mag_dbcan_out_dir):
            print(f"[INFO] dbCAN results already exist for {mag_base_name}. Skipping.")
            continue
            
        # Create output directory
        os.makedirs(mag_dbcan_out_dir, exist_ok=True)
        
        # Run dbCAN
        print(f"[INFO] Processing {protein_file} with dbCAN")
        
        # Get conda activation command from config
        dbcan_activate = get_tool_path("dbcan_activate")
        
        # Build environment and command
        env = os.environ.copy()
        dbcan_env_path = "/usr/home/qgg/maralta/.local/bin/miniconda3/envs/dbcan"
        
        # Add dbCAN environment to PATH
        if os.path.exists(dbcan_env_path):
            env['PATH'] = f"{dbcan_env_path}/bin:{env.get('PATH', '')}"
            env['CONDA_PREFIX'] = dbcan_env_path
            env['CONDA_DEFAULT_ENV'] = 'dbcan'
        
        if dbcan_activate:
            # Use configured activation command
            dbcan_cmd = f"source {dbcan_activate} && {dbcan_path} {protein_file} protein --db_dir {dbcan_db_dir} --out_dir {mag_dbcan_out_dir} --dia_cpu {cpus} --hmm_cpu {cpus} -t hmmer"
        else:
            # Fallback: try to activate manually
            conda_sh = "/usr/home/qgg/maralta/.local/bin/miniconda3/etc/profile.d/conda.sh"
            if os.path.exists(conda_sh):
                dbcan_cmd = f"source {conda_sh} && conda activate dbcan && {dbcan_path} {protein_file} protein --db_dir {dbcan_db_dir} --out_dir {mag_dbcan_out_dir} --dia_cpu {cpus} --hmm_cpu {cpus} -t hmmer"
            else:
                # Last resort: run with modified environment
                dbcan_cmd = f"{dbcan_path} {protein_file} protein --db_dir {dbcan_db_dir} --out_dir {mag_dbcan_out_dir} --dia_cpu {cpus} --hmm_cpu {cpus} -t hmmer"
        
        print(f"[INFO] Running command: {dbcan_cmd}")
        
        try:
            subprocess.run(dbcan_cmd, shell=True, check=True, executable="/bin/bash", env=env)
            print(f"[INFO] Successfully processed {mag_base_name}")
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to process {mag_base_name}: {e}")
            
            # Check if the issue is missing dependencies
            hmmscan_path = os.path.join(dbcan_env_path, "bin", "hmmscan")
            if not os.path.exists(hmmscan_path):
                print(f"[ERROR] hmmscan not found at {hmmscan_path}")
                print(f"[ERROR] Please ensure the dbCAN conda environment is properly installed with all dependencies")
                print(f"[ERROR] Try: conda activate dbcan && conda install hmmer")
            else:
                print(f"[ERROR] dbCAN execution failed despite hmmscan being available")
            
            continue
        
    print(f"[INFO] dbCAN annotation completed. Results available in {dbcan_output_dir}")
    
    # Summarize CAZyme annotations
    print(f"[INFO] Summarizing CAZyme annotations")
    cazyme_dir = summarize_cazyme_annotations(dbcan_output_dir, project_output)
    
    # Create visualizations
    try:
        print(f"[INFO] Creating CAZyme visualization plots")
        plots_dir = create_visualization_plots(cazyme_dir)
        if plots_dir:
            print(f"[INFO] Visualization plots created and available in {plots_dir}")
    except Exception as e:
        print(f"[WARNING] Failed to create visualization plots: {e}")
        print(f"[WARNING] Make sure matplotlib, seaborn, and pandas are installed")
        import traceback
        traceback.print_exc()

def summarize_cazyme_annotations(dbcan_output_dir, project_output):
    """
    Summarize CAZyme annotations from dbCAN output.
    Create counts of CAZyme families and categories for each MAG.
    
    Parameters:
    -----------
    dbcan_output_dir : str
        Path to directory containing dbCAN output files
    project_output : str
        Path to project output directory
    """
    # Create output directories
    cazyme_dir = os.path.join(project_output, "cazyme_annotations")
    os.makedirs(cazyme_dir, exist_ok=True)
    
    cazyme_counts_dir = os.path.join(cazyme_dir, "cazyme_counts")
    os.makedirs(cazyme_counts_dir, exist_ok=True)
    
    # Process each dbCAN result directory
    mag_count = 0
    for result_dir_name in os.listdir(dbcan_output_dir):
        if not result_dir_name.endswith("_dbcan_results"):
            continue
            
        result_dir = os.path.join(dbcan_output_dir, result_dir_name)
        if not os.path.isdir(result_dir):
            continue
            
        mag_base_name = result_dir_name.replace("_dbcan_results", "")
        
        # Check for the overview.txt file
        overview_file = os.path.join(result_dir, "overview.txt")
        if not os.path.exists(overview_file):
            print(f"[WARNING] overview.txt not found for {mag_base_name}")
            continue
            
        # Create a counts file for this MAG
        cazyme_counts_file = os.path.join(cazyme_counts_dir, f"{mag_base_name}.counts")
        
        # Skip if counts file already exists and is not empty
        if os.path.exists(cazyme_counts_file) and os.path.getsize(cazyme_counts_file) > 0:
            mag_count += 1
            continue
            
        # Parse the overview.txt file to extract CAZymes
        cazymes = {}
        categories = {"GT": 0, "CE": 0, "GH": 0, "CBM": 0, "PL": 0, "AA": 0}
        
        try:
            with open(overview_file, 'r') as f:
                lines = f.readlines()
                print(f"[DEBUG] Reading overview.txt for {mag_base_name}")
                print(f"[DEBUG] File has {len(lines)} lines")
                
                for line_num, line in enumerate(lines):
                    # Skip header line and empty lines
                    if line_num == 0 or line.startswith("#") or line.strip() == '':
                        continue
                        
                    parts = line.strip().split('\t')
                    
                    if len(parts) < 2:
                        continue
                    
                    # Debug first few data lines
                    if line_num < 10:
                        print(f"[DEBUG] Line {line_num}: {line.strip()}")
                        print(f"[DEBUG] Parts count: {len(parts)}")
                    
                    # CAZyme info is in column 1 (HMMER column)
                    hmmer_hit = parts[1]
                    
                    if hmmer_hit == '-':
                        continue
                    
                    # Debug CAZyme extraction
                    if line_num < 10:
                        print(f"[DEBUG] Found CAZyme: {hmmer_hit}")
                    
                    # Parse all CAZy families identified by HMMER
                    # Handle multiple families separated by +
                    for hit in hmmer_hit.split('+'):
                        # Extract just the family name (before parentheses)
                        family = hit.split('(')[0].strip()
                        
                        # Debug the family extraction
                        if line_num < 10:
                            print(f"[DEBUG] Extracted family: {family} from hit: {hit}")
                        
                        # Update category counts
                        category_found = False
                        for cat in categories.keys():
                            if family.startswith(cat):
                                categories[cat] += 1
                                category_found = True
                                break
                        
                        if not category_found:
                            print(f"[DEBUG] Unknown CAZyme family: {family}")
                                
                        # Update family counts
                        if family in cazymes:
                            cazymes[family] += 1
                        else:
                            cazymes[family] = 1
                            
            # Calculate total count
            total = sum(categories.values())
            categories["Total"] = total
            
            # Debug final counts
            print(f"[DEBUG] For {mag_base_name}:")
            print(f"[DEBUG] - Total CAZymes found: {total}")
            print(f"[DEBUG] - Categories: {categories}")
            print(f"[DEBUG] - Number of unique families: {len(cazymes)}")
            
            # Write cazyme counts to file
            with open(cazyme_counts_file, 'w') as f:
                # Write CAZyme family counts
                for family, count in sorted(cazymes.items(), key=lambda x: x[0]):
                    f.write(f"{count:5d} {family}\n")
                
                # Write category counts
                f.write('\n')
                for category, count in categories.items():
                    f.write(f"{category} {count}\n")
                    
            mag_count += 1
            
        except Exception as e:
            print(f"[ERROR] Failed to process {overview_file}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"[INFO] CAZyme counts generated for {mag_count} MAGs. Results available in {cazyme_counts_dir}")
    
    # Generate combined CAZyme reports
    try:
        generate_cazyme_reports(cazyme_counts_dir, cazyme_dir)
    except Exception as e:
        print(f"[ERROR] Failed to generate CAZyme reports: {e}")
        
    # Return cazyme_dir for visualization functions
    return cazyme_dir
def generate_cazyme_reports(cazyme_counts_dir, cazyme_dir):
    """
    Generate combined CAZyme reports from individual MAG count files.
    Creates Excel files with counts of GH families and CAZyme categories.
    
    Parameters:
    -----------
    cazyme_counts_dir : str
        Path to directory containing CAZyme count files
    cazyme_dir : str
        Path to directory for storing output reports
    """
    try:
        import pandas as pd
    except ImportError:
        print("[ERROR] Pandas is required for generating combined reports. Please install it with 'pip install pandas'.")
        return
    
    # Initialize dictionaries to store data
    gh_data = {}
    category_data = {}
    all_cazyme_classes = set()
    categories = ['GT', 'CE', 'GH', 'CBM', 'PL', 'AA', 'Total']

    # Read each count file
    for filename in os.listdir(cazyme_counts_dir):
        file_path = os.path.join(cazyme_counts_dir, filename)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            
            # Find the empty line separating CAZyme counts from category counts
            separator_idx = next((i for i, line in enumerate(lines) if line.strip() == ''), len(lines))
            
            # Process CAZyme family counts
            cazyme_counts = {}
            for line in lines[:separator_idx]:
                parts = line.strip().split()
                if len(parts) >= 2:
                    count = int(parts[0])
                    cazyme_class = parts[1]
                    cazyme_counts[cazyme_class] = count
                    all_cazyme_classes.add(cazyme_class)
            
            # Process category counts
            cat_counts = {}
            for line in lines[separator_idx+1:]:
                parts = line.strip().split()
                if len(parts) >= 2:
                    category = parts[0]
                    count = int(parts[1])
                    cat_counts[category] = count
            
            gh_data[filename] = cazyme_counts
            category_data[filename] = cat_counts

    # Create DataFrame for all CAZyme counts
    cazyme_df = pd.DataFrame({filename: {cazy_class: gh_data[filename].get(cazy_class, 0) 
                                        for cazy_class in all_cazyme_classes} 
                             for filename in gh_data})

    # Create DataFrame for category counts
    category_df = pd.DataFrame(category_data).T

    # Create GH-specific DataFrame
    gh_classes = [cls for cls in all_cazyme_classes if cls.startswith('GH')]
    if gh_classes:
        gh_df = cazyme_df.loc[gh_classes]
        
        # Write to Excel files
        cazyme_df.to_excel(os.path.join(cazyme_dir, 'all_cazyme_counts.xlsx'))
        category_df.to_excel(os.path.join(cazyme_dir, 'category_counts.xlsx'))
        gh_df.to_excel(os.path.join(cazyme_dir, 'gh_cazyme_counts.xlsx'))
        
        # Group and sum similar GH classes
        try:
            grouped_gh_df = gh_df.groupby(gh_df.index.str.split('_').str[0]).sum()
            grouped_gh_df.to_excel(os.path.join(cazyme_dir, 'grouped_gh_counts.xlsx'))
        except Exception as e:
            print(f"[WARNING] Could not create grouped GH counts: {e}")
    else:
        # No GH classes found
        print("[WARNING] No GH families found in the data")
        cazyme_df.to_excel(os.path.join(cazyme_dir, 'all_cazyme_counts.xlsx'))
        category_df.to_excel(os.path.join(cazyme_dir, 'category_counts.xlsx'))
    
    print(f"[INFO] Combined CAZyme reports generated in {cazyme_dir}")

def create_visualization_plots(cazyme_dir):
    """
    Create visualization plots for CAZyme data.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results
        
    Returns:
    --------
    str
        Path to the plots directory
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        import numpy as np
    except ImportError as e:
        print(f"[ERROR] Required packages for visualization not found: {e}")
        print("[ERROR] Please install missing packages with: pip install matplotlib seaborn pandas numpy")
        return None
    
    # Create plots directory
    plots_dir = os.path.join(cazyme_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    print(f"[INFO] Creating visualization plots in {plots_dir}")
    
    # Load necessary data files
    try:
        gh_file = os.path.join(cazyme_dir, 'gh_cazyme_counts.xlsx')
        cat_file = os.path.join(cazyme_dir, 'category_counts.xlsx')
        
        # Check if files exist
        if not os.path.exists(cat_file):
            print(f"[ERROR] Category counts file not found: {cat_file}")
            return None
            
        cat_df = pd.read_excel(cat_file, index_col=0)
        
        # Check if dataframe is empty
        if cat_df.empty:
            print(f"[WARNING] Category counts dataframe is empty. No data to visualize.")
            return plots_dir
        
        # Load or create GH dataframe
        if not os.path.exists(gh_file):
            print(f"[WARNING] GH CAZyme counts file not found: {gh_file}")
            if os.path.exists(os.path.join(cazyme_dir, 'all_cazyme_counts.xlsx')):
                print(f"[INFO] Using all_cazyme_counts.xlsx to extract GH families")
                all_df = pd.read_excel(os.path.join(cazyme_dir, 'all_cazyme_counts.xlsx'), index_col=0)
                
                # Filter GH families
                gh_rows = [idx for idx in all_df.index if str(idx).startswith('GH')]
                if gh_rows:
                    gh_df = all_df.loc[gh_rows]
                    gh_df.to_excel(gh_file)
                else:
                    print(f"[WARNING] No GH families found in the data")
                    gh_df = pd.DataFrame()
            else:
                print(f"[ERROR] No CAZyme data files found. Cannot create visualizations.")
                return None
        else:
            gh_df = pd.read_excel(gh_file, index_col=0)
        
        # Debug information
        print(f"[DEBUG] Category DataFrame shape: {cat_df.shape}")
        print(f"[DEBUG] GH DataFrame shape: {gh_df.shape}")
        
        # 1. Create GH heatmap (only if we have data)
        if not gh_df.empty:
            plt.figure(figsize=(16, 10))
            sns.heatmap(gh_df, annot=False, cmap="YlGnBu")
            plt.title('Heatmap of GH Families Across MAGs', fontsize=14)
            plt.xlabel('MAGs', fontsize=12)
            plt.ylabel('GH Family', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'gh_families_heatmap.pdf'))
            plt.savefig(os.path.join(plots_dir, 'gh_families_heatmap.png'), dpi=300)
            plt.close()
        else:
            print(f"[WARNING] No GH families found for heatmap visualization")
        
        # 2. Create category barplot
        if 'Total' in cat_df.columns:
            cat_df = cat_df.drop('Total', axis=1)
            
        if not cat_df.empty:
            plt.figure(figsize=(12, 8))
            cat_df.plot(kind='bar', stacked=True, figsize=(12, 8), colormap='viridis')
            plt.title('CAZyme Category Distribution by MAG', fontsize=14)
            plt.xlabel('MAG', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.legend(title='CAZyme Category')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'cazyme_category_distribution.pdf'))
            plt.savefig(os.path.join(plots_dir, 'cazyme_category_distribution.png'), dpi=300)
            plt.close()
        
        # 3. Create average category bar chart
        if not cat_df.empty:
            category_avg = cat_df.mean()
            if not category_avg.empty:
                plt.figure(figsize=(10, 6))
                category_avg.plot(kind='bar', figsize=(10, 6), color=sns.color_palette('viridis', len(category_avg)))
                plt.title('Average CAZyme Category Counts Across MAGs', fontsize=14)
                plt.xlabel('CAZyme Category', fontsize=12)
                plt.ylabel('Average Count', fontsize=12)
                for i, v in enumerate(category_avg):
                    plt.text(i, v + 0.5, f'{v:.1f}', ha='center')
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, 'average_cazyme_categories.pdf'))
                plt.savefig(os.path.join(plots_dir, 'average_cazyme_categories.png'), dpi=300)
                plt.close()
        
        # 4. Create GH family distribution plot (only if we have GH data)
        if not gh_df.empty:
            gh_totals = gh_df.sum(axis=1).sort_values(ascending=False)
            top_n = min(20, len(gh_totals))
            if top_n > 0:
                top_gh_families = gh_totals.head(top_n)
                
                plt.figure(figsize=(14, 8))
                top_gh_families.plot(kind='bar', figsize=(14, 8), color=sns.color_palette('viridis', len(top_gh_families)))
                plt.title(f'Top {top_n} GH Families Across All MAGs', fontsize=14)
                plt.xlabel('GH Family', fontsize=12)
                plt.ylabel('Total Count', fontsize=12)
                for i, v in enumerate(top_gh_families):
                    plt.text(i, v + 0.5, str(v), ha='center')
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, f'top_{top_n}_gh_families.pdf'))
                plt.savefig(os.path.join(plots_dir, f'top_{top_n}_gh_families.png'), dpi=300)
                plt.close()
        
        # 5. Create prevalence plot (only if we have GH data)
        if not gh_df.empty:
            prevalence = (gh_df > 0).mean(axis=1).sort_values(ascending=False)
            top_n = min(20, len(prevalence))
            if top_n > 0:
                top_prevalent_gh = prevalence.head(top_n)
                
                plt.figure(figsize=(14, 8))
                top_prevalent_gh.plot(kind='bar', figsize=(14, 8), color=sns.color_palette('viridis', len(top_prevalent_gh)))
                plt.title(f'Top {top_n} Most Prevalent GH Families', fontsize=14)
                plt.xlabel('GH Family', fontsize=12)
                plt.ylabel('Fraction of MAGs', fontsize=12)
                for i, v in enumerate(top_prevalent_gh):
                    plt.text(i, v + 0.02, f'{v:.2f}', ha='center')
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, f'top_{top_n}_prevalent_gh_families.pdf'))
                plt.savefig(os.path.join(plots_dir, f'top_{top_n}_prevalent_gh_families.png'), dpi=300)
                plt.close()
            
        # 6. Create dendrogram clustering (only if we have enough data)
        if not gh_df.empty and gh_df.shape[1] > 2 and gh_df.shape[0] > 2:  # Need at least 3 MAGs and 3 CAZymes
            try:
                # Transpose so MAGs are rows and CAZymes are columns
                gh_df_t = gh_df.T
                
                plt.figure(figsize=(16, 10))
                g = sns.clustermap(gh_df_t, method='ward', metric='euclidean', 
                                   cmap="YlGnBu", figsize=(16, 10),
                                   xticklabels=True, yticklabels=True)
                plt.savefig(os.path.join(plots_dir, 'cazyme_dendrogram.pdf'))
                plt.savefig(os.path.join(plots_dir, 'cazyme_dendrogram.png'), dpi=300)
                plt.close()
            except Exception as e:
                print(f"[WARNING] Could not create dendrogram clustering: {e}")
        
        print(f"[INFO] Visualization plots created successfully in {plots_dir}")
        return plots_dir
        
    except Exception as e:
        print(f"[ERROR] Failed to create plots: {e}")
        import traceback
        traceback.print_exc()
        return None

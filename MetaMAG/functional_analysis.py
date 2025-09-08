#!/usr/bin/env python3
"""
Enhanced Functional Annotation Analysis Pipeline

This script combines and enhances functionality from functional_analysis.py and 
species_all_functional_annotation_analysis.py to create a comprehensive suite of
visualizations for functional annotation data.
"""

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import Counter, defaultdict
import glob
import logging
import traceback
import argparse
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.path as mpath
import matplotlib.patches as mpatches

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='[%(asctime)s] [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# COG category descriptions
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'X': 'Mobilome: prophages, transposons',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown'
}

# COG functional categories
COG_MAIN_CATEGORIES = {
    "INFORMATION STORAGE AND PROCESSING": ["A", "B", "J", "K", "L"],
    "CELLULAR PROCESSES AND SIGNALING": ["D", "M", "N", "O", "T", "U", "Y", "Z", "V", "W"],
    "METABOLISM": ["C", "E", "F", "G", "H", "I", "P", "Q"],
    "POORLY CHARACTERIZED": ["R", "S", "X"]
}

# COG category colors (matching common conventions)
COG_COLORS = {
    "A": "#E41A1C",  # Red
    "B": "#E41A1C",  # Red
    "J": "#E41A1C",  # Red
    "K": "#E41A1C",  # Red
    "L": "#E41A1C",  # Red
    
    "D": "#377EB8",  # Blue
    "M": "#377EB8",  # Blue
    "N": "#377EB8",  # Blue
    "O": "#377EB8",  # Blue
    "T": "#377EB8",  # Blue
    "U": "#377EB8",  # Blue
    "Y": "#377EB8",  # Blue
    "Z": "#377EB8",  # Blue
    "V": "#377EB8",  # Blue
    "W": "#377EB8",  # Blue
    
    "C": "#4DAF4A",  # Green
    "E": "#4DAF4A",  # Green
    "F": "#4DAF4A",  # Green
    "G": "#4DAF4A",  # Green
    "H": "#4DAF4A",  # Green
    "I": "#4DAF4A",  # Green
    "P": "#4DAF4A",  # Green
    "Q": "#4DAF4A",  # Green
    
    "R": "#984EA3",  # Purple
    "S": "#984EA3",  # Purple
    "X": "#984EA3"   # Purple
}
def run_functional_analysis(project_output, kegg_db_file=None, species_based=False):
    """
    Main function to run the functional analysis on EggNOG-mapper results.
    
    Parameters:
    -----------
    project_output : str
        Path to the project output directory
    kegg_db_file : str, optional
        Path to the KEGG database file (ko.txt)
    species_based : bool, optional
        Whether to perform species-based analysis
    """
    logger.info("Starting functional annotation analysis")
    
    # Call the enhanced version with all parameters passed through
    return run_enhanced_functional_analysis(project_output, kegg_db_file, species_based)

def create_folder_structure(analysis_dir):
    """
    Create the folder structure for organizing files and plots.
    
    Parameters:
    -----------
    analysis_dir : str
        Path to the analysis directory
    
    Returns:
    --------
    dict
        Dictionary containing paths to the created folders
    """
    # Define folder structure
    folders = {
        "files": {
            "ko": os.path.join(analysis_dir, "files", "ko"),
            "cog": os.path.join(analysis_dir, "files", "cog"),
            "kegg": os.path.join(analysis_dir, "files", "kegg"),
            "species": os.path.join(analysis_dir, "files", "species"),
            "network": os.path.join(analysis_dir, "files", "network")
        },
        "plots": {
            "ko": os.path.join(analysis_dir, "plots", "ko"),
            "cog": os.path.join(analysis_dir, "plots", "cog"),
            "kegg": os.path.join(analysis_dir, "plots", "kegg"),
            "species": os.path.join(analysis_dir, "plots", "species"),
            "network": os.path.join(analysis_dir, "plots", "network"),
            "sankey": os.path.join(analysis_dir, "plots", "sankey"),
            "heatmap": os.path.join(analysis_dir, "plots", "heatmap"),
            "circular": os.path.join(analysis_dir, "plots", "circular"),
            "combined": os.path.join(analysis_dir, "plots", "combined")
        },
        "reports": os.path.join(analysis_dir, "reports")
    }
    
    # Create the folders
    for category in folders:
        if isinstance(folders[category], dict):
            for subcategory in folders[category]:
                os.makedirs(folders[category][subcategory], exist_ok=True)
                logger.info(f"Created directory: {folders[category][subcategory]}")
        else:
            os.makedirs(folders[category], exist_ok=True)
            logger.info(f"Created directory: {folders[category]}")
    
    return folders

def count_eggnog_terms(input_dir, analysis_dir, folders):
    """
    Count KO and COG terms from EggNOG annotation files.
    
    Parameters:
    -----------
    input_dir : str
        Path to directory containing EggNOG annotation files
    analysis_dir : str
        Path to analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    
    Returns:
    --------
    tuple
        DataFrames containing COG counts, individual COG matrix, and KO counts
    """
    cog_counts_df = pd.DataFrame()
    ko_counts_df = pd.DataFrame()
    individual_cog_matrix = pd.DataFrame()  # For storing individual COG letters
    
    # Iterate over annotation files in the input folder
    annotation_files = glob.glob(os.path.join(input_dir, "*.emapper.annotations"))
    
    if not annotation_files:
        logger.warning(f"No annotation files found in {input_dir}")
        return cog_counts_df, individual_cog_matrix, ko_counts_df
    
    logger.info(f"Found {len(annotation_files)} annotation files")
    
    for annotation_file in annotation_files:
        try:
            base_name = os.path.basename(annotation_file).split('.emapper.annotations')[0]
            logger.info(f"Processing {base_name}")
            
            # Read the annotation file, skipping comment lines
            with open(annotation_file, 'r') as f:
                lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
            
            if not lines:
                logger.warning(f"No data found in {annotation_file}")
                continue
            
            # Process each line manually
            cog_counter = Counter()
            ko_counter = Counter()
            multi_cog_counter = Counter()  # For storing multi-letter COGs
            
            for line in lines:
                parts = line.split('\t')
                
                # Make sure we have enough columns
                if len(parts) < 12:
                    continue
                
                # Get COG and KO columns (7th and 12th columns, 0-indexed)
                cog_column = 6
                ko_column = 11
                
                # COG categories
                if cog_column < len(parts) and parts[cog_column] != '-' and parts[cog_column]:
                    cogs = parts[cog_column].strip()
                    
                    # Store multi-letter COGs
                    if len(cogs) > 1:
                        multi_cog_counter[cogs] += 1
                    
                    # Count each individual COG letter separately
                    for cog in cogs:
                        if cog.isalpha():  # Only count alphabet characters
                            cog_counter[cog] += 1
                
                # KO terms
                if ko_column < len(parts) and parts[ko_column] != '-' and parts[ko_column]:
                    kos = parts[ko_column].strip().split(',')
                    for ko in kos:
                        ko = ko.strip()
                        if ko.startswith('ko:'):
                            ko = ko[3:]  # Remove 'ko:' prefix
                        ko_counter[ko] += 1
            
            # Add counts to respective dataframes
            cog_counts_df[base_name] = pd.Series(multi_cog_counter)  # Store multi-letter COGs
            individual_cog_matrix[base_name] = pd.Series(cog_counter)  # Store individual COG letters
            ko_counts_df[base_name] = pd.Series(ko_counter)
            
            logger.info(f"Found {len(cog_counter)} COG categories and {len(ko_counter)} KO terms in {base_name}")
            
        except Exception as e:
            logger.error(f"Failed to process annotation file {annotation_file}: {str(e)}")
            logger.error(traceback.format_exc())
    
    # Replace empty cells with zero
    cog_counts_df = cog_counts_df.fillna(0)
    individual_cog_matrix = individual_cog_matrix.fillna(0)
    ko_counts_df = ko_counts_df.fillna(0)
    
    # Reorganize columns alphabetically
    cog_counts_df = cog_counts_df.reindex(sorted(cog_counts_df.columns), axis=1)
    individual_cog_matrix = individual_cog_matrix.reindex(sorted(individual_cog_matrix.columns), axis=1)
    ko_counts_df = ko_counts_df.reindex(sorted(ko_counts_df.columns), axis=1)
    
    # Write dataframes to Excel files in their respective folders
    cog_counts_output = os.path.join(folders["files"]["cog"], "cog_counts.xlsx")
    individual_cog_output = os.path.join(folders["files"]["cog"], "individual_cog_counts.xlsx")
    ko_counts_output = os.path.join(folders["files"]["ko"], "ko_counts.xlsx")
    
    cog_counts_df.to_excel(cog_counts_output, index_label="COG_term")
    individual_cog_matrix.to_excel(individual_cog_output, index_label="COG_term")
    ko_counts_df.to_excel(ko_counts_output, index_label="KO_term")
    
    # Also save as CSV for better compatibility
    cog_counts_df.to_csv(os.path.join(folders["files"]["cog"], "cog_counts.csv"), index_label="COG_term")
    individual_cog_matrix.to_csv(os.path.join(folders["files"]["cog"], "individual_cog_counts.csv"), index_label="COG_term")
    ko_counts_df.to_csv(os.path.join(folders["files"]["ko"], "ko_counts.csv"), index_label="KO_term")
    
    # Also save copies in the main analysis directory for backward compatibility
    cog_counts_df.to_excel(os.path.join(analysis_dir, "cog_counts.xlsx"), index_label="COG_term")
    individual_cog_matrix.to_excel(os.path.join(analysis_dir, "individual_cog_counts.xlsx"), index_label="COG_term")
    ko_counts_df.to_excel(os.path.join(analysis_dir, "ko_counts.xlsx"), index_label="KO_term")
    
    logger.info(f"COG counts saved to: {cog_counts_output}")
    logger.info(f"Individual COG counts saved to: {individual_cog_output}")
    logger.info(f"KO counts saved to: {ko_counts_output}")
    
    # Create a more readable version of COG counts with descriptions
    try:
        # Use individual COG matrix for descriptions
        if not individual_cog_matrix.empty:
            cog_with_desc = individual_cog_matrix.copy()
            cog_with_desc = cog_with_desc.reset_index()
            
            # Ensure the index is named properly
            if 'COG_term' in cog_with_desc.columns:
                # Add descriptions
                cog_with_desc['Description'] = cog_with_desc['COG_term'].apply(
                    lambda x: COG_CATEGORIES.get(x, 'Unknown') if isinstance(x, str) and len(x) == 1 else 'Unknown'
                )
                desc_file = os.path.join(folders["files"]["cog"], "cog_counts_with_descriptions.xlsx")
                cog_with_desc.to_excel(desc_file, index=False)
                # Also save to original location for backward compatibility
                cog_with_desc.to_excel(os.path.join(analysis_dir, "cog_counts_with_descriptions.xlsx"), index=False)
                logger.info(f"COG counts with descriptions saved to: {desc_file}")
            else:
                logger.warning("COG_term column not found in cog_with_desc dataframe")
        else:
            logger.warning("individual_cog_matrix is empty, cannot create descriptions")
    except Exception as e:
        logger.error(f"Failed to create COG counts with descriptions: {str(e)}")
        logger.error(traceback.format_exc())
        
    return cog_counts_df, individual_cog_matrix, ko_counts_df

def process_ko_annotations(ko_dir, analysis_dir, folders, ko_matrix=None):
    """
    Process KO annotation files to create abundance matrices.
    
    Parameters:
    -----------
    ko_dir : str
        Path to directory containing KO annotation files
    analysis_dir : str
        Path to analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    ko_matrix : DataFrame, optional
        Existing KO counts matrix
    """
    # First combine all KO files into a single matrix if not already done
    ko_matrix_file = os.path.join(folders["files"]["ko"], "ko_counts.xlsx")
    
    if ko_matrix is not None:
        # Save the provided matrix
        ko_matrix.to_excel(ko_matrix_file, index_label="KO_term")
        ko_matrix.to_csv(os.path.join(folders["files"]["ko"], "ko_counts.csv"), index_label="KO_term")
        logger.info(f"Using provided KO matrix saved to: {ko_matrix_file}")
    elif not os.path.exists(ko_matrix_file):
        logger.info("KO matrix file not found. Creating from individual KO files.")
        
        # Try to find it in the original location for backward compatibility
        original_ko_file = os.path.join(analysis_dir, "ko_counts.xlsx")
        if os.path.exists(original_ko_file):
            import shutil
            shutil.copy2(original_ko_file, ko_matrix_file)
            logger.info(f"Copied {original_ko_file} to {ko_matrix_file}")
        else:
            # Get all KO files
            ko_files = glob.glob(os.path.join(ko_dir, "*_KOs.tsv"))
            
            if not ko_files:
                logger.warning(f"No KO files found in {ko_dir}")
                return
            
            logger.info(f"Found {len(ko_files)} KO files")
            
            # Create a dataframe to store all KO counts
            ko_matrix = pd.DataFrame()
            
            for ko_file in ko_files:
                try:
                    sample_name = os.path.basename(ko_file).split('_KOs.tsv')[0]
                    
                    # Read the KO file
                    ko_df = pd.read_csv(ko_file, sep='\t', header=None, names=['gene', 'ko'])
                    
                    # Process KO IDs to standardize them
                    ko_df['ko'] = ko_df['ko'].apply(lambda x: x.replace('ko:', '') if isinstance(x, str) and x.startswith('ko:') else x)
                    
                    # Count KO terms
                    ko_counts = ko_df['ko'].value_counts()
                    
                    # Add to matrix
                    ko_matrix[sample_name] = ko_counts
                    
                    logger.info(f"Processed {sample_name}: found {len(ko_counts)} KO terms")
                    
                except Exception as e:
                    logger.error(f"Failed to process KO file {ko_file}: {str(e)}")
                    logger.error(traceback.format_exc())
            
            # Fill NaN with 0
            ko_matrix = ko_matrix.fillna(0)
            
            # Save the matrix
            ko_matrix.to_excel(ko_matrix_file, index_label="KO_term")
            ko_matrix.to_csv(os.path.join(folders["files"]["ko"], "ko_counts.csv"), index_label="KO_term")
            # Also save to original location for backward compatibility
            ko_matrix.to_excel(os.path.join(analysis_dir, "ko_counts.xlsx"), index_label="KO_term")
            logger.info(f"KO matrix saved to: {ko_matrix_file}")
    
    # Load the KO matrix if not provided
    if ko_matrix is None:
        try:
            ko_matrix = pd.read_excel(ko_matrix_file, index_col="KO_term")
        except Exception as e:
            logger.error(f"Failed to load KO matrix: {str(e)}")
            return
    
    # Process the KO matrix to split combined terms
    modified_ko_file = os.path.join(folders["files"]["ko"], "modified_ko_counts.xlsx")
    
    try:
        # Create a new dataframe for the modified data
        new_ko_df = pd.DataFrame()
        
        # Process each row to split combined terms
        for index, row in ko_matrix.iterrows():
            try:
                # Skip if not a string or doesn't contain a comma
                if not isinstance(index, str) or ',' not in index:
                    new_ko_df = pd.concat([new_ko_df, pd.DataFrame([row], index=[index])])
                    continue
                
                # Split combined KO terms
                ko_terms = [term.strip() for term in index.split(',')]
                
                # Divide counts evenly
                divided_counts = row / len(ko_terms)
                
                # Add each term with divided counts
                for term in ko_terms:
                    # Remove 'ko:' prefix if present
                    if isinstance(term, str) and term.startswith('ko:'):
                        term = term[3:]
                        
                    if term in new_ko_df.index:
                        # If term already exists, add the counts
                        new_ko_df.loc[term] += divided_counts
                    else:
                        # Otherwise create a new row
                        new_ko_df = pd.concat([new_ko_df, pd.DataFrame([divided_counts], index=[term])])
            except Exception as e:
                logger.warning(f"Error processing KO term {index}: {str(e)}")
                # Add the original row if there's an error
                new_ko_df = pd.concat([new_ko_df, pd.DataFrame([row], index=[index])])
        
        # Save the modified matrix
        new_ko_df.to_excel(modified_ko_file)
        new_ko_df.to_csv(os.path.join(folders["files"]["ko"], "modified_ko_counts.csv"))
        # Also save to original location for backward compatibility
        new_ko_df.to_excel(os.path.join(analysis_dir, "modified_ko_counts.xlsx"))
        logger.info(f"Modified KO matrix saved to: {modified_ko_file}")
        
        # Create a summarized KO counts file by summing all samples
        try:
            summarized_ko_file = os.path.join(folders["files"]["ko"], "summarized_ko_counts.xlsx")
            # Sum across all samples
            total_ko = new_ko_df.sum(axis=1).sort_values(ascending=False)
            total_ko_df = pd.DataFrame({'Total': total_ko})
            total_ko_df.to_excel(summarized_ko_file)
            # Also save to original location for backward compatibility
            total_ko_df.to_excel(os.path.join(analysis_dir, "summarized_ko_counts.xlsx"))
            logger.info(f"Summarized KO counts saved to: {summarized_ko_file}")
        except Exception as e:
            logger.error(f"Failed to create summarized KO counts: {str(e)}")
            logger.error(traceback.format_exc())
    except Exception as e:
        logger.error(f"Failed to process KO matrix: {str(e)}")
        logger.error(traceback.format_exc())

def process_cog_annotations(cog_dir, analysis_dir, folders, cog_matrix=None):
    """
    Process COG annotation files to create abundance matrices.
    
    Parameters:
    -----------
    cog_dir : str
        Path to directory containing COG annotation files
    analysis_dir : str
        Path to analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    cog_matrix : DataFrame, optional
        Existing individual COG matrix
    """
    # Check if individual COG counts exist, if not create them
    individual_cog_matrix_file = os.path.join(folders["files"]["cog"], "individual_cog_counts.xlsx")
    
    if cog_matrix is not None:
        # Save the provided matrix
        cog_matrix.to_excel(individual_cog_matrix_file, index_label="COG_term")
        cog_matrix.to_csv(os.path.join(folders["files"]["cog"], "individual_cog_counts.csv"), index_label="COG_term")
        logger.info(f"Using provided COG matrix saved to: {individual_cog_matrix_file}")
    elif not os.path.exists(individual_cog_matrix_file):
        logger.info("Individual COG matrix file not found. Creating from COG files.")
        
        # Try to find in original location for backward compatibility
        original_cog_file = os.path.join(analysis_dir, "individual_cog_counts.xlsx")
        if os.path.exists(original_cog_file):
            import shutil
            shutil.copy2(original_cog_file, individual_cog_matrix_file)
            logger.info(f"Copied {original_cog_file} to {individual_cog_matrix_file}")
        else:
            # Get all COG files
            cog_files = glob.glob(os.path.join(cog_dir, "*_COGs.tsv"))
            
            if not cog_files:
                logger.warning(f"No COG files found in {cog_dir}")
                return
            
            logger.info(f"Found {len(cog_files)} COG files")
            
            # Create a dataframe to store all individual COG letter counts
            individual_cog_matrix = pd.DataFrame()
            
            for cog_file in cog_files:
                try:
                    sample_name = os.path.basename(cog_file).split('_COGs.tsv')[0]
                    
                    # Read the COG file
                    cog_df = pd.read_csv(cog_file, sep='\t', header=None, names=['gene', 'cog'])
                    
                    # Process each COG entry to split multi-letter entries
                    individual_cogs = []
                    for cog in cog_df['cog']:
                        if isinstance(cog, str) and cog != '-':
                            for single_cog in cog:
                                if single_cog.isalpha():  # Only include alphabet characters
                                    individual_cogs.append(single_cog)
                    
                    # Count individual COG letters
                    cog_counts = pd.Series(Counter(individual_cogs))
                    
                    # Add to matrix
                    individual_cog_matrix[sample_name] = cog_counts
                    
                    logger.info(f"Processed {sample_name}: found {len(cog_counts)} COG categories")
                    
                except Exception as e:
                    logger.error(f"Failed to process COG file {cog_file}: {str(e)}")
                    logger.error(traceback.format_exc())
            
            # Fill NaN with 0
            individual_cog_matrix = individual_cog_matrix.fillna(0)
            
            # Save the matrix
            individual_cog_matrix.to_excel(individual_cog_matrix_file, index_label="COG_term")
            individual_cog_matrix.to_csv(os.path.join(folders["files"]["cog"], "individual_cog_counts.csv"), index_label="COG_term")
            # Also save to original location for backward compatibility
            individual_cog_matrix.to_excel(os.path.join(analysis_dir, "individual_cog_counts.xlsx"), index_label="COG_term")
            logger.info(f"Individual COG matrix saved to: {individual_cog_matrix_file}")
    
    # Load the COG matrix if not provided
    if cog_matrix is None:
        try:
            cog_matrix = pd.read_excel(individual_cog_matrix_file, index_col="COG_term")
        except Exception as e:
            logger.error(f"Failed to load COG matrix: {str(e)}")
            return
    
    # Also process the original multi-letter COG categories
    try:
        # Get all COG files
        cog_files = glob.glob(os.path.join(cog_dir, "*_COGs.tsv"))
        
        # Create a dataframe to store multi-letter COG counts
        multi_cog_matrix = pd.DataFrame()
        
        for cog_file in cog_files:
            try:
                sample_name = os.path.basename(cog_file).split('_COGs.tsv')[0]
                
                # Read the COG file
                cog_df = pd.read_csv(cog_file, sep='\t', header=None, names=['gene', 'cog'])
                
                # Only keep multi-letter COGs
                multi_cogs = cog_df[cog_df['cog'].apply(lambda x: isinstance(x, str) and len(x) > 1)]
                
                if not multi_cogs.empty:
                    # Count multi-letter COG terms
                    cog_counts = multi_cogs['cog'].value_counts()
                    
                    # Add to matrix
                    multi_cog_matrix[sample_name] = cog_counts
                    
                    logger.info(f"Found {len(cog_counts)} multi-letter COG categories in {sample_name}")
                
            except Exception as e:
                logger.error(f"Failed to process multi-letter COGs in {cog_file}: {str(e)}")
                logger.error(traceback.format_exc())
        
        if not multi_cog_matrix.empty:
            # Fill NaN with 0
            multi_cog_matrix = multi_cog_matrix.fillna(0)
            
            # Save the matrix
            multi_cog_matrix_file = os.path.join(folders["files"]["cog"], "multi_letter_cog_counts.xlsx")
            multi_cog_matrix.to_excel(multi_cog_matrix_file, index_label="COG_term")
            multi_cog_matrix.to_csv(os.path.join(folders["files"]["cog"], "multi_letter_cog_counts.csv"), index_label="COG_term")
            # Also save to original location for backward compatibility
            multi_cog_matrix.to_excel(os.path.join(analysis_dir, "multi_letter_cog_counts.xlsx"), index_label="COG_term")
            logger.info(f"Multi-letter COG matrix saved to: {multi_cog_matrix_file}")
    except Exception as e:
        logger.error(f"Failed to process multi-letter COG matrix: {str(e)}")
        logger.error(traceback.format_exc())
        
    # Create a summarized COG counts file
    try:
        # Make sure it has a COG_term column (or is already indexed by COG_term)
        if isinstance(cog_matrix.index, pd.Index) and cog_matrix.index.name == "COG_term":
            individual_cog_matrix = cog_matrix.copy()
        else:
            # Try to set COG_term as index if it's a column
            if 'COG_term' in cog_matrix.columns:
                individual_cog_matrix = cog_matrix.set_index('COG_term')
            else:
                logger.warning("Cannot find COG_term in the COG matrix")
                return
        
        # Sum across all samples
        total_cog = individual_cog_matrix.sum(axis=1).sort_values(ascending=False)
        total_cog_df = pd.DataFrame({'Total': total_cog})
        
        # Reset index for adding descriptions
        total_cog_df = total_cog_df.reset_index()
        
        # Add descriptions
        total_cog_df['Description'] = total_cog_df['COG_term'].apply(
            lambda x: COG_CATEGORIES.get(x, 'Unknown') if isinstance(x, str) else 'Unknown'
        )
        
        # Save the summarized matrix
        summarized_cog_file = os.path.join(folders["files"]["cog"], "summarized_cog_counts.xlsx")
        total_cog_df.to_excel(summarized_cog_file, index=False)
        total_cog_df.to_csv(os.path.join(folders["files"]["cog"], "summarized_cog_counts.csv"), index=False)
        # Also save to original location for backward compatibility
        total_cog_df.to_excel(os.path.join(analysis_dir, "summarized_cog_counts.xlsx"), index=False)
        logger.info(f"Summarized COG counts saved to: {summarized_cog_file}")
    except Exception as e:
        logger.error(f"Failed to create summarized COG counts: {str(e)}")
        logger.error(traceback.format_exc())

def map_ko_to_kegg_hierarchy(analysis_dir, folders):
    """
    Map KO terms to KEGG hierarchy levels.
    
    Parameters:
    -----------
    analysis_dir : str
        Path to the analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    
    Returns:
    --------
    DataFrame
        Merged KO abundance with KEGG hierarchy levels
    """
    ko_file = os.path.join(folders["files"]["ko"], "modified_ko_counts.xlsx")
    kegg_db_file = os.path.join(folders["files"]["kegg"], "ko.txt")
    
    # Check if files exist at new paths, otherwise try original paths
    if not os.path.exists(ko_file):
        original_ko_file = os.path.join(analysis_dir, "modified_ko_counts.xlsx")
        if os.path.exists(original_ko_file):
            import shutil
            shutil.copy2(original_ko_file, ko_file)
            logger.info(f"Copied {original_ko_file} to {ko_file}")
        else:
            logger.warning(f"Modified KO file not found: {ko_file}")
            return None
    
    if not os.path.exists(kegg_db_file):
        original_kegg_file = os.path.join(analysis_dir, "ko.txt")
        if os.path.exists(original_kegg_file):
            import shutil
            shutil.copy2(original_kegg_file, kegg_db_file)
            logger.info(f"Copied {original_kegg_file} to {kegg_db_file}")
        else:
            logger.warning(f"KEGG database file not found: {kegg_db_file}")
            return None
    
    # Extract KEGG hierarchy levels from ko.txt
    logger.info(f"Extracting KEGG hierarchy levels from {kegg_db_file}")
    
    # Initialize lists to store extracted data
    level_1, level_2, level_3, ko_ids = [], [], [], []
    
    # Variables to store current hierarchy levels
    current_level_1 = None
    current_level_2 = None
    current_level_3 = None
    
    try:
        # Read the KO hierarchy file and extract the levels
        with open(kegg_db_file, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("A"):  # Level 1
                    current_level_1 = line[1:].strip()
                elif line.startswith("B"):  # Level 2
                    current_level_2 = line[2:].strip()
                elif line.startswith("C"):  # Level 3
                    current_level_3 = re.sub(r"\[.*?\]", "", line[4:]).strip()  # Remove [PATH:koXXXXX]
                elif line.startswith("D"):  # KO ID
                    match = re.match(r"D\s+(K\d+)", line)
                    if match:
                        ko_id = match.group(1)
                        level_1.append(current_level_1)
                        level_2.append(current_level_2)
                        level_3.append(current_level_3)
                        ko_ids.append(ko_id)
        
        logger.info(f"Extracted {len(ko_ids)} KO IDs with hierarchy information")
        
        # Create a DataFrame with extracted KO hierarchy
        ko_hierarchy_df = pd.DataFrame({
            "KO_ID": ko_ids,
            "Level_1": level_1,
            "Level_2": level_2,
            "Level_3": level_3
        })
        
        # Save the KO hierarchy mapping
        hierarchy_file = os.path.join(folders["files"]["kegg"], "ko_hierarchy_mapping.tsv")
        ko_hierarchy_df.to_csv(hierarchy_file, sep="\t", index=False)
        # Also save to original location for backward compatibility
        ko_hierarchy_df.to_csv(os.path.join(analysis_dir, "ko_hierarchy_mapping.tsv"), sep="\t", index=False)
        logger.info(f"KO hierarchy mapping saved to: {hierarchy_file}")
        
        # Load KO abundance file
        try:
            ko_abundance_df = pd.read_excel(ko_file)
            
            # Ensure KO_term column exists
            if ko_abundance_df.index.name == "KO_term":
                ko_abundance_df = ko_abundance_df.reset_index()
            
            # If KO_term is not in columns, try to add it
            if "KO_term" not in ko_abundance_df.columns:
                # If the dataframe has an index column, rename it
                if ko_abundance_df.index.name is not None:
                    ko_abundance_df = ko_abundance_df.reset_index()
                    ko_abundance_df = ko_abundance_df.rename(columns={ko_abundance_df.columns[0]: "KO_term"})
                # Otherwise, use the first column as KO_term
                elif len(ko_abundance_df.columns) > 0:
                    first_col = ko_abundance_df.columns[0]
                    ko_abundance_df = ko_abundance_df.rename(columns={first_col: "KO_term"})
                else:
                    logger.error("Cannot find or create KO_term column in KO abundance file")
                    return None
            
            # Clean up KO IDs and convert to strings for matching
            ko_abundance_df["KO_ID"] = ko_abundance_df["KO_term"].apply(
                lambda x: str(x).replace("ko:", "") if isinstance(x, str) and x.startswith("ko:") else str(x)
            )
            
            # Merge KO abundance with KEGG hierarchy
            merged_abundance_df = ko_abundance_df.merge(ko_hierarchy_df, on="KO_ID", how="left")
            
            # Save the final merged file
            merged_file = os.path.join(folders["files"]["kegg"], "ko_abundance_with_levels.tsv")
            merged_abundance_df.to_csv(merged_file, sep="\t", index=False)
            # Also save as Excel for easier viewing
            merged_abundance_df.to_excel(os.path.join(folders["files"]["kegg"], "ko_abundance_with_levels.xlsx"), index=False)
            # Also save to original location for backward compatibility
            merged_abundance_df.to_csv(os.path.join(analysis_dir, "ko_abundance_with_levels.tsv"), sep="\t", index=False)
            merged_abundance_df.to_excel(os.path.join(analysis_dir, "ko_abundance_with_levels.xlsx"), index=False)
            logger.info(f"KO abundance with KEGG hierarchy levels saved to: {merged_file}")
            
            # Log the match rate
            total_kos = len(ko_abundance_df)
            matched_kos = merged_abundance_df["Level_1"].notnull().sum()
            match_percentage = (matched_kos/total_kos*100) if total_kos > 0 else 0
            logger.info(f"Matched {matched_kos} out of {total_kos} KO terms ({match_percentage:.2f}%)")
            
            # Create aggregated files by KEGG levels
            try:
                # Remove rows without hierarchy information
                merged_for_agg = merged_abundance_df.dropna(subset=["Level_1"])
                
                if len(merged_for_agg) == 0:
                    logger.warning("No KO terms matched with KEGG hierarchy")
                    return merged_abundance_df
                    
                # Get numeric columns (sample abundance data)
                numeric_cols = merged_for_agg.select_dtypes(include=[np.number]).columns
                numeric_cols = [col for col in numeric_cols if col != "KO_ID"]
                
                if len(numeric_cols) == 0:
                    logger.warning("No numeric abundance data found in merged KO file")
                    return merged_abundance_df
                    
                # Aggregate by Level 1
                level1_df = merged_for_agg.groupby("Level_1")[numeric_cols].sum().reset_index()
                level1_file = os.path.join(folders["files"]["kegg"], "ko_abundance_by_level1.xlsx")
                level1_df.to_excel(level1_file, index=False)
                # Also save to original location for backward compatibility
                level1_df.to_excel(os.path.join(analysis_dir, "ko_abundance_by_level1.xlsx"), index=False)
                logger.info(f"KO abundance aggregated by Level 1 saved to: {level1_file}")
                
                # Aggregate by Level 2
                level2_df = merged_for_agg.groupby("Level_2")[numeric_cols].sum().reset_index()
                level2_file = os.path.join(folders["files"]["kegg"], "ko_abundance_by_level2.xlsx")
                level2_df.to_excel(level2_file, index=False)
                # Also save to original location for backward compatibility
                level2_df.to_excel(os.path.join(analysis_dir, "ko_abundance_by_level2.xlsx"), index=False)
                logger.info(f"KO abundance aggregated by Level 2 saved to: {level2_file}")
                
                # Aggregate by Level 3
                level3_df = merged_for_agg.groupby("Level_3")[numeric_cols].sum().reset_index()
                level3_file = os.path.join(folders["files"]["kegg"], "ko_abundance_by_level3.xlsx")
                level3_df.to_excel(level3_file, index=False)
                # Also save to original location for backward compatibility
                level3_df.to_excel(os.path.join(analysis_dir, "ko_abundance_by_level3.xlsx"), index=False)
                logger.info(f"KO abundance aggregated by Level 3 saved to: {level3_file}")
                
                # Return the merged dataframe for further processing
                return merged_abundance_df
                
            except Exception as e:
                logger.error(f"Failed to aggregate KO abundance by KEGG levels: {str(e)}")
                logger.error(traceback.format_exc())
                return merged_abundance_df
                
        except Exception as e:
            logger.error(f"Failed to load or process KO abundance file: {str(e)}")
            logger.error(traceback.format_exc())
            return None
        
    except Exception as e:
        logger.error(f"Failed to map KO terms to KEGG hierarchy: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def plot_cog_by_main_categories(analysis_dir, folders):
    """
    Generate additional COG plots organized by main functional categories
    with full names instead of just symbols.
    
    Parameters:
    -----------
    analysis_dir : str
        Path to the analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    """
    individual_cog_file = os.path.join(folders["files"]["cog"], "individual_cog_counts.xlsx")
    
    if not os.path.exists(individual_cog_file):
        # If not found in the new location, try the original location
        original_file = os.path.join(analysis_dir, "individual_cog_counts.xlsx")
        if os.path.exists(original_file):
            import shutil
            shutil.copy2(original_file, individual_cog_file)
            logger.info(f"Copied {original_file} to {individual_cog_file}")
        else:
            logger.warning(f"Individual COG counts file not found: {individual_cog_file}")
            return
    
    try:
        # Read the COG data
        cog_df = pd.read_excel(individual_cog_file)
        
        if "COG_term" not in cog_df.columns:
            logger.warning("COG_term column not found in COG file")
            return
        
        # Set COG_term as index
        cog_df = cog_df.set_index("COG_term")
        
        # Remove any non-standard COG categories
        valid_cogs = [cog for cog in COG_CATEGORIES.keys()]
        cog_df = cog_df.loc[cog_df.index.isin(valid_cogs)]
        
        # Add full names and main categories
        cog_data = []
        for cog, row in cog_df.iterrows():
            if not isinstance(cog, str):
                continue
                
            # Find which main category this COG belongs to
            main_category = "OTHER"
            for category, cogs in COG_MAIN_CATEGORIES.items():
                if cog in cogs:
                    main_category = category
                    break
                    
            # Get the full description
            description = COG_CATEGORIES.get(cog, "Unknown")
            
            # Add data for each sample
            for sample in cog_df.columns:
                cog_data.append({
                    "COG": cog,
                    "Main_Category": main_category,
                    "Description": description,
                    "Full_Name": f"{cog} - {description}",
                    "Sample": sample,
                    "Count": row[sample]
                })
        
        # Create a dataframe from the collected data
        full_cog_df = pd.DataFrame(cog_data)
        
        if full_cog_df.empty:
            logger.warning("No valid COG data found for enhanced plotting")
            return
        
        # 1. Plot grouped by main category
        plt.figure(figsize=(15, 10))
        
        # Sum counts by main category and sample
        main_cat_counts = full_cog_df.groupby(['Main_Category', 'Sample'])['Count'].sum().unstack()
        
        # Create a stacked bar chart
        ax = main_cat_counts.plot(kind='bar', stacked=True, figsize=(15, 10))
        
        plt.title('COG Main Functional Categories Distribution')
        plt.ylabel('Count')
        plt.xlabel('Functional Category')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        # Save the plot
        main_cat_plot = os.path.join(folders["plots"]["cog"], "cog_main_categories.pdf")
        plt.savefig(main_cat_plot)
        # Also save as PNG for web viewing
        main_cat_plot_png = os.path.join(folders["plots"]["cog"], "cog_main_categories.png")
        plt.savefig(main_cat_plot_png, dpi=300)
        plt.close()
        logger.info(f"COG main categories plot saved to: {main_cat_plot}")
        
        # 2. Create detailed barplot with all COGs and their descriptions
        # Group by COG and get total counts across samples
        cog_totals = full_cog_df.groupby(['COG', 'Full_Name', 'Main_Category'])['Count'].sum().reset_index()
        cog_totals = cog_totals.sort_values('Count', ascending=False)
        
        # Assign colors based on main category
        cog_colors = [COG_COLORS.get(cog, "#999999") for cog in cog_totals['COG']]
        
        plt.figure(figsize=(20, 12))
        bars = plt.barh(cog_totals['Full_Name'], cog_totals['Count'], color=cog_colors)
        
        plt.title('COG Categories with Descriptions (Total Counts)')
        plt.xlabel('Count')
        plt.ylabel('COG Category')
        plt.tight_layout()
        
        # Save the plot
        detailed_cog_plot = os.path.join(folders["plots"]["cog"], "cog_detailed_barplot.pdf")
        plt.savefig(detailed_cog_plot)
        # Also save as PNG for web viewing
        detailed_cog_plot_png = os.path.join(folders["plots"]["cog"], "cog_detailed_barplot.png")
        plt.savefig(detailed_cog_plot_png, dpi=300)
        plt.close()
        logger.info(f"Detailed COG plot saved to: {detailed_cog_plot}")
        
        # 3. Create a heatmap of all COGs with descriptions
        pivot_df = full_cog_df.pivot_table(index='Full_Name', columns='Sample', values='Count', fill_value=0)
        
        # Sort by main category and then by COG
        cog_order = full_cog_df.sort_values(['Main_Category', 'COG']).drop_duplicates('Full_Name')['Full_Name']
        pivot_df = pivot_df.reindex(cog_order)
        
        plt.figure(figsize=(15, 20))
        plt.title('COG Categories Heatmap with Descriptions')
        
        # Create the heatmap
        heatmap = plt.pcolor(pivot_df.values, cmap='viridis')
        plt.colorbar(heatmap)
        
        # Set axis labels
        plt.yticks(np.arange(0.5, len(pivot_df.index)), pivot_df.index, fontsize=8)
        plt.xticks(np.arange(0.5, len(pivot_df.columns)), pivot_df.columns, rotation=90)
        
        plt.xlabel('Sample')
        plt.ylabel('COG Category')
        
        plt.tight_layout()
        
        # Save the heatmap
        cog_heatmap = os.path.join(folders["plots"]["cog"], "cog_detailed_heatmap.pdf")
        plt.savefig(cog_heatmap)
        # Also save as PNG for web viewing
        cog_heatmap_png = os.path.join(folders["plots"]["cog"], "cog_detailed_heatmap.png")
        plt.savefig(cog_heatmap_png, dpi=300)
        plt.close()
        logger.info(f"Detailed COG heatmap saved to: {cog_heatmap}")
        
        # 4. Create a summary table of COGs by main category
        summary_df = full_cog_df.pivot_table(
            index=['Main_Category', 'COG', 'Description'], 
            columns='Sample', 
            values='Count', 
            aggfunc='sum',
            fill_value=0
        ).reset_index()
        
        summary_file = os.path.join(folders["files"]["cog"], "cog_summary_by_category.xlsx")
        summary_df.to_excel(summary_file, index=False)
        logger.info(f"COG summary by category saved to: {summary_file}")
        
        # 5. Create a stacked area chart for COGs by main category
        plt.figure(figsize=(15, 10))
        
        # Prepare data - transpose to get samples on x-axis
        area_data = main_cat_counts.T
        
        # Create stacked area chart
        area_data.plot.area(figsize=(15, 10), alpha=0.6, stacked=True)
        
        plt.title('COG Main Categories Across Samples')
        plt.xlabel('Sample')
        plt.ylabel('Count')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save the plot
        area_plot = os.path.join(folders["plots"]["cog"], "cog_main_categories_area.pdf")
        plt.savefig(area_plot)
        # Also save as PNG for web viewing
        area_plot_png = os.path.join(folders["plots"]["cog"], "cog_main_categories_area.png")
        plt.savefig(area_plot_png, dpi=300)
        plt.close()
        logger.info(f"COG main categories area plot saved to: {area_plot}")
        
        # 6. Create a bubble plot showing distribution of COGs
        # Get the top COGs
        top_cogs = cog_totals.head(20)
        
        plt.figure(figsize=(16, 10))
        
        # Create a categorical colormap based on main categories
        unique_categories = top_cogs['Main_Category'].unique()
        category_colors = plt.cm.tab10(np.linspace(0, 1, len(unique_categories)))
        color_map = dict(zip(unique_categories, category_colors))
        
        # Create bubble chart
        for i, (_, row) in enumerate(top_cogs.iterrows()):
            plt.scatter(
                i, row['Count'], 
                s=row['Count']*5, 
                color=color_map[row['Main_Category']], 
                alpha=0.7, 
                edgecolor='black'
            )
            plt.text(i, row['Count'], row['COG'], ha='center', va='center')
        
        # Add legend for main categories
        for category, color in zip(unique_categories, category_colors):
            plt.scatter([], [], color=color, label=category, s=100)
        
        plt.legend(title="Main Categories", loc='upper right')
        
        plt.title('Top 20 COG Categories by Count')
        plt.xlabel('COG Rank')
        plt.ylabel('Count')
        plt.xticks(range(len(top_cogs)), top_cogs['COG'])
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        
        # Save the plot
        bubble_plot = os.path.join(folders["plots"]["cog"], "cog_bubble_plot.pdf")
        plt.savefig(bubble_plot)
        # Also save as PNG for web viewing
        bubble_plot_png = os.path.join(folders["plots"]["cog"], "cog_bubble_plot.png")
        plt.savefig(bubble_plot_png, dpi=300)
        plt.close()
        logger.info(f"COG bubble plot saved to: {bubble_plot}")
        
    except Exception as e:
        logger.error(f"Failed to generate enhanced COG plots: {str(e)}")
        logger.error(traceback.format_exc())

def generate_plots(analysis_dir, folders):
    """
    Generate basic plots for functional annotation analysis.
    
    Parameters:
    -----------
    analysis_dir : str
        Path to the analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    """
    # Generate plots for COG categories
    cog_files = [
        os.path.join(folders["files"]["cog"], "individual_cog_counts.xlsx"),
        os.path.join(folders["files"]["cog"], "summarized_cog_counts.xlsx"),
        os.path.join(folders["files"]["cog"], "cog_counts_with_descriptions.xlsx")
    ]
    
    # Check if files exist, try original locations if not found
    for i, file_path in enumerate(cog_files):
        if not os.path.exists(file_path):
            original_file = os.path.join(analysis_dir, os.path.basename(file_path))
            if os.path.exists(original_file):
                import shutil
                shutil.copy2(original_file, file_path)
                logger.info(f"Copied {original_file} to {file_path}")
                # Update the list with the existing file
                cog_files[i] = file_path
            else:
                # Replace with original path in case it exists but wasn't detected
                cog_files[i] = original_file
    
    for cog_file in cog_files:
        if os.path.exists(cog_file):
            try:
                logger.info(f"Generating plot for {cog_file}")
                
                # Read the COG data
                cog_df = pd.read_excel(cog_file)
                
                # If it's a summarized file with descriptions, handle differently
                if "Description" in cog_df.columns:
                    # Set COG_term as index
                    if "COG_term" in cog_df.columns:
                        cog_df = cog_df.set_index("COG_term")
                    
                    # Remove Description column for plotting
                    plot_df = cog_df.drop(columns=["Description"])
                    
                    # Sort by total counts descending
                    if "Total" in plot_df.columns:
                        plot_df = plot_df.sort_values(by="Total", ascending=False)
                    
                    # Create a bar plot
                    plt.figure(figsize=(15, 10))
                    ax = plot_df.plot(kind='bar', figsize=(15, 10))
                    
                    # Add COG descriptions to y-axis labels
                    if not plot_df.empty:
                        plt.yticks(range(len(plot_df.index)), 
                                  [f"{idx} - {cog_df.loc[idx, 'Description']}" for idx in plot_df.index])
                    
                    # Add title and labels
                    plt.title('COG Categories Distribution')
                    plt.ylabel('COG Category')
                    plt.xlabel('Count')
                    
                    # Adjust layout
                    plt.tight_layout()
                    
                    # Save the plot
                    plot_file = os.path.join(
                        folders["plots"]["cog"], 
                        os.path.basename(os.path.splitext(cog_file)[0]) + "_barplot.png"
                    )
                    plt.savefig(plot_file, dpi=300)
                    plt.close()
                    
                    logger.info(f"COG bar plot saved to: {plot_file}")
                    
                else:
                    # Standard case - stacked bar plot for samples
                    
                    # If the index column is named, set it as index
                    if "COG_term" in cog_df.columns:
                        cog_df = cog_df.set_index("COG_term")
                    
                    # Generate a stacked bar plot
                    plt.figure(figsize=(15, 10))
                    ax = cog_df.T.plot(kind='bar', stacked=True, figsize=(15, 10))
                    
                    # Add title and labels
                    plt.title('COG Categories Distribution')
                    plt.ylabel('Count')
                    plt.xlabel('Sample')
                    
                    # Adjust legend
                    plt.legend(title='COG Categories', bbox_to_anchor=(1.05, 1), loc='upper left')
                    
                    # Save the plot
                    plot_file = os.path.join(
                        folders["plots"]["cog"], 
                        os.path.basename(os.path.splitext(cog_file)[0]) + "_stacked_barplot.png"
                    )
                    plt.tight_layout()
                    plt.savefig(plot_file, dpi=300)
                    plt.close()
                    
                    logger.info(f"COG stacked bar plot saved to: {plot_file}")
                    
                    # Also generate a heatmap for better visualization
                    plt.figure(figsize=(15, 15))
                    plt.title('COG Categories Heatmap')
                    
                    # Convert to numpy array for pcolor
                    heatmap_data = cog_df.T.values
                    heatmap = plt.pcolor(heatmap_data, cmap='viridis')
                    plt.colorbar(heatmap)
                    
                    # Set x and y labels
                    plt.xticks(np.arange(0.5, len(cog_df.index)), cog_df.index, rotation=90)
                    plt.yticks(np.arange(0.5, len(cog_df.columns)), cog_df.columns)
                    
                    plt.xlabel('COG Categories')
                    plt.ylabel('Sample')
                    
                    # Save the heatmap
                    heatmap_file = os.path.join(
                        folders["plots"]["cog"], 
                        os.path.basename(os.path.splitext(cog_file)[0]) + "_heatmap.png"
                    )
                    plt.tight_layout()
                    plt.savefig(heatmap_file, dpi=300)
                    plt.close()
                    
                    logger.info(f"COG heatmap saved to: {heatmap_file}")
                
            except Exception as e:
                logger.error(f"Failed to generate COG plot for {cog_file}: {str(e)}")
                logger.error(traceback.format_exc())
    
    # Generate plots for KO terms aggregated by KEGG levels
    kegg_files = [
        os.path.join(folders["files"]["kegg"], "ko_abundance_by_level1.xlsx"),
        os.path.join(folders["files"]["kegg"], "ko_abundance_by_level2.xlsx"),
        os.path.join(folders["files"]["kegg"], "ko_abundance_by_level3.xlsx")
    ]
    
    # Check if files exist, try original locations if not found
    for i, file_path in enumerate(kegg_files):
        if not os.path.exists(file_path):
            original_file = os.path.join(analysis_dir, os.path.basename(file_path))
            if os.path.exists(original_file):
                import shutil
                shutil.copy2(original_file, file_path)
                logger.info(f"Copied {original_file} to {file_path}")
                # Update the list with the existing file
                kegg_files[i] = file_path
            else:
                # Replace with original path in case it exists but wasn't detected
                kegg_files[i] = original_file
    
    for kegg_file in kegg_files:
        if os.path.exists(kegg_file):
            try:
                logger.info(f"Generating plot for {kegg_file}")
                
                # Read the KEGG data
                level = "Level_1" if "level1" in kegg_file else "Level_2" if "level2" in kegg_file else "Level_3"
                kegg_df = pd.read_excel(kegg_file)
                
                # Set the level column as index
                if level in kegg_df.columns:
                    kegg_df = kegg_df.set_index(level)
                
                # Remove non-numeric columns
                numeric_cols = kegg_df.select_dtypes(include=[np.number]).columns
                kegg_df = kegg_df[numeric_cols]
                
                # Generate a stacked bar plot if there's data
                if not kegg_df.empty:
                    # For Level 3, we need to handle differently due to potentially large number of categories
                    if level == "Level_3" and len(kegg_df) > 30:
                        # Sum across all samples
                        total_counts = kegg_df.sum(axis=1).sort_values(ascending=False)
                        
                        # Take top 30 categories
                        top_categories = total_counts.head(30).index
                        top_df = kegg_df.loc[top_categories]
                        
                        # Plot top categories
                        plt.figure(figsize=(15, 10))
                        ax = top_df.T.plot(kind='bar', stacked=True, figsize=(15, 10))
                        
                        plt.title(f'Top 30 KEGG {level} Categories')
                        plt.ylabel('Count')
                        plt.xlabel('Sample')
                        
                        # Adjust legend
                        plt.legend(title=f'KEGG {level}', bbox_to_anchor=(1.05, 1), loc='upper left')
                        
                        # Save the plot
                        plot_file = os.path.join(
                            folders["plots"]["kegg"], 
                            os.path.basename(os.path.splitext(kegg_file)[0]) + "_top30_plot.png"
                        )
                        plt.tight_layout()
                        plt.savefig(plot_file, dpi=300)
                        plt.close()
                        
                        logger.info(f"Top 30 KEGG {level} plot saved to: {plot_file}")
                        
                        # Also create a summary bar plot
                        plt.figure(figsize=(15, 10))
                        total_counts.head(30).plot(kind='bar', figsize=(15, 10))
                        
                        plt.title(f'Top 30 KEGG {level} Categories (Total Counts)')
                        plt.ylabel('Count')
                        plt.xlabel('KEGG Category')
                        plt.xticks(rotation=90)
                        
                        # Save the plot
                        summary_plot_file = os.path.join(
                            folders["plots"]["kegg"], 
                            os.path.basename(os.path.splitext(kegg_file)[0]) + "_top30_summary.png"
                        )
                        plt.tight_layout()
                        plt.savefig(summary_plot_file, dpi=300)
                        plt.close()
                        
                        logger.info(f"Top 30 KEGG {level} summary plot saved to: {summary_plot_file}")
                        
                    else:
                        # Standard plot for Level 1 and 2 (or small Level 3)
                        plt.figure(figsize=(15, 10))
                        ax = kegg_df.T.plot(kind='bar', stacked=True, figsize=(15, 10))
                        
                        # Add title and labels
                        plt.title(f'KEGG {level} Distribution')
                        plt.ylabel('Count')
                        plt.xlabel('Sample')
                        
                        # Adjust legend
                        plt.legend(title=f'KEGG {level}', bbox_to_anchor=(1.05, 1), loc='upper left')
                        
                        # Save the plot
                        plot_file = os.path.join(
                            folders["plots"]["kegg"], 
                            os.path.basename(os.path.splitext(kegg_file)[0]) + "_plot.png"
                        )
                        plt.tight_layout()
                        plt.savefig(plot_file, dpi=300)
                        plt.close()
                        
                        logger.info(f"KEGG {level} plot saved to: {plot_file}")
                        
                        # For level 1 and 2, also generate a heatmap
                        if level in ["Level_1", "Level_2"]:
                            plt.figure(figsize=(15, 15))
                            plt.title(f'KEGG {level} Heatmap')
                            
                            # Convert to numpy array for pcolor
                            heatmap_data = kegg_df.T.values
                            heatmap = plt.pcolor(heatmap_data, cmap='viridis')
                            plt.colorbar(heatmap)
                            
                            # Set x and y labels
                            plt.xticks(np.arange(0.5, len(kegg_df.index)), kegg_df.index, rotation=90)
                            plt.yticks(np.arange(0.5, len(kegg_df.columns)), kegg_df.columns)
                            
                            plt.xlabel(f'KEGG {level}')
                            plt.ylabel('Sample')
                            
                            # Save the heatmap
                            heatmap_file = os.path.join(
                                folders["plots"]["kegg"], 
                                os.path.basename(os.path.splitext(kegg_file)[0]) + "_heatmap.png"
                            )
                            plt.tight_layout()
                            plt.savefig(heatmap_file, dpi=300)
                            plt.close()
                            
                            logger.info(f"KEGG {level} heatmap saved to: {heatmap_file}")
                            
                        # Also create a summary bar plot
                        plt.figure(figsize=(15, 10))
                        kegg_df.sum(axis=1).sort_values(ascending=False).plot(kind='bar', figsize=(15, 10))
                        
                        plt.title(f'KEGG {level} Categories (Total Counts)')
                        plt.ylabel('Count')
                        plt.xlabel('KEGG Category')
                        plt.xticks(rotation=90)
                        
                        # Save the plot
                        summary_plot_file = os.path.join(
                            folders["plots"]["kegg"], 
                            os.path.basename(os.path.splitext(kegg_file)[0]) + "_summary.png"
                        )
                        plt.tight_layout()
                        plt.savefig(summary_plot_file, dpi=300)
                        plt.close()
                        
                        logger.info(f"KEGG {level} summary plot saved to: {summary_plot_file}")
                else:
                    logger.warning(f"No numeric data in {kegg_file} to plot")
                
            except Exception as e:
                logger.error(f"Failed to generate KEGG plot for {kegg_file}: {str(e)}")
                logger.error(traceback.format_exc())
    
    # Generate plots for full KO dataset
    ko_files = [
        os.path.join(folders["files"]["ko"], "modified_ko_counts.xlsx"),
        os.path.join(folders["files"]["ko"], "ko_counts.xlsx")
    ]
    
    # Check if files exist, try original locations if not found
    for i, file_path in enumerate(ko_files):
        if not os.path.exists(file_path):
            original_file = os.path.join(analysis_dir, os.path.basename(file_path))
            if os.path.exists(original_file):
                import shutil
                shutil.copy2(original_file, file_path)
                logger.info(f"Copied {original_file} to {file_path}")
                # Update the list with the existing file
                ko_files[i] = file_path
            else:
                # Replace with original path in case it exists but wasn't detected
                ko_files[i] = original_file
    
    for ko_file in ko_files:
        if os.path.exists(ko_file):
            try:
                logger.info(f"Generating summary plot for {ko_file}")
                
                # Read the KO data
                ko_df = pd.read_excel(ko_file)
                
                # Debug information
                logger.info(f"DataFrame shape: {ko_df.shape}")
                logger.info(f"DataFrame columns: {ko_df.columns.tolist()}")
                logger.info(f"DataFrame dtypes: {ko_df.dtypes}")
            
                # Set KO_term as index if present
                if "KO_term" in ko_df.columns:
                    ko_df = ko_df.set_index("KO_term")
                    logger.info("Set KO_term as index")
            
                # First, ensure all data columns are numeric
                # Get all columns except the index
                data_columns = ko_df.columns.tolist()
            
                # Convert non-numeric columns to numeric where possible
                for col in data_columns:
                     if not pd.api.types.is_numeric_dtype(ko_df[col]):
                         logger.warning(f"Column {col} is not numeric. Attempting to convert.")
                         ko_df[col] = pd.to_numeric(ko_df[col], errors='coerce')
            
                # Filter out any rows with all NaN values after conversion
                ko_df = ko_df.dropna(how='all')
            
                # Check if dataframe is empty after cleaning
                if ko_df.empty:
                    logger.warning(f"Dataframe is empty after cleaning: {ko_file}")
                    continue
            
                # Create a manual sum of rows to avoid type errors
                # Initialize with zeros
                total_counts = pd.Series(0.0, index=ko_df.index)
            
                # Sum each column individually
                for col in data_columns:
                    # Use only numeric values for the sum
                    total_counts = total_counts + ko_df[col].fillna(0)
            
                # Sort the results
                total_counts = total_counts.sort_values(ascending=False)
            
                # Check if there's data to plot
                if total_counts.empty:
                    logger.warning(f"No data to plot in {ko_file}")
                    continue
            
                # Take top 50 KO terms for plotting
                top_ko = total_counts.head(50)
            
                plt.figure(figsize=(15, 10))
                top_ko.plot(kind='bar', figsize=(15, 10))
            
                plt.title('Top 50 KO Terms by Total Count')
                plt.ylabel('Count')
                plt.xlabel('KO Term')
                plt.xticks(rotation=90)
            
                # Save the plot
                plot_file = os.path.join(
                    folders["plots"]["ko"], 
                    os.path.basename(os.path.splitext(ko_file)[0]) + "_top50.png"
                )
                plt.tight_layout()
                plt.savefig(plot_file, dpi=300)
                plt.close()
            
                logger.info(f"KO top 50 plot saved to: {plot_file}")
            
                # Save the full sorted data
                sorted_ko_file = os.path.join(
                    folders["files"]["ko"], 
                    os.path.basename(os.path.splitext(ko_file)[0]) + "_sorted_by_total.xlsx"
                )
                total_df = pd.DataFrame({'Total': total_counts})
                total_df.to_excel(sorted_ko_file)
                # Also save to original location for backward compatibility
                total_df.to_excel(os.path.join(analysis_dir, os.path.basename(os.path.splitext(ko_file)[0]) + "_sorted_by_total.xlsx"))
                logger.info(f"KO sorted by total count saved to: {sorted_ko_file}")
            
            except Exception as e:
                logger.error(f"Failed to generate KO plot for {ko_file}: {str(e)}")
                logger.error(traceback.format_exc())

def create_ko_treemap(ko_df, output_dir):
    """
    Create a treemap visualization for KEGG Orthology (KO) terms.
    
    Parameters:
    -----------
    ko_df : DataFrame
        DataFrame containing KO term counts
    output_dir : str
        Directory to save the output file
    
    Returns:
    --------
    str
        Path to the saved figure
    """
    try:
        from matplotlib.patches import Rectangle
        import matplotlib.colors as mcolors
        
        # Create figure with explicit axes
        fig, ax = plt.subplots(figsize=(16, 12))
        
        # Limit to top 100 for visibility
        top_ko_df = ko_df.copy()
        if len(top_ko_df) > 100:
            top_ko_df = top_ko_df.head(100)
        
        # Normalize data for size and color
        sizes = top_ko_df['Gene_Count'].values
        sizes = sizes / sizes.sum() * 100  # Normalize to percentage
        
        # Color mapping
        norm = mcolors.Normalize(vmin=top_ko_df['Species_Count'].min(), vmax=top_ko_df['Species_Count'].max())
        colors = plt.cm.viridis(norm(top_ko_df['Species_Count']))
        
        # Calculate treemap layout (simple algorithm)
        labels = top_ko_df['KO'].values
        
        # Simple treemap layout algorithm
        def squarify(sizes, x, y, width, height):
            areas = sizes * width * height / 100
            rects = []
            
            # Sort by size
            order = np.argsort(sizes)[::-1]
            sizes = sizes[order]
            labels_ordered = [labels[i] for i in order]
            colors_ordered = [colors[i] for i in order]
            
            # Simple division algorithm
            remaining_x, remaining_y = x, y
            remaining_width, remaining_height = width, height
            
            for i, (size, label, color) in enumerate(zip(sizes, labels_ordered, colors_ordered)):
                if i < len(sizes) - 1:
                    # For all but the last rectangle
                    if remaining_width > remaining_height:
                        # Divide horizontally
                        rect_width = size * remaining_width / sum(sizes[i:])
                        rect = [remaining_x, remaining_y, rect_width, remaining_height]
                        remaining_x += rect_width
                        remaining_width -= rect_width
                    else:
                        # Divide vertically
                        rect_height = size * remaining_height / sum(sizes[i:])
                        rect = [remaining_x, remaining_y, remaining_width, rect_height]
                        remaining_y += rect_height
                        remaining_height -= rect_height
                else:
                    # Last rectangle takes all remaining space
                    rect = [remaining_x, remaining_y, remaining_width, remaining_height]
                    
                rects.append((rect, label, color))
            
            return rects
        
        # Create the treemap
        rects = squarify(sizes, 0, 0, 1, 1)
        
        # Plot rectangles and labels
        for (rect, label, color), size in zip(rects, sizes):
            x, y, width, height = rect
            ax.add_patch(Rectangle((x, y), width, height, facecolor=color, edgecolor='white', linewidth=1))
            
            # Add text if rectangle is big enough
            if width * height > 0.01:  # Only label rectangles that are large enough
                rx, ry = x + width/2, y + height/2
                fontsize = min(10, 100 * width * height)
                
                # Truncate label if too long
                if len(label) > 20:
                    short_label = label[:17] + "..."
                else:
                    short_label = label
                    
                ax.text(rx, ry, f"{short_label}\n({size:.1f}%)", 
                      ha='center', va='center', fontsize=fontsize,
                      color='white' if np.mean(color[:3]) < 0.5 else 'black')
        
        # Remove axes
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        # Title
        ax.set_title('KEGG Orthology (KO) Terms Treemap', fontweight='bold', fontsize=16)
        
        # Add a colorbar for species count
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label('Number of Species')
        
        # Save figure
        plot_file = os.path.join(output_dir, "treemap_kegg_orthology.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return plot_file
    except Exception as e:
        logger.error(f"Failed to create KO treemap: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_ko_bubble_plot(ko_df, output_dir):
    """
    Create a bubble plot for KEGG Orthology (KO) distribution.
    
    Parameters:
    -----------
    ko_df : DataFrame
        DataFrame containing KO term counts
    output_dir : str
        Directory to save the output file
    
    Returns:
    --------
    str
        Path to the saved figure
    """
    try:
        # Limit to top 50 KOs for visibility
        top_ko_df = ko_df.copy()
        if len(top_ko_df) > 50:
            top_ko_df = top_ko_df.head(50)
        
        plt.figure(figsize=(14, 10))
        
        # Create bubble plot
        scatter = plt.scatter(
            top_ko_df['Gene_Count'], 
            top_ko_df['Species_Count'], 
            s=top_ko_df['Gene_Count'] * 10,  # Bubble size proportional to gene count
            alpha=0.6, 
            c=range(len(top_ko_df)),  # Color by position in the sorted list
            cmap='viridis'
        )
        
        # Add KO labels
        for i, row in top_ko_df.iterrows():
            plt.text(row['Gene_Count'], row['Species_Count'], row['KO'], 
                   fontsize=8, ha='center', va='center',
                   bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        
        # Styling
        plt.title('Top 50 KEGG Orthology (KO) Terms - Gene Count vs Species Count', fontweight='bold', fontsize=14)
        plt.xlabel('Number of Genes', fontsize=12)
        plt.ylabel('Number of Species', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save figure
        plot_file = os.path.join(output_dir, "ko_bubble_plot.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return plot_file
    except Exception as e:
        logger.error(f"Failed to create KO bubble plot: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_kegg_sankey_diagram(merged_kegg_df, output_dir):
    """
    Create a Sankey diagram showing flow from KEGG Level 1 -> 2 -> 3.
    
    Parameters:
    -----------
    merged_kegg_df : DataFrame
        DataFrame containing KEGG hierarchy information
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    dict
        Dictionary with paths to saved HTML and PNG files
    """
    try:
        import plotly.graph_objects as go
        from plotly.offline import plot
    except ImportError:
        logger.warning("Plotly not available. Installing required packages...")
        try:
            import subprocess
            subprocess.check_call(["pip", "install", "plotly", "kaleido"])
            import plotly.graph_objects as go
            from plotly.offline import plot
        except Exception as e:
            logger.error(f"Failed to install required packages: {str(e)}")
            return None
    
    try:
        # Check if we have the required columns
        if not all(col in merged_kegg_df.columns for col in ['Level_1', 'Level_2', 'Level_3']):
            logger.warning("Missing required KEGG hierarchy columns for Sankey diagram")
            return None
        
        # Filter out rows with missing hierarchy information
        valid_df = merged_kegg_df.dropna(subset=['Level_1', 'Level_2', 'Level_3'])
        
        if valid_df.empty:
            logger.warning("No valid data for Sankey diagram after filtering")
            return None
        
        # Get numeric columns for counts
        numeric_cols = valid_df.select_dtypes(include=[np.number]).columns
        
        # Sum counts across samples
        if len(numeric_cols) > 0:
            valid_df['total_count'] = valid_df[numeric_cols].sum(axis=1)
        else:
            logger.warning("No numeric count columns found for Sankey diagram")
            return None
        
        # Create lists for Sankey diagram
        labels = []
        source = []
        target = []
        value = []
        
        # Dictionary to map KEGG terms to node indices
        node_indices = {}
        
        # Process Level 1 -> Level 2 flows
        level1_level2 = valid_df.groupby(['Level_1', 'Level_2'])['total_count'].sum().reset_index()
        
        for _, row in level1_level2.iterrows():
            l1 = row['Level_1']
            l2 = row['Level_2']
            
            # Add nodes if they don't exist
            if l1 not in node_indices:
                node_indices[l1] = len(node_indices)
                labels.append(l1)
            
            if l2 not in node_indices:
                node_indices[l2] = len(node_indices)
                labels.append(l2)
            
            # Add flow
            source.append(node_indices[l1])
            target.append(node_indices[l2])
            value.append(row['total_count'])
        
        # Process Level 2 -> Level 3 flows
        level2_level3 = valid_df.groupby(['Level_2', 'Level_3'])['total_count'].sum().reset_index()
        
        for _, row in level2_level3.iterrows():
            l2 = row['Level_2']
            l3 = row['Level_3']
            
            # Add nodes if they don't exist
            if l2 not in node_indices:
                node_indices[l2] = len(node_indices)
                labels.append(l2)
            
            if l3 not in node_indices:
                node_indices[l3] = len(node_indices)
                labels.append(l3)
            
            # Add flow
            source.append(node_indices[l2])
            target.append(node_indices[l3])
            value.append(row['total_count'])
        
        # Create the Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node = dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=labels,
                color="blue"
            ),
            link = dict(
                source=source,
                target=target,
                value=value
            )
        )])
        
        # Update layout
        fig.update_layout(
            title_text="KEGG Pathway Hierarchy Flow",
            font_size=10,
            width=1200,
            height=900
        )
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Save as interactive HTML
        html_file = os.path.join(output_dir, "kegg_sankey_diagram.html")
        plot(fig, filename=html_file, auto_open=False)
        
        # Save as static image
        png_file = os.path.join(output_dir, "kegg_sankey_diagram.png")
        fig.write_image(png_file)
        
        return {'html': html_file, 'png': png_file}
    
    except Exception as e:
        logger.error(f"Failed to create KEGG Sankey diagram: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_kegg_level_network(merged_kegg_df, output_dir):
    """
    Create network visualizations of KEGG pathway hierarchies.
    
    Parameters:
    -----------
    merged_kegg_df : DataFrame
        DataFrame containing KEGG hierarchy information
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    dict
        Dictionary with paths to saved network visualization files
    """
    try:
        import networkx as nx
    except ImportError:
        logger.warning("NetworkX not available. Installing required packages...")
        try:
            import subprocess
            subprocess.check_call(["pip", "install", "networkx"])
            import networkx as nx
        except Exception as e:
            logger.error(f"Failed to install required packages: {str(e)}")
            return None
    
    try:
        # Check if we have the required columns
        if not all(col in merged_kegg_df.columns for col in ['Level_1', 'Level_2', 'Level_3']):
            logger.warning("Missing required KEGG hierarchy columns for network visualization")
            return None
        
        # Filter out rows with missing hierarchy information
        valid_df = merged_kegg_df.dropna(subset=['Level_1', 'Level_2', 'Level_3'])
        
        if valid_df.empty:
            logger.warning("No valid data for network visualization after filtering")
            return None
        
        # Get numeric columns for counts
        numeric_cols = valid_df.select_dtypes(include=[np.number]).columns
        
        # Sum counts across samples
        if len(numeric_cols) > 0:
            valid_df['total_count'] = valid_df[numeric_cols].sum(axis=1)
        else:
            logger.warning("No numeric count columns found for network visualization")
            return None
        
        # Create a directed graph
        G = nx.DiGraph()
        
        # Add nodes and edges
        for level_num in range(1, 3):  # Level 1->2 and Level 2->3
            src_level = f'Level_{level_num}'
            dst_level = f'Level_{level_num+1}'
            
            # Add nodes with their respective level as attribute
            for level in [src_level, dst_level]:
                for term in valid_df[level].unique():
                    if not pd.isna(term):
                        # Get the count for this term
                        term_count = valid_df[valid_df[level] == term]['total_count'].sum()
                        G.add_node(term, level=level, count=term_count)
            
            # Add edges with weights
            edge_weights = valid_df.groupby([src_level, dst_level])['total_count'].sum().reset_index()
            
            for _, row in edge_weights.iterrows():
                src = row[src_level]
                dst = row[dst_level]
                weight = row['total_count']
                
                if not pd.isna(src) and not pd.isna(dst):
                    G.add_edge(src, dst, weight=weight)
        
        # Create the network visualization
        plt.figure(figsize=(20, 16))
        
        # Try to use hierarchical layout
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            # Fall back to spring layout if graphviz is not available
            pos = nx.spring_layout(G, seed=42)
        
        # Get node attributes for visualization
        node_levels = nx.get_node_attributes(G, 'level')
        node_counts = nx.get_node_attributes(G, 'count')
        
        # Define colors for different levels
        level_colors = {
            'Level_1': 'skyblue',
            'Level_2': 'lightgreen',
            'Level_3': 'salmon'
        }
        
        # Draw nodes with size based on count and color based on level
        for level, color in level_colors.items():
            nodelist = [node for node, node_level in node_levels.items() if node_level == level]
            if nodelist:
                node_sizes = [node_counts.get(node, 1) * 10 for node in nodelist]
                nx.draw_networkx_nodes(
                    G, pos, 
                    nodelist=nodelist,
                    node_size=node_sizes,
                    node_color=color,
                    alpha=0.8,
                    label=f"KEGG {level}"
                )
        
        # Draw edges with width based on weight
        edgelist = [(u, v) for u, v in G.edges()]
        edge_weights = [G[u][v]['weight'] / 10 for u, v in edgelist]
        
        nx.draw_networkx_edges(
            G, pos, 
            edgelist=edgelist,
            width=edge_weights,
            alpha=0.6,
            edge_color='gray',
            arrows=True,
            arrowsize=15,
            connectionstyle='arc3,rad=0.1'
        )
        
        # Draw labels for nodes above a certain count threshold
        count_threshold = np.percentile(list(node_counts.values()), 70)  # Label top 30% of nodes by count
        labels = {node: node for node, count in node_counts.items() if count >= count_threshold}
        
        nx.draw_networkx_labels(
            G, pos, 
            labels=labels,
            font_size=8,
            font_color='black',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )
        
        plt.title('KEGG Pathway Hierarchy Network', fontsize=16, fontweight='bold')
        plt.legend(scatterpoints=1)
        plt.axis('off')
        plt.tight_layout()
        
        # Save the visualization
        output_file = os.path.join(output_dir, "kegg_hierarchy_network.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create a more focused network visualization showing top pathways and their relationships
        try:
            # Create a subgraph of the top nodes by count
            top_nodes = sorted(node_counts.items(), key=lambda x: x[1], reverse=True)[:50]
            top_node_names = [node for node, _ in top_nodes]
            
            # Create subgraph
            SG = G.subgraph(top_node_names)
            
            plt.figure(figsize=(20, 16))
            
            # Use a different layout for the subgraph
            pos_sg = nx.spring_layout(SG, seed=42, k=0.5)
            
            # Get node attributes for the subgraph
            sg_node_levels = {node: node_levels.get(node, 'Unknown') for node in SG.nodes()}
            sg_node_counts = {node: node_counts.get(node, 1) for node in SG.nodes()}
            
            # Draw nodes with size based on count and color based on level
            for level, color in level_colors.items():
                sg_nodelist = [node for node, node_level in sg_node_levels.items() if node_level == level]
                if sg_nodelist:
                    sg_node_sizes = [sg_node_counts.get(node, 1) * 20 for node in sg_nodelist]
                    nx.draw_networkx_nodes(
                        SG, pos_sg, 
                        nodelist=sg_nodelist,
                        node_size=sg_node_sizes,
                        node_color=color,
                        alpha=0.8,
                        label=f"KEGG {level}"
                    )
            
            # Draw edges with width based on weight
            sg_edgelist = [(u, v) for u, v in SG.edges()]
            sg_edge_weights = [SG[u][v]['weight'] / 5 for u, v in sg_edgelist]
            
            nx.draw_networkx_edges(
                SG, pos_sg, 
                edgelist=sg_edgelist,
                width=sg_edge_weights,
                alpha=0.6,
                edge_color='gray',
                arrows=True,
                arrowsize=20,
                connectionstyle='arc3,rad=0.1'
            )
            
            # Draw labels for all nodes in the subgraph
            sg_labels = {node: node for node in SG.nodes()}
            
            nx.draw_networkx_labels(
                SG, pos_sg, 
                labels=sg_labels,
                font_size=10,
                font_color='black',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
            )
            
            plt.title('Top 50 KEGG Pathways Network', fontsize=16, fontweight='bold')
            plt.legend(scatterpoints=1)
            plt.axis('off')
            plt.tight_layout()
            
            # Save the visualization
            top_output_file = os.path.join(output_dir, "top_kegg_pathways_network.png")
            plt.savefig(top_output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            return {
                'full_network': output_file,
                'top_network': top_output_file
            }
            
        except Exception as e:
            logger.error(f"Failed to create top pathways network: {str(e)}")
            logger.error(traceback.format_exc())
            return {'full_network': output_file}
    
    except Exception as e:
        logger.error(f"Failed to create KEGG network visualization: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_circular_kegg_plot(merged_kegg_df, output_dir):
    """
    Create a circular chord diagram showing relationships between KEGG levels.
    
    Parameters:
    -----------
    merged_kegg_df : DataFrame
        DataFrame containing KEGG hierarchy information
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    str
        Path to the saved circular plot
    """
    try:
        # Check if we have the required columns
        if not all(col in merged_kegg_df.columns for col in ['Level_1', 'Level_2', 'Level_3']):
            logger.warning("Missing required KEGG hierarchy columns for circular plot")
            return None
        
        # Filter out rows with missing hierarchy information
        valid_df = merged_kegg_df.dropna(subset=['Level_1', 'Level_2'])
        
        if valid_df.empty:
            logger.warning("No valid data for circular plot after filtering")
            return None
        
        # Get numeric columns for counts
        numeric_cols = valid_df.select_dtypes(include=[np.number]).columns
        
        # Sum counts across samples
        if len(numeric_cols) > 0:
            valid_df['total_count'] = valid_df[numeric_cols].sum(axis=1)
        else:
            logger.warning("No numeric count columns found for circular plot")
            return None
        
        # Create a figure
        plt.figure(figsize=(20, 20))
        
        # Get top Level 1 and Level 2 categories by count
        top_level1 = valid_df.groupby('Level_1')['total_count'].sum().nlargest(8).index.tolist()
        filtered_level2 = valid_df[valid_df['Level_1'].isin(top_level1)]
        top_level2 = filtered_level2.groupby('Level_2')['total_count'].sum().nlargest(20).index.tolist()
        
        # Filter data to only include top categories
        plot_data = valid_df[
            (valid_df['Level_1'].isin(top_level1)) & 
            (valid_df['Level_2'].isin(top_level2))
        ]
        
        # Create a circular layout
        ax = plt.subplot(111, polar=True)
        
        # Set up the figure
        ax.set_theta_direction(-1)  # Clockwise
        ax.set_theta_offset(np.pi/2)  # Start at top
        
        # Get unique categories
        level1_categories = sorted(plot_data['Level_1'].unique())
        level2_categories = sorted(plot_data['Level_2'].unique())
        
        # Calculate positions
        level1_positions = {}
        level2_positions = {}
        
        # Calculate positions for Level 1 (outer circle)
        n_level1 = len(level1_categories)
        level1_width = 2 * np.pi / n_level1
        
        for i, category in enumerate(level1_categories):
            start_angle = i * level1_width
            end_angle = (i + 0.8) * level1_width  # Leave a small gap
            level1_positions[category] = (start_angle, end_angle)
        
        # Calculate positions for Level 2 (inner circle)
        n_level2 = len(level2_categories)
        level2_width = 2 * np.pi / n_level2
        
        for i, category in enumerate(level2_categories):
            start_angle = i * level2_width
            end_angle = (i + 0.8) * level2_width  # Leave a small gap
            level2_positions[category] = (start_angle, end_angle)
        
        # Draw Level 1 categories (outer circle)
        outer_radius = 1.0
        inner_radius = 0.8
        
        cmap = plt.cm.tab20
        colors = cmap(np.linspace(0, 1, len(level1_categories)))
        
        for i, (category, (start, end)) in enumerate(level1_positions.items()):
            # Draw bar
            bars = ax.bar(
                [np.mean([start, end])],  # Center position
                [outer_radius - inner_radius],  # Height (radial)
                width=(end - start),  # Width (angular)
                bottom=inner_radius,  # Start radius
                color=colors[i],
                edgecolor='white',
                linewidth=1,
                alpha=0.7
            )
            
            # Add label
            angle = np.mean([start, end])
            x = np.cos(angle) * (outer_radius + 0.1)
            y = np.sin(angle) * (outer_radius + 0.1)
            
            # Convert to text coordinates
            ha = 'left' if x > 0 else 'right'
            va = 'bottom' if y > 0 else 'top'
            
            plt.text(
                angle, outer_radius + 0.05,
                category,
                ha=ha, va='center',
                rotation=np.degrees(angle) - 90 if x < 0 else np.degrees(angle) + 90,
                fontsize=12,
                fontweight='bold'
            )
        
        # Draw Level 2 categories (inner circle)
        outer_radius_inner = 0.7
        inner_radius_inner = 0.5
        
        cmap_inner = plt.cm.tab20c
        colors_inner = cmap_inner(np.linspace(0, 1, len(level2_categories)))
        
        for i, (category, (start, end)) in enumerate(level2_positions.items()):
            # Draw bar
            bars = ax.bar(
                [np.mean([start, end])],  # Center position
                [outer_radius_inner - inner_radius_inner],  # Height (radial)
                width=(end - start),  # Width (angular)
                bottom=inner_radius_inner,  # Start radius
                color=colors_inner[i],
                edgecolor='white',
                linewidth=1,
                alpha=0.7
            )
        
        # Draw connections between Level 1 and Level 2
        connections = plot_data.groupby(['Level_1', 'Level_2'])['total_count'].sum().reset_index()
        
        max_count = connections['total_count'].max()
        min_count = connections['total_count'].min()
        
        for _, row in connections.iterrows():
            level1 = row['Level_1']
            level2 = row['Level_2']
            count = row['total_count']
            
            # Skip very small connections for clarity
            if count < max_count * 0.05:
                continue
            
            # Get positions
            angle1 = np.mean(level1_positions[level1])
            angle2 = np.mean(level2_positions[level2])
            
            # Calculate curve points
            theta = np.linspace(angle1, angle2, 50)
            
            # Use count to determine curve height and width
            normalized_count = (count - min_count) / (max_count - min_count + 0.001)
            curve_height = 0.75 - normalized_count * 0.15
            
            # Create a radial curve
            r = np.linspace(inner_radius, outer_radius_inner, 50)
            r = r - (r - 0.75) * (r - 0.75) * 4 * normalized_count
            
            # Set curve width proportional to count
            linewidth = 0.5 + 3 * normalized_count
            
            # Draw the curve
            ax.plot(theta, r, linewidth=linewidth, alpha=0.5, color='gray')
        
        # Add title and remove axis ticks
        plt.title('KEGG Pathway Hierarchy - Circular Visualization', y=1.08, fontsize=18, fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['polar'].set_visible(False)
        
        # Save the visualization
        output_file = os.path.join(output_dir, "kegg_circular_plot.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    except Exception as e:
        logger.error(f"Failed to create circular KEGG plot: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_advanced_kegg_heatmap(merged_kegg_df, output_dir):
    """
    Create advanced hierarchical clustered heatmaps for KEGG pathways.
    
    Parameters:
    -----------
    merged_kegg_df : DataFrame
        DataFrame containing KEGG hierarchy information
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    dict
        Dictionary with paths to saved heatmap visualizations
    """
    try:
        import seaborn as sns
    except ImportError:
        logger.warning("Seaborn not available. Installing required packages...")
        try:
            import subprocess
            subprocess.check_call(["pip", "install", "seaborn"])
            import seaborn as sns
        except Exception as e:
            logger.error(f"Failed to install required packages: {str(e)}")
            return None
    
    try:
        outputs = {}
        
        # Check if we have the required columns
        if not all(col in merged_kegg_df.columns for col in ['Level_1', 'Level_2', 'Level_3']):
            logger.warning("Missing required KEGG hierarchy columns for advanced heatmap")
            return None
        
        # Get numeric columns for counts
        numeric_cols = merged_kegg_df.select_dtypes(include=[np.number]).columns
        
        if len(numeric_cols) == 0:
            logger.warning("No numeric count columns found for advanced heatmap")
            return None
        
        # Create heatmaps for each KEGG level
        for level_num, level_col in enumerate(['Level_1', 'Level_2', 'Level_3'], 1):
            # Filter out rows with missing values for the current level
            valid_df = merged_kegg_df.dropna(subset=[level_col])
            
            if valid_df.empty:
                logger.warning(f"No valid data for KEGG {level_col} heatmap")
                continue
            
            # Create a pivot table with the level as rows and numeric columns as values
            pivot_data = valid_df.groupby(level_col)[numeric_cols].sum()
            
            # Normalize the data for better visualization (log scale)
            pivot_data_log = np.log1p(pivot_data)
            
            # Create the heatmap
            plt.figure(figsize=(max(12, len(pivot_data.columns)), max(10, len(pivot_data) * 0.4)))
            
            # Use clustermap for hierarchical clustering
            cg = sns.clustermap(
                pivot_data_log,
                cmap='YlOrRd',
                figsize=(max(12, len(pivot_data.columns)), max(10, len(pivot_data) * 0.4)),
                dendrogram_ratio=(0.2, 0.3),
                cbar_pos=(0.02, 0.8, 0.05, 0.18),
                cbar_kws={'label': 'Log Count'},
                xticklabels=1,
                yticklabels=1,
                linewidths=0.5
            )
            
            # Rotate x-axis labels for readability
            plt.setp(cg.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
            
            # Set title
            cg.fig.suptitle(f'Hierarchical Clustering of KEGG {level_col} Categories', 
                         fontsize=16, fontweight='bold', y=0.98)
            
            # Adjust layout
            cg.fig.tight_layout(rect=[0, 0, 1, 0.96])
            
            # Save the heatmap
            heatmap_file = os.path.join(output_dir, f"kegg_{level_col}_clustered_heatmap.png")
            cg.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            outputs[f'{level_col}_clustered_heatmap'] = heatmap_file
            
            # Also create a standard heatmap with additional annotations
            plt.figure(figsize=(max(12, len(pivot_data.columns)), max(10, len(pivot_data) * 0.4)))
            
            ax = sns.heatmap(
                pivot_data,
                cmap='YlOrRd',
                annot=True,
                fmt='.0f',
                linewidths=0.5,
                cbar_kws={'label': 'Count'}
            )
            
            plt.title(f'KEGG {level_col} Categories Heatmap', fontsize=16, fontweight='bold')
            plt.xlabel('Sample')
            plt.ylabel(f'KEGG {level_col}')
            plt.tight_layout()
            
            # Save the heatmap
            standard_heatmap_file = os.path.join(output_dir, f"kegg_{level_col}_heatmap.png")
            plt.savefig(standard_heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            outputs[f'{level_col}_standard_heatmap'] = standard_heatmap_file
        
        # Create a combined multi-level heatmap if possible
        try:
            # Get top categories from each level based on total counts
            level1_counts = merged_kegg_df.groupby('Level_1')[numeric_cols[0]].sum().nlargest(5).reset_index()
            level2_counts = merged_kegg_df.groupby('Level_2')[numeric_cols[0]].sum().nlargest(10).reset_index()
            level3_counts = merged_kegg_df.groupby('Level_3')[numeric_cols[0]].sum().nlargest(15).reset_index()
            
            # Create a combined DataFrame with hierarchical index
            filtered_df = merged_kegg_df[
                (merged_kegg_df['Level_1'].isin(level1_counts['Level_1'])) &
                (merged_kegg_df['Level_2'].isin(level2_counts['Level_2'])) &
                (merged_kegg_df['Level_3'].isin(level3_counts['Level_3']))
            ]
            
            if not filtered_df.empty:
                # Sum counts by the hierarchical levels
                hierarchical_data = filtered_df.groupby(['Level_1', 'Level_2', 'Level_3'])[numeric_cols].sum()
                
                # Get the first numeric column for a simple heatmap
                if len(numeric_cols) > 1:
                    # Create multiple columns by summing them
                    hierarchical_data['total'] = hierarchical_data.sum(axis=1)
                    heatmap_data = hierarchical_data.total.unstack(level=2)
                else:
                    heatmap_data = hierarchical_data[numeric_cols[0]].unstack(level=2)
                
                # Fill NAs
                heatmap_data = heatmap_data.fillna(0)
                
                # Create the heatmap
                plt.figure(figsize=(20, 12))
                ax = sns.heatmap(
                    heatmap_data,
                    cmap='YlOrRd',
                    linewidth=0.5,
                    cbar_kws={'label': 'Count'}
                )
                
                plt.title('Multi-level KEGG Pathway Hierarchy Heatmap', fontsize=16, fontweight='bold')
                plt.tight_layout()
                
                # Save the heatmap
                multi_level_file = os.path.join(output_dir, "kegg_multi_level_heatmap.png")
                plt.savefig(multi_level_file, dpi=300, bbox_inches='tight')
                plt.close()
                
                outputs['multi_level_heatmap'] = multi_level_file
        except Exception as e:
            logger.error(f"Failed to create multi-level heatmap: {str(e)}")
            logger.error(traceback.format_exc())
        
        return outputs
    
    except Exception as e:
        logger.error(f"Failed to create advanced KEGG heatmaps: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_species_distribution_plots(species_df, output_dir):
    """
    Create distribution plots of genes across species.
    
    Parameters:
    -----------
    species_df : DataFrame
        DataFrame containing genes for specified species
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    dict
        Dictionary with paths to saved species plots
    """
    try:
        outputs = {}
        
        # Check if Species column exists
        if 'Species' not in species_df.columns:
            logger.warning("Species column not found in the data")
            return outputs
        
        # Count genes per species
        species_counts = species_df['Species'].value_counts()
        
        # Create species distribution barplot
        plt.figure(figsize=(12, max(8, len(species_counts) * 0.4)))
        
        # Create barplot with gradient color
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(species_counts)))
        bars = plt.barh(species_counts.index, species_counts.values, color=colors)
        
        # Add count labels
        for i, v in enumerate(species_counts.values):
            plt.text(v + 0.1, i, str(v), va='center')
        
        # Styling
        plt.title('Gene Distribution Across Species', fontweight='bold', fontsize=14)
        plt.xlabel('Number of Genes', fontsize=12)
        plt.ylabel('Species', fontsize=12)
        plt.gca().invert_yaxis()  # To have highest count at the top
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save figure
        barplot_file = os.path.join(output_dir, "species_gene_distribution.png")
        plt.savefig(barplot_file, dpi=300, bbox_inches='tight')
        plt.close()
        outputs['barplot'] = barplot_file
        
        # Create a pie chart for the top species
        plt.figure(figsize=(12, 12))
        
        # Limit to top 10 species for readability
        top_species = species_counts.head(10)
        other_count = species_counts.iloc[10:].sum() if len(species_counts) > 10 else 0
        
        # Add "Other" category if needed
        if other_count > 0:
            top_species = pd.concat([top_species, pd.Series([other_count], index=['Other'])])
        
        # Create pie chart
        wedges, texts, autotexts = plt.pie(
            top_species,
            labels=top_species.index,
            autopct='%1.1f%%',
            startangle=90,
            shadow=False,
            wedgeprops={'edgecolor': 'white'},
            textprops={'fontsize': 14}
        )
        
        # Styling
        plt.title('Species Distribution (Top 10)', fontweight='bold', fontsize=16)
        plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
        
        # Save figure
        pie_file = os.path.join(output_dir, "species_distribution_pie.png")
        plt.savefig(pie_file, dpi=300, bbox_inches='tight')
        plt.close()
        outputs['pie'] = pie_file
        
        # Create a donut chart with annotations
        plt.figure(figsize=(14, 14))
        
        # Create donut chart
        wedges, texts = plt.pie(
            top_species,
            startangle=90,
            wedgeprops=dict(width=0.5, edgecolor='white'),
            colors=plt.cm.tab20(range(len(top_species)))
        )
        
        # Add species names and gene counts
        bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(arrowprops=dict(arrowstyle="-"),
                bbox=bbox_props, zorder=0, va="center")
        
        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1)/2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = f"angle,angleA=0,angleB={ang}"
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            
            plt.annotate(
                f"{top_species.index[i]}: {top_species.iloc[i]} genes",
                xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw
            )
        
        plt.title('Species Distribution with Gene Counts', fontweight='bold', fontsize=16)
        
        # Add center circle for donut
        centre_circle = plt.Circle((0,0),0.70,fc='white')
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)
        
        # Add centered text in the middle
        plt.text(0, 0, f'Total: {species_counts.sum()} genes', 
               ha='center', va='center', fontsize=16, fontweight='bold')
        
        plt.axis('equal')
        
        # Save figure
        donut_file = os.path.join(output_dir, "species_distribution_donut.png")
        plt.savefig(donut_file, dpi=300, bbox_inches='tight')
        plt.close()
        outputs['donut'] = donut_file
        
        # Create a treemap for species distribution
        try:
            from matplotlib.patches import Rectangle
            
            # Create figure
            fig, ax = plt.subplots(figsize=(14, 10))
            
            # Calculate percentages for each species
            total_genes = species_counts.sum()
            percentages = species_counts / total_genes * 100
            
            # Sort by count
            percentages = percentages.sort_values(ascending=False)
            
            # Limit to top 20 species for readability
            top_percentages = percentages.head(20)
            other_percentage = percentages.iloc[20:].sum() if len(percentages) > 20 else 0
            
            # Add "Other" category if needed
            if other_percentage > 0:
                top_percentages = pd.concat([top_percentages, pd.Series([other_percentage], index=['Other'])])
            
            # Create treemap layout
            values = top_percentages.values
            labels = top_percentages.index
            
            # Simple treemap layout algorithm
            def squarify(sizes, x, y, width, height):
                areas = sizes * width * height / 100
                rects = []
                
                # Sort by size
                order = np.argsort(sizes)[::-1]
                sizes = sizes[order]
                labels_ordered = [labels[i] for i in order]
                
                # Simple division algorithm
                remaining_x, remaining_y = x, y
                remaining_width, remaining_height = width, height
                
                for i, (size, label) in enumerate(zip(sizes, labels_ordered)):
                    if i < len(sizes) - 1:
                        # For all but the last rectangle
                        if remaining_width > remaining_height:
                            # Divide horizontally
                            rect_width = size * remaining_width / sum(sizes[i:])
                            rect = [remaining_x, remaining_y, rect_width, remaining_height]
                            remaining_x += rect_width
                            remaining_width -= rect_width
                        else:
                            # Divide vertically
                            rect_height = size * remaining_height / sum(sizes[i:])
                            rect = [remaining_x, remaining_y, remaining_width, rect_height]
                            remaining_y += rect_height
                            remaining_height -= rect_height
                    else:
                        # Last rectangle takes all remaining space
                        rect = [remaining_x, remaining_y, remaining_width, remaining_height]
                        
                    rects.append((rect, label, size))
                
                return rects
            
            # Create the treemap
            rects = squarify(top_percentages.values, 0, 0, 1, 1)
            
            # Plot rectangles and labels
            cmap = plt.cm.viridis
            
            for i, (rect, label, size) in enumerate(rects):
                x, y, width, height = rect
                color = cmap(i / len(rects))
                
                ax.add_patch(Rectangle((x, y), width, height, facecolor=color, edgecolor='white', linewidth=1))
                
                # Add text if rectangle is big enough
                if width * height > 0.02:  # Only label rectangles that are large enough
                    rx, ry = x + width/2, y + height/2
                    fontsize = min(12, 200 * width * height)
                    
                    ax.text(rx, ry, f"{label}\n{size:.1f}%", 
                          ha='center', va='center', fontsize=fontsize,
                          color='white' if i > len(rects)/2 else 'black')
            
            # Remove axes
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect('equal')
            ax.axis('off')
            
            plt.title('Species Distribution Treemap (% of Total Genes)', fontweight='bold', fontsize=16)
            
            # Save figure
            treemap_file = os.path.join(output_dir, "species_distribution_treemap.png")
            plt.savefig(treemap_file, dpi=300, bbox_inches='tight')
            plt.close()
            outputs['treemap'] = treemap_file
            
        except Exception as e:
            logger.error(f"Failed to create species treemap: {str(e)}")
            logger.error(traceback.format_exc())
        
        return outputs
    
    except Exception as e:
        logger.error(f"Failed to create species distribution plots: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def create_html_report(analysis_dir, folders, outputs):
    """
    Create an HTML summary report of all analyses.
    
    Parameters:
    -----------
    analysis_dir : str
        Path to the analysis directory
    folders : dict
        Dictionary containing paths to the created folders
    outputs : dict
        Dictionary of output file paths from all analyses
    
    Returns:
    --------
    str
        Path to the saved report
    """
    try:
        report_file = os.path.join(folders["reports"], "functional_analysis_report.html")
        
        # Create HTML content
        html_content = """<!DOCTYPE html>
<html>
<head>
    <title>Enhanced Functional Annotation Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 30px; line-height: 1.6; }
        h1 { color: #2C3E50; }
        h2 { color: #3498DB; margin-top: 30px; }
        h3 { color: #E74C3C; }
        h4 { color: #9B59B6; }
        img { max-width: 100%; border: 1px solid #ddd; margin: 20px 0; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
        .plot-container { margin-bottom: 40px; }
        .summary-box { background-color: #f8f9fa; border: 1px solid #ddd; padding: 15px; border-radius: 5px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .tab { overflow: hidden; border: 1px solid #ccc; background-color: #f1f1f1; }
        .tab button { background-color: inherit; float: left; border: none; outline: none; cursor: pointer; padding: 14px 16px; transition: 0.3s; }
        .tab button:hover { background-color: #ddd; }
        .tab button.active { background-color: #ccc; }
        .tabcontent { display: none; padding: 6px 12px; border: 1px solid #ccc; border-top: none; }
        .container { display: flex; flex-wrap: wrap; }
        .plot-item { margin: 10px; flex: 1 1 45%; min-width: 300px; max-width: 650px; }
        .plot-item img { width: 100%; height: auto; }
        .plot-item-full { margin: 10px; flex: 1 1 100%; }
    </style>
    <script>
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tabcontent");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tablinks = document.getElementsByClassName("tablinks");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
    </script>
</head>
<body>
    <h1>Enhanced Functional Annotation Analysis Report</h1>
"""
        
        # Add summary information
        html_content += """
    <div class="summary-box">
        <h2>Analysis Summary</h2>
        <p>This report presents the results of an enhanced functional annotation analysis pipeline that combines multiple visualization and analysis techniques.</p>
        <p>The pipeline processes functional annotations from EggNOG-mapper results, with a focus on KEGG Orthology (KO) terms, KEGG Pathways, and COG categories.</p>
        <p>Advanced visualizations include treemaps, network diagrams, Sankey diagrams, clustered heatmaps, and circular plots.</p>
    </div>
"""
        
        # Add tabs for different sections
        html_content += """
    <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'COG')" id="defaultOpen">COG Analysis</button>
        <button class="tablinks" onclick="openTab(event, 'KO')">KEGG Orthology</button>
        <button class="tablinks" onclick="openTab(event, 'KEGGPathways')">KEGG Pathways</button>
        <button class="tablinks" onclick="openTab(event, 'NetworkViz')">Network Visualizations</button>
        <button class="tablinks" onclick="openTab(event, 'Heatmaps')">Heatmaps</button>
        <button class="tablinks" onclick="openTab(event, 'Species')">Species Analysis</button>
        <button class="tablinks" onclick="openTab(event, 'AdvancedViz')">Advanced Visualizations</button>
        <button class="tablinks" onclick="openTab(event, 'DataFiles')">Data Files</button>
    </div>
"""
        
        # COG Analysis tab
        html_content += """
    <div id="COG" class="tabcontent">
        <h2>COG Functional Categories Analysis</h2>
        <p>Clusters of Orthologous Groups (COG) provide functional classification of proteins. This section shows the distribution of COG categories in the dataset.</p>
        
        <div class="container">
"""
        
        # Add COG plots
        cog_plots = [f for f in os.listdir(folders["plots"]["cog"]) if f.endswith('.png')]
        for plot in sorted(cog_plots):
            plot_path = os.path.join(folders["plots"]["cog"], plot)
            plot_rel_path = os.path.relpath(plot_path, analysis_dir)
            
            html_content += f"""
            <div class="plot-item">
                <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                <img src="{plot_rel_path}" alt="{plot}">
            </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # KEGG Orthology tab
        html_content += """
    <div id="KO" class="tabcontent">
        <h2>KEGG Orthology (KO) Analysis</h2>
        <p>KEGG Orthology (KO) terms represent molecular functions in the KEGG database. This section shows the distribution of KO terms in the dataset.</p>
        
        <div class="container">
"""
        
        # Add KO plots
        ko_plots = [f for f in os.listdir(folders["plots"]["ko"]) if f.endswith('.png')]
        for plot in sorted(ko_plots):
            plot_path = os.path.join(folders["plots"]["ko"], plot)
            plot_rel_path = os.path.relpath(plot_path, analysis_dir)
            
            html_content += f"""
            <div class="plot-item">
                <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                <img src="{plot_rel_path}" alt="{plot}">
            </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # KEGG Pathways tab
        html_content += """
    <div id="KEGGPathways" class="tabcontent">
        <h2>KEGG Pathways Analysis</h2>
        <p>KEGG Pathways represent molecular interaction networks for cellular processes. This section shows the distribution of KEGG pathways in the dataset.</p>
        
        <div class="container">
"""
        
        # Add KEGG plots
        kegg_plots = [f for f in os.listdir(folders["plots"]["kegg"]) if f.endswith('.png')]
        for plot in sorted(kegg_plots):
            plot_path = os.path.join(folders["plots"]["kegg"], plot)
            plot_rel_path = os.path.relpath(plot_path, analysis_dir)
            
            html_content += f"""
            <div class="plot-item">
                <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                <img src="{plot_rel_path}" alt="{plot}">
            </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # Network Visualizations tab
        html_content += """
    <div id="NetworkViz" class="tabcontent">
        <h2>Network Visualizations</h2>
        <p>Network visualizations show relationships between different functional categories. This section includes network diagrams for KEGG pathways.</p>
        
        <div class="container">
"""
        
        # Add Network plots if they exist
        if os.path.exists(folders["plots"]["network"]):
            network_plots = [f for f in os.listdir(folders["plots"]["network"]) if f.endswith('.png')]
            for plot in sorted(network_plots):
                plot_path = os.path.join(folders["plots"]["network"], plot)
                plot_rel_path = os.path.relpath(plot_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item-full">
                    <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                    <img src="{plot_rel_path}" alt="{plot}">
                </div>
"""
        
        # Check if networks exist in outputs dictionary
        if 'networks' in outputs and outputs['networks']:
            for network_name, network_path in outputs['networks'].items():
                network_rel_path = os.path.relpath(network_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item-full">
                    <h3>{network_name.replace('_', ' ').title()}</h3>
                    <img src="{network_rel_path}" alt="{network_name}">
                </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # Heatmaps tab
        html_content += """
    <div id="Heatmaps" class="tabcontent">
        <h2>Heatmap Visualizations</h2>
        <p>Heatmaps provide a visual representation of data density. This section includes various heatmap visualizations of the functional annotations.</p>
        
        <div class="container">
"""
        
        # Add Heatmap plots if they exist
        if os.path.exists(folders["plots"]["heatmap"]):
            heatmap_plots = [f for f in os.listdir(folders["plots"]["heatmap"]) if f.endswith('.png')]
            for plot in sorted(heatmap_plots):
                plot_path = os.path.join(folders["plots"]["heatmap"], plot)
                plot_rel_path = os.path.relpath(plot_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item">
                    <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                    <img src="{plot_rel_path}" alt="{plot}">
                </div>
"""
        
        # Check if heatmaps exist in outputs dictionary
        if 'heatmaps' in outputs and outputs['heatmaps']:
            for heatmap_name, heatmap_path in outputs['heatmaps'].items():
                heatmap_rel_path = os.path.relpath(heatmap_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item">
                    <h3>{heatmap_name.replace('_', ' ').title()}</h3>
                    <img src="{heatmap_rel_path}" alt="{heatmap_name}">
                </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # Species Analysis tab
        html_content += """
    <div id="Species" class="tabcontent">
        <h2>Species Analysis</h2>
        <p>Species analysis visualizes the distribution of genes across different species in the dataset.</p>
        
        <div class="container">
"""
        
        # Add Species plots if they exist
        if os.path.exists(folders["plots"]["species"]):
            species_plots = [f for f in os.listdir(folders["plots"]["species"]) if f.endswith('.png')]
            for plot in sorted(species_plots):
                plot_path = os.path.join(folders["plots"]["species"], plot)
                plot_rel_path = os.path.relpath(plot_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item">
                    <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                    <img src="{plot_rel_path}" alt="{plot}">
                </div>
"""
        
        # Check if species plots exist in outputs dictionary
        if 'species' in outputs and outputs['species']:
            for species_name, species_path in outputs['species'].items():
                species_rel_path = os.path.relpath(species_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item">
                    <h3>{species_name.replace('_', ' ').title()}</h3>
                    <img src="{species_rel_path}" alt="{species_name}">
                </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # Advanced Visualizations tab
        html_content += """
    <div id="AdvancedViz" class="tabcontent">
        <h2>Advanced Visualizations</h2>
        <p>This section includes advanced visualizations such as Sankey diagrams, circular plots, and treemaps.</p>
        
        <div class="container">
"""
        
        # Add Sankey plots if they exist
        if os.path.exists(folders["plots"]["sankey"]):
            sankey_plots = [f for f in os.listdir(folders["plots"]["sankey"]) if f.endswith('.png')]
            for plot in sorted(sankey_plots):
                plot_path = os.path.join(folders["plots"]["sankey"], plot)
                plot_rel_path = os.path.relpath(plot_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item-full">
                    <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                    <img src="{plot_rel_path}" alt="{plot}">
                </div>
"""
        
        # Add Circular plots if they exist
        if os.path.exists(folders["plots"]["circular"]):
            circular_plots = [f for f in os.listdir(folders["plots"]["circular"]) if f.endswith('.png')]
            for plot in sorted(circular_plots):
                plot_path = os.path.join(folders["plots"]["circular"], plot)
                plot_rel_path = os.path.relpath(plot_path, analysis_dir)
                
                html_content += f"""
                <div class="plot-item-full">
                    <h3>{plot.replace('.png', '').replace('_', ' ').title()}</h3>
                    <img src="{plot_rel_path}" alt="{plot}">
                </div>
"""
        
        # Check if advanced visualizations exist in outputs dictionary
        # Sankey diagrams
        if 'sankey' in outputs and outputs['sankey']:
            for key, value in outputs['sankey'].items():
                if key == 'png' and value:
                    sankey_rel_path = os.path.relpath(value, analysis_dir)
                    
                    html_content += f"""
                    <div class="plot-item-full">
                        <h3>KEGG Pathway Sankey Diagram</h3>
                        <img src="{sankey_rel_path}" alt="KEGG Pathway Sankey Diagram">
                    </div>
"""
                elif key == 'html' and value:
                    sankey_html_rel_path = os.path.relpath(value, analysis_dir)
                    
                    html_content += f"""
                    <div class="plot-item-full">
                        <h3>Interactive KEGG Pathway Sankey Diagram</h3>
                        <p><a href="{sankey_html_rel_path}" target="_blank">Open Interactive Sankey Diagram</a></p>
                    </div>
"""
        
        # Circular plots
        if 'circular' in outputs and outputs['circular']:
            circular_rel_path = os.path.relpath(outputs['circular'], analysis_dir)
            
            html_content += f"""
            <div class="plot-item-full">
                <h3>KEGG Pathway Circular Visualization</h3>
                <img src="{circular_rel_path}" alt="KEGG Pathway Circular Visualization">
            </div>
"""
        
        html_content += """
        </div>
    </div>
"""
        
        # Data Files tab
        html_content += """
    <div id="DataFiles" class="tabcontent">
        <h2>Data Files</h2>
        <p>This section provides links to the data files generated during the analysis.</p>
        
        <h3>COG Files</h3>
        <table>
            <tr>
                <th>File Name</th>
                <th>Description</th>
            </tr>
"""
        
        # Add COG files
        cog_files = [f for f in os.listdir(folders["files"]["cog"]) if f.endswith('.xlsx') or f.endswith('.csv')]
        for file in sorted(cog_files):
            file_path = os.path.join(folders["files"]["cog"], file)
            file_rel_path = os.path.relpath(file_path, analysis_dir)
            
            html_content += f"""
            <tr>
                <td><a href="{file_rel_path}">{file}</a></td>
                <td>{file.replace('.xlsx', '').replace('.csv', '').replace('_', ' ').title()}</td>
            </tr>
"""
        
        html_content += """
        </table>
        
        <h3>KO Files</h3>
        <table>
            <tr>
                <th>File Name</th>
                <th>Description</th>
            </tr>
"""
        
        # Add KO files
        ko_files = [f for f in os.listdir(folders["files"]["ko"]) if f.endswith('.xlsx') or f.endswith('.csv')]
        for file in sorted(ko_files):
            file_path = os.path.join(folders["files"]["ko"], file)
            file_rel_path = os.path.relpath(file_path, analysis_dir)
            
            html_content += f"""
            <tr>
                <td><a href="{file_rel_path}">{file}</a></td>
                <td>{file.replace('.xlsx', '').replace('.csv', '').replace('_', ' ').title()}</td>
            </tr>
"""
        
        html_content += """
        </table>
        
        <h3>KEGG Files</h3>
        <table>
            <tr>
                <th>File Name</th>
                <th>Description</th>
            </tr>
"""
        
        # Add KEGG files
        kegg_files = [f for f in os.listdir(folders["files"]["kegg"]) if f.endswith('.xlsx') or f.endswith('.csv') or f.endswith('.tsv')]
        for file in sorted(kegg_files):
            file_path = os.path.join(folders["files"]["kegg"], file)
            file_rel_path = os.path.relpath(file_path, analysis_dir)
            
            html_content += f"""
            <tr>
                <td><a href="{file_rel_path}">{file}</a></td>
                <td>{file.replace('.xlsx', '').replace('.csv', '').replace('.tsv', '').replace('_', ' ').title()}</td>
            </tr>
"""
        
        html_content += """
        </table>
    </div>
"""
        
        # Add script to open the default tab and close the HTML
        html_content += """
    <script>
        // Open the first tab by default
        document.getElementById("defaultOpen").click();
    </script>
</body>
</html>
"""
        
        # Write the HTML content to file
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"HTML report saved to: {report_file}")
        
        return report_file
    
    except Exception as e:
        logger.error(f"Failed to create HTML report: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def extract_species_from_eggnog(eggnog_annotations, output_dir):
    """
    Extract species information from EggNOG annotation files.
    
    Parameters:
    -----------
    eggnog_annotations : str
        Path to directory containing EggNOG annotation files
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    DataFrame
        DataFrame containing genes with their species information
    """
    try:
        species_data = []
        
        # Iterate over annotation files
        annotation_files = glob.glob(os.path.join(eggnog_annotations, "*.emapper.annotations"))
        
        if not annotation_files:
            logger.warning(f"No annotation files found in {eggnog_annotations}")
            return pd.DataFrame()
        
        logger.info(f"Found {len(annotation_files)} annotation files")
        
        for annotation_file in annotation_files:
            try:
                base_name = os.path.basename(annotation_file).split('.emapper.annotations')[0]
                logger.info(f"Processing {base_name} for species information")
                
                # Read the annotation file, skipping comment lines
                with open(annotation_file, 'r') as f:
                    lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]
                
                if not lines:
                    logger.warning(f"No data found in {annotation_file}")
                    continue
                
                # Process each line manually
                for line in lines:
                    parts = line.split('\t')
                    
                    # Make sure we have enough columns
                    if len(parts) < 5:
                        continue
                    
                    # Get required columns
                    gene_id = parts[0]
                    species = parts[1] if len(parts) > 1 else "Unknown"
                    description = parts[5] if len(parts) > 5 else "-"
                    cog_category = parts[6] if len(parts) > 6 else "-"
                    ko_terms = parts[11] if len(parts) > 11 else "-"
                    
                    # Extract KEGG pathway information if available
                    kegg_pathway = "-"
                    if len(parts) > 12:
                        kegg_pathway = parts[12]
                    
                    # Add to data
                    species_data.append({
                        'Geneid': gene_id,
                        'Species': species,
                        'Description': description,
                        'COG_category': cog_category,
                        'KEGG_ko': ko_terms,
                        'KEGG_Pathway': kegg_pathway
                    })
                
                logger.info(f"Extracted species information for {len(species_data)} genes from {base_name}")
                
            except Exception as e:
                logger.error(f"Failed to process annotation file {annotation_file} for species: {str(e)}")
                logger.error(traceback.format_exc())
        
        # Convert to DataFrame
        species_df = pd.DataFrame(species_data)
        
        # Save the species data
        species_file = os.path.join(output_dir, "species_annotations.xlsx")
        species_df.to_excel(species_file, index=False)
        species_df.to_csv(os.path.join(output_dir, "species_annotations.csv"), index=False)
        
        logger.info(f"Species annotations saved to: {species_file}")
        
        return species_df
    
    except Exception as e:
        logger.error(f"Failed to extract species information: {str(e)}")
        logger.error(traceback.format_exc())
        return pd.DataFrame()

def run_enhanced_functional_analysis(project_output, kegg_db_file=None, species_based=False):
    """
    Run enhanced functional annotation analysis on EggNOG-mapper results.
    
    Parameters:
    -----------
    project_output : str
        Path to the project output directory
    kegg_db_file : str, optional
        Path to the KEGG database file (ko.txt)
    species_based : bool
        Whether to perform species-based analysis
    
    Returns:
    --------
    str
        Path to the HTML report
    """
    logger.info("Starting enhanced functional annotation analysis")
    
    # Set up directories
    eggnog_output_dir = os.path.join(project_output, "eggnog_output")
    functional_dir = os.path.join(project_output, "functional_annotations")
    analysis_dir = os.path.join(project_output, "functional_analysis")
    ko_dir = os.path.join(functional_dir, "KOs")
    cog_dir = os.path.join(functional_dir, "COGs")
    
    # Create analysis directory
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Create folder structure for organizing files and plots
    folders = create_folder_structure(analysis_dir)
    
    # Initialize outputs dictionary
    outputs = {}
    
    # Check if required directories exist
    if not os.path.exists(eggnog_output_dir):
        logger.error(f"EggNOG output directory not found: {eggnog_output_dir}")
        return None
    
    if not os.path.exists(functional_dir):
        logger.error(f"Functional annotations directory not found: {functional_dir}")
        return None
    
    # Process EggNOG annotation files to count terms
    logger.info("Counting functional terms from EggNOG annotations")
    cog_counts_df, individual_cog_matrix, ko_counts_df = count_eggnog_terms(eggnog_output_dir, analysis_dir, folders)
    
    # Process KO files
    logger.info("Processing KO annotations")
    if os.path.exists(ko_dir):
        process_ko_annotations(ko_dir, analysis_dir, folders, ko_counts_df)
    else:
        logger.warning(f"KO directory not found: {ko_dir}")
        process_ko_annotations(None, analysis_dir, folders, ko_counts_df)
    
    # Process COG files
    logger.info("Processing COG annotations")
    if os.path.exists(cog_dir):
        process_cog_annotations(cog_dir, analysis_dir, folders, individual_cog_matrix)
    else:
        logger.warning(f"COG directory not found: {cog_dir}")
        process_cog_annotations(None, analysis_dir, folders, individual_cog_matrix)
    
    # Handle KEGG database file
    kegg_db_dest = os.path.join(folders["files"]["kegg"], "ko.txt")
    
    # If KEGG database file is provided and exists, use it
    if kegg_db_file and os.path.exists(kegg_db_file):
        # Copy the KEGG database file to the analysis directory
        import shutil
        shutil.copy2(kegg_db_file, kegg_db_dest)
        logger.info(f"Copied KEGG database file from {kegg_db_file} to {kegg_db_dest}")
    else:
        # If not provided, try to find ko.txt in the project directory
        found_kegg_db = False
        if not kegg_db_file:
            logger.info("No KEGG database file provided, searching in project directory")
            # Look for ko.txt in project directory
            for root, dirs, files in os.walk(project_output):
                if "ko.txt" in files:
                    found_kegg_db = True
                    kegg_db_file = os.path.join(root, "ko.txt")
                    # Copy to analysis directory
                    import shutil
                    shutil.copy2(kegg_db_file, kegg_db_dest)
                    logger.info(f"Found and copied KEGG database file from {kegg_db_file} to {kegg_db_dest}")
                    break
            
            if not found_kegg_db:
                logger.warning("No KEGG database file (ko.txt) found in project directory")
                # Also look in common locations
                common_locations = [
                    os.path.join(project_output, "database"),
                    os.path.join(project_output, "databases"),
                    os.path.join(project_output, "kegg"),
                    os.path.join(project_output, "data"),
                    os.path.join(project_output, "reference")
                ]
                
                for location in common_locations:
                    if os.path.exists(location) and os.path.exists(os.path.join(location, "ko.txt")):
                        found_kegg_db = True
                        kegg_db_file = os.path.join(location, "ko.txt")
                        # Copy to analysis directory
                        import shutil
                        shutil.copy2(kegg_db_file, kegg_db_dest)
                        logger.info(f"Found and copied KEGG database file from {kegg_db_file} to {kegg_db_dest}")
                        break
                
                if not found_kegg_db:
                    logger.warning("KEGG database file not found in common locations")
                    logger.warning("KEGG hierarchy mapping will be skipped")
    
    # Try to map KO terms to KEGG hierarchy if database file exists
    merged_kegg_df = None
    if os.path.exists(kegg_db_dest):
        logger.info("Mapping KO terms to KEGG hierarchy")
        merged_kegg_df = map_ko_to_kegg_hierarchy(analysis_dir, folders)
    else:
        logger.warning(f"KEGG database file not found: {kegg_db_dest}")
        logger.warning("Download or provide a KEGG hierarchy file (ko.txt) for complete analysis")
    
    # Generate basic plots
    logger.info("Generating basic plots")
    generate_plots(analysis_dir, folders)
    
    # Generate enhanced COG plots with main categories and full names
    logger.info("Generating enhanced COG plots with main categories")
    plot_cog_by_main_categories(analysis_dir, folders)
    
    # If species-based analysis is requested, extract species information and create species plots
    species_df = None
    if species_based:
        logger.info("Performing species-based analysis")
        species_df = extract_species_from_eggnog(eggnog_output_dir, folders["files"]["species"])
        
        if not species_df.empty:
            logger.info(f"Extracted species information for {len(species_df)} genes from {len(species_df['Species'].unique())} species")
            
            # Create species distribution plots
            logger.info("Creating species distribution plots")
            species_plots = create_species_distribution_plots(species_df, folders["plots"]["species"])
            outputs['species'] = species_plots
            
            # Create KO analysis by species if merged_kegg_df is available
            if merged_kegg_df is not None:
                # Extract KO terms with species info
                ko_species_df = species_df[['Geneid', 'Species', 'KEGG_ko']].dropna(subset=['KEGG_ko'])
                ko_species_df = ko_species_df[ko_species_df['KEGG_ko'] != '-']
                
                if not ko_species_df.empty:
                    # Clean up KO terms
                    ko_species_df['KO'] = ko_species_df['KEGG_ko'].apply(
                        lambda x: x.replace('ko:', '') if isinstance(x, str) and x.startswith('ko:') else x
                    )
                    
                    # Get species count for each KO
                    ko_species_counts = ko_species_df.groupby('KO')['Species'].nunique().reset_index()
                    ko_species_counts.columns = ['KO', 'Species_Count']
                    
                    # Get gene count for each KO
                    ko_gene_counts = ko_species_df.groupby('KO').size().reset_index()
                    ko_gene_counts.columns = ['KO', 'Gene_Count']
                    
                    # Merge counts
                    ko_counts_species = pd.merge(ko_gene_counts, ko_species_counts, on='KO')
                    
                    # Generate enhanced KO visualizations
                    logger.info("Creating enhanced KO visualizations by species")
                    
                    # Create KO treemap
                    ko_treemap = create_ko_treemap(ko_counts_species, folders["plots"]["ko"])
                    if ko_treemap:
                        outputs.setdefault('ko', {})['treemap'] = ko_treemap
                    
                    # Create KO bubble plot
                    ko_bubble = create_ko_bubble_plot(ko_counts_species, folders["plots"]["ko"])
                    if ko_bubble:
                        outputs.setdefault('ko', {})['bubble'] = ko_bubble
    
    # Create advanced visualizations if merged_kegg_df is available
    if merged_kegg_df is not None:
        # Create KEGG Sankey diagram
        logger.info("Creating KEGG Sankey diagram")
        sankey_output = create_kegg_sankey_diagram(merged_kegg_df, folders["plots"]["sankey"])
        if sankey_output:
            outputs['sankey'] = sankey_output
        
        # Create KEGG network visualization
        logger.info("Creating KEGG network visualization")
        network_output = create_kegg_level_network(merged_kegg_df, folders["plots"]["network"])
        if network_output:
            outputs['networks'] = network_output
        
        # Create KEGG circular plot
        logger.info("Creating KEGG circular plot")
        circular_output = create_circular_kegg_plot(merged_kegg_df, folders["plots"]["circular"])
        if circular_output:
            outputs['circular'] = circular_output
        
        # Create advanced KEGG heatmaps
        logger.info("Creating advanced KEGG heatmaps")
        heatmap_output = create_advanced_kegg_heatmap(merged_kegg_df, folders["plots"]["heatmap"])
        if heatmap_output:
            outputs['heatmaps'] = heatmap_output
    
    # Create HTML report
    logger.info("Creating HTML report")
    report_file = create_html_report(analysis_dir, folders, outputs)
    
    logger.info(f"Functional annotation analysis completed. Results available in {analysis_dir}")
    logger.info(f"HTML report: {report_file}")
    
    return report_file

def main():
    """Main function to run the enhanced functional annotation analysis."""
    parser = argparse.ArgumentParser(description='Enhanced Functional Annotation Analysis Pipeline')
    parser.add_argument('--project_dir', required=True, help='Path to the project output directory')
    parser.add_argument('--kegg_db', required=False, help='Path to the KEGG database file (ko.txt)')
    parser.add_argument('--species_based', action='store_true', help='Whether to perform species-based analysis')
    
    args = parser.parse_args()
    
    # Run the analysis
    report_file = run_enhanced_functional_analysis(args.project_dir, args.kegg_db, args.species_based)
    
    if report_file:
        print(f"Analysis completed successfully. HTML report: {report_file}")
    else:
        print("Analysis failed. Check log for details.")

if __name__ == "__main__":
    main()
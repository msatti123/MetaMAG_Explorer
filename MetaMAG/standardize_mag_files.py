#!/usr/bin/env python3
"""
Standalone script to standardize MAG names across all input files.
Creates new standardized versions of all files in the same directories.
"""

import os
import sys
import pandas as pd
import numpy as np
import json
import logging
from pathlib import Path
import argparse
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MAGFileStandardizer:
    """Standardize MAG names across all annotation files."""
    
    def __init__(self):
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.name_mapping = {}
        self.reverse_mapping = {}
        
    def extract_mag_core(self, mag_name):
        """Extract core MAG identity."""
        if not mag_name or pd.isna(mag_name):
            return None
            
        core = str(mag_name)
        
        # Remove suffixes
        suffixes = ['_fa_faa', '_fa.counts', '.counts', '_fa', '.fa', 
                   '.fasta', '.fna', '_faa', '.faa', '_contigs', '_scaffolds']
        for suffix in suffixes:
            if core.endswith(suffix):
                core = core[:-len(suffix)]
        
        # Standardize bin numbering
        import re
        core = re.sub(r'_bin\.(\d+)', r'_bin_\1', core)
        core = re.sub(r'\.bin\.(\d+)', r'_bin_\1', core)
        
        return core
    
    def build_name_mapping(self, tree_file, cazyme_file, cog_file, kegg_file, taxonomy_file, novel_mags_file):
        """Build comprehensive name mapping from all files."""
        logger.info("Building MAG name mapping from all files...")
        
        all_mags = set()
        core_to_variations = {}
        
        # Extract MAG names from tree file
        if tree_file and os.path.exists(tree_file):
            logger.info(f"Reading tree file: {tree_file}")
            with open(tree_file, 'r') as f:
                tree_content = f.read()
                import re
                # Extract MAG names from tree
                mag_pattern = r'([A-Za-z0-9_]+_bin[._]\d+[A-Za-z0-9_]*)'
                tree_mags = re.findall(mag_pattern, tree_content)
                for mag in tree_mags:
                    core = self.extract_mag_core(mag)
                    if core not in core_to_variations:
                        core_to_variations[core] = []
                    core_to_variations[core].append(mag)
                    all_mags.add(mag)
            logger.info(f"  Found {len(tree_mags)} MAGs in tree")
        
        # Extract from CAZyme file
        if cazyme_file and os.path.exists(cazyme_file):
            logger.info(f"Reading CAZyme file: {cazyme_file}")
            df = pd.read_excel(cazyme_file)
            # First column contains MAG names
            if len(df.columns) > 0:
                mag_col = df.iloc[:, 0]
                for mag in mag_col:
                    if mag and not pd.isna(mag):
                        core = self.extract_mag_core(str(mag))
                        if core not in core_to_variations:
                            core_to_variations[core] = []
                        core_to_variations[core].append(str(mag))
                        all_mags.add(str(mag))
            logger.info(f"  Found {len(mag_col)} MAGs in CAZyme file")
        
        # Extract from COG file
        if cog_file and os.path.exists(cog_file):
            logger.info(f"Reading COG file: {cog_file}")
            df = pd.read_excel(cog_file)
            # MAG names are in columns starting from index 3
            if len(df.columns) > 3:
                mag_cols = df.columns[3:]
                for mag in mag_cols:
                    if mag and not pd.isna(mag):
                        core = self.extract_mag_core(str(mag))
                        if core not in core_to_variations:
                            core_to_variations[core] = []
                        core_to_variations[core].append(str(mag))
                        all_mags.add(str(mag))
            logger.info(f"  Found {len(mag_cols)} MAGs in COG file")
        
        # Extract from KEGG file
        if kegg_file and os.path.exists(kegg_file):
            logger.info(f"Reading KEGG file: {kegg_file}")
            df = pd.read_excel(kegg_file)
            # MAG names are in columns starting from index 1
            if len(df.columns) > 1:
                mag_cols = df.columns[1:]
                for mag in mag_cols:
                    if mag and not pd.isna(mag):
                        core = self.extract_mag_core(str(mag))
                        if core not in core_to_variations:
                            core_to_variations[core] = []
                        core_to_variations[core].append(str(mag))
                        all_mags.add(str(mag))
            logger.info(f"  Found {len(mag_cols)} MAGs in KEGG file")
        
        # Extract from taxonomy file
        if taxonomy_file and os.path.exists(taxonomy_file):
            logger.info(f"Reading taxonomy file: {taxonomy_file}")
            sep = '\t' if taxonomy_file.endswith('.tsv') else ','
            df = pd.read_csv(taxonomy_file, sep=sep)
            # Look for user_genome column
            id_col = None
            for col in ['user_genome', 'genome', 'bin_id', 'name']:
                if col in df.columns:
                    id_col = col
                    break
            if id_col:
                for mag in df[id_col]:
                    if mag and not pd.isna(mag):
                        core = self.extract_mag_core(str(mag))
                        if core not in core_to_variations:
                            core_to_variations[core] = []
                        core_to_variations[core].append(str(mag))
                        all_mags.add(str(mag))
                logger.info(f"  Found {len(df[id_col])} MAGs in taxonomy file")
        
        # Extract from novel MAGs file
        if novel_mags_file and os.path.exists(novel_mags_file):
            logger.info(f"Reading novel MAGs file: {novel_mags_file}")
            with open(novel_mags_file, 'r') as f:
                novel_mags = [line.strip() for line in f if line.strip()]
                for mag in novel_mags:
                    core = self.extract_mag_core(mag)
                    if core not in core_to_variations:
                        core_to_variations[core] = []
                    core_to_variations[core].append(mag)
                    all_mags.add(mag)
            logger.info(f"  Found {len(novel_mags)} novel MAGs")
        
        # Create standardized mapping
        logger.info(f"\nCreating standardized names for {len(core_to_variations)} unique MAGs...")
        sorted_cores = sorted(core_to_variations.keys())
        
        for idx, core in enumerate(sorted_cores, 1):
            std_name = f"MAG{idx:04d}"  # MAG0001, MAG0002, etc.
            
            # Map all variations to standard name
            for variation in set(core_to_variations[core]):
                self.name_mapping[variation] = std_name
                self.reverse_mapping[std_name] = core
            
            # Also map the core itself
            self.name_mapping[core] = std_name
        
        logger.info(f"Created mapping for {len(self.name_mapping)} name variations")
        return self.name_mapping
    
    def standardize_cazyme_file(self, input_file, output_file=None):
        """Standardize CAZyme Excel file."""
        if not os.path.exists(input_file):
            return None
            
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_standardized.xlsx"
        
        logger.info(f"Standardizing CAZyme file: {input_file}")
        df = pd.read_excel(input_file)
        
        # First column contains MAG names
        if len(df.columns) > 0:
            first_col = df.columns[0]
            df['original_name'] = df[first_col].copy()
            
            # Standardize names
            df[first_col] = df[first_col].apply(lambda x: self.name_mapping.get(str(x), str(x)) if not pd.isna(x) else x)
            
            # Save
            df.to_excel(output_file, index=False)
            logger.info(f"  Saved to: {output_file}")
            return output_file
        
        return None
    
    def standardize_cog_file(self, input_file, output_file=None):
        """Standardize COG Excel file."""
        if not os.path.exists(input_file):
            return None
            
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_standardized.xlsx"
        
        logger.info(f"Standardizing COG file: {input_file}")
        df = pd.read_excel(input_file)
        
        # MAG names are in columns starting from index 3
        if len(df.columns) > 3:
            new_columns = list(df.columns[:3])  # Keep first 3 columns as is
            
            # Standardize MAG column names
            for col in df.columns[3:]:
                std_name = self.name_mapping.get(str(col), str(col))
                new_columns.append(std_name)
            
            df.columns = new_columns
            
            # Save
            df.to_excel(output_file, index=False)
            logger.info(f"  Saved to: {output_file}")
            return output_file
        
        return None
    
    def standardize_kegg_file(self, input_file, output_file=None):
        """Standardize KEGG Excel file."""
        if not os.path.exists(input_file):
            return None
            
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_standardized.xlsx"
        
        logger.info(f"Standardizing KEGG file: {input_file}")
        df = pd.read_excel(input_file)
        
        # MAG names are in columns starting from index 1
        if len(df.columns) > 1:
            new_columns = [df.columns[0]]  # Keep first column as is
            
            # Standardize MAG column names
            for col in df.columns[1:]:
                std_name = self.name_mapping.get(str(col), str(col))
                new_columns.append(std_name)
            
            df.columns = new_columns
            
            # Save
            df.to_excel(output_file, index=False)
            logger.info(f"  Saved to: {output_file}")
            return output_file
        
        return None
    
    def standardize_taxonomy_file(self, input_file, output_file=None):
        """Standardize taxonomy TSV file."""
        if not os.path.exists(input_file):
            return None
            
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            ext = '.tsv' if input_file.endswith('.tsv') else '.csv'
            output_file = f"{base}_standardized{ext}"
        
        logger.info(f"Standardizing taxonomy file: {input_file}")
        sep = '\t' if input_file.endswith('.tsv') else ','
        df = pd.read_csv(input_file, sep=sep)
        
        # Find and standardize genome ID column
        id_col = None
        for col in ['user_genome', 'genome', 'bin_id', 'name']:
            if col in df.columns:
                id_col = col
                break
        
        if id_col:
            df['original_name'] = df[id_col].copy()
            df[id_col] = df[id_col].apply(lambda x: self.name_mapping.get(str(x), str(x)) if not pd.isna(x) else x)
        
        # Save
        df.to_csv(output_file, sep=sep, index=False)
        logger.info(f"  Saved to: {output_file}")
        return output_file
    
    def standardize_novel_mags_file(self, input_file, output_file=None):
        """Standardize novel MAGs text file."""
        if not os.path.exists(input_file):
            return None
            
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_standardized.txt"
        
        logger.info(f"Standardizing novel MAGs file: {input_file}")
        
        standardized = []
        with open(input_file, 'r') as f:
            for line in f:
                mag = line.strip()
                if mag:
                    std_name = self.name_mapping.get(mag, mag)
                    standardized.append(std_name)
        
        # Save
        with open(output_file, 'w') as f:
            for mag in standardized:
                f.write(f"{mag}\n")
        
        logger.info(f"  Saved to: {output_file}")
        return output_file

    def save_mapping_table(self, output_dir):
        """Save the mapping table for reference."""
        # Ensure output_dir is a directory, not a file path
        if output_dir and os.path.isfile(output_dir):
            output_dir = os.path.dirname(output_dir)
        elif not output_dir:
            output_dir = '.'
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        mapping_file = os.path.join(output_dir, f"mag_name_mapping_{self.timestamp}.csv")
        
        data = []
        for original, standardized in self.name_mapping.items():
            data.append({
                'original_name': original,
                'standardized_name': standardized,
                'core_identity': self.reverse_mapping.get(standardized, '')
            })
        
        df = pd.DataFrame(data)
        df = df.sort_values('standardized_name')
        df.to_csv(mapping_file, index=False)
        
        logger.info(f"Saved mapping table to: {mapping_file}")
        return mapping_file 
   
    def standardize_tree_file(self, input_file, output_file=None):
        """Standardize tree file by replacing MAG names."""
        if not input_file or not os.path.exists(input_file):
            logger.warning(f"Tree file not found: {input_file}")
            return None
        
        if output_file is None:
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_standardized.tree"
        
        logger.info(f"Standardizing tree file: {input_file}")
        
        # Read tree content
        with open(input_file, 'r') as f:
            tree_content = f.read()
        
        # Count replacements
        replacements = 0
        
        # Sort by length (longest first) to avoid partial replacements
        sorted_mappings = sorted(self.name_mapping.items(), key=lambda x: len(x[0]), reverse=True)
        
        # Replace each original name with standardized name
        for original, standardized in sorted_mappings:
            if original in tree_content:
                # Use word boundaries to avoid partial replacements
                import re
                # Match the name when it's followed by : or ) or , (typical tree format delimiters)
                pattern = re.escape(original) + r'(?=[:,)])'
                tree_content = re.sub(pattern, standardized, tree_content)
                replacements += 1
        
        # Save standardized tree
        with open(output_file, 'w') as f:
            f.write(tree_content)
        
        logger.info(f"  Replaced {replacements} MAG names in tree")
        logger.info(f"  Saved to: {output_file}")
        return output_file
    
    def find_tree_file(self, project_dir):
        """Find the tree file in the project directory."""
        tree_locations = [
            os.path.join(project_dir, "Phylogeny", "gtdbtk", "taxonomy2", "trees"),
            os.path.join(project_dir, "Phylogeny", "gtdbtk"),
            os.path.join(project_dir, "Novel_Mags", "gtdbtk"),
            os.path.join(project_dir, "gtdbtk"),
        ]
        
        for location in tree_locations:
            if os.path.exists(location):
                for file in os.listdir(location):
                    if file.endswith('.tree') or file.endswith('.nwk'):
                        tree_path = os.path.join(location, file)
                        logger.info(f"Found tree file: {tree_path}")
                        return tree_path
        
        logger.warning("No tree file found in standard locations")
        return None    

def main():
    parser = argparse.ArgumentParser(description='Standardize MAG names across all annotation files')
    
    parser.add_argument('--tree-file', help='Path to phylogenetic tree file')
    parser.add_argument('--project-dir', help='Project directory to search for tree file')
    parser.add_argument('--cazyme-file', required=True, help='Path to CAZyme Excel file')
    parser.add_argument('--cog-file', required=True, help='Path to COG Excel file')
    parser.add_argument('--kegg-file', required=True, help='Path to KEGG Excel file')
    parser.add_argument('--taxonomy-file', required=True, help='Path to taxonomy TSV file')
    parser.add_argument('--novel-mags-file', required=True, help='Path to novel MAGs text file')
    parser.add_argument('--output-dir', help='Output directory for mapping table', default='.')
    
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("MAG NAME STANDARDIZATION TOOL")
    logger.info("=" * 80)
    
    # Initialize standardizer
    standardizer = MAGFileStandardizer()
    
    # Find tree file if not provided
    tree_file = args.tree_file
    if not tree_file and args.project_dir:
        tree_file = standardizer.find_tree_file(args.project_dir)
    
    # Build mapping from all files
    standardizer.build_name_mapping(
        tree_file,
        args.cazyme_file,
        args.cog_file,
        args.kegg_file,
        args.taxonomy_file,
        args.novel_mags_file
    )
    
    # Standardize each file
    logger.info("\n" + "=" * 80)
    logger.info("STANDARDIZING FILES")
    logger.info("=" * 80)
    
    new_files = {}
    
    # Standardize tree file if found
    if tree_file:
        new_tree = standardizer.standardize_tree_file(tree_file)
        if new_tree:
            new_files['tree'] = new_tree
    
    new_cazyme = standardizer.standardize_cazyme_file(args.cazyme_file)
    if new_cazyme:
        new_files['cazyme'] = new_cazyme
    
    new_cog = standardizer.standardize_cog_file(args.cog_file)
    if new_cog:
        new_files['cog'] = new_cog
    
    new_kegg = standardizer.standardize_kegg_file(args.kegg_file)
    if new_kegg:
        new_files['kegg'] = new_kegg
    
    new_taxonomy = standardizer.standardize_taxonomy_file(args.taxonomy_file)
    if new_taxonomy:
        new_files['taxonomy'] = new_taxonomy
    
    new_novel = standardizer.standardize_novel_mags_file(args.novel_mags_file)
    if new_novel:
        new_files['novel_mags'] = new_novel
    
    # Save mapping table
    mapping_file = standardizer.save_mapping_table(args.output_dir)
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("STANDARDIZATION COMPLETE")
    logger.info("=" * 80)
    logger.info("\nNew standardized files created:")
    for file_type, file_path in new_files.items():
        logger.info(f"  {file_type}: {file_path}")
    logger.info(f"\nMapping table: {mapping_file}")

if __name__ == "__main__":
    main()
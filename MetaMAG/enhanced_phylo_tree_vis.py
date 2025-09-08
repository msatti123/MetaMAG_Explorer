"""
Enhanced Phylogenetic Tree Visualization - COMPLETE VERSION WITH ALL FUNCTIONS
- Support for multiple functional annotation types (CAZyme, COG, KEGG)
- Automatically creates all visualization types:
  * Multi-functional circular tree with pie charts
  * Single functional circular tree
  * Simplified taxonomic circular tree (phylum & genus only)
  * Rectangular tree with colored strips
- Fixed MAG name matching between functional data and tree tips
- Enhanced publication-quality formatting
- Backward compatible with single annotation type
"""

import os
import subprocess
import pandas as pd
import numpy as np
import re
import time
import json
import logging
import glob
from collections import defaultdict
from MetaMAG import enhanced_phylo_tree_vis_pie as pie_viz

# Try to import from MetaMAG, but handle the case where it might not be available
# Import pie module and set availability flag
try:
    from MetaMAG import enhanced_phylo_tree_vis_pie as pie_viz
    PIE_MODULE_AVAILABLE = True
    print("SUCCESS: Pie module imported")  # Use print to see it immediately
except ImportError as e:
    PIE_MODULE_AVAILABLE = False
    print(f"ERROR: Could not import pie module: {e}")
try:
    from MetaMAG.utils import ensure_directory_exists, run_command
    from MetaMAG.config import setup_r_environment as config_setup_r, get_r_config, detect_r_environment
except ImportError:
    # Fallback implementations
    def ensure_directory_exists(directory):
        """Create directory if it doesn't exist."""
        os.makedirs(directory, exist_ok=True)
# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# NOW log the module status
logger.info(f"PIE_MODULE_AVAILABLE = {PIE_MODULE_AVAILABLE}")  

def run_command(cmd):
    """Run a shell command."""
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
def config_setup_r():
    """Dummy R setup."""
    return True
    
def get_r_config():
    """Dummy R config."""
    return {'env_name': 'base', 'r_bin': 'R'}
    
def detect_r_environment():
    """Dummy R environment detection."""
    return 'base'

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_enhanced_tree_visualization(cpus, memory, time_limit, project_input, project_output, 
                                  taxonomy_source=None, novel_mags_file=None, 
                                  functional_annotations=None, annotation_type="auto",
                                  cazyme_annotations=None, cog_annotations=None, 
                                  kegg_annotations=None):
    """Enhanced phylogenetic tree visualization - automatically creates all visualization types."""
    start_time = time.time()
    logger.info("=" * 60)
    logger.info("STARTING ENHANCED PHYLOGENETIC TREE VISUALIZATION")
    logger.info("=" * 60)
    
    # Validate inputs
    if not project_output or not os.path.exists(project_output):
        logger.error(f"Project output directory does not exist: {project_output}")
        return False
    
    # Set up R environment
    logger.info("Setting up R environment...")
    success = config_setup_r()
    if not success:
        logger.error("Failed to set up R environment")
        return False
    
    r_config = get_r_config()
    if r_config:
        logger.info(f"R environment configured: {r_config['env_name']}")
    
    # Define paths
    phylogeny_dir = os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy2")
    trees_dir = os.path.join(phylogeny_dir, "trees")
    vis_output_dir = os.path.join(phylogeny_dir, "enhanced_visualizations")
    
    # Ensure output directory exists
    try:
        ensure_directory_exists(vis_output_dir)
        logger.info(f"Output directory: {vis_output_dir}")
    except Exception as e:
        logger.error(f"Failed to create output directory: {e}")
        return False
    
    # Find tree files
    logger.info("Searching for phylogenetic tree files...")
    tree_files = find_tree_files(trees_dir, project_output)
    
    if not tree_files:
        logger.error("No tree files found!")
        return False
    
    logger.info(f"Found {len(tree_files)} tree file(s):")
    for tree_file in tree_files:
        logger.info(f"  - {tree_file}")
    
    # Process novel MAGs list
    novel_mags = []
    if novel_mags_file and os.path.exists(novel_mags_file):
        novel_mags = load_novel_mags_list_enhanced(novel_mags_file)
        logger.info(f"Loaded {len(novel_mags)} novel MAGs for highlighting")
        if len(novel_mags) > 0:
            logger.info(f"First few novel MAGs: {novel_mags[:5]}")
    
    # Get tree tip names for matching
    tree_tips = get_tree_tip_names(tree_files[0])
    logger.info(f"Found {len(tree_tips)} tree tip names")
    logger.info(f"Sample tree tips: {tree_tips[:5]}")
    
    # Process functional annotations - UPDATED to handle multiple types
    functional_data = {}
    
    # Handle backward compatibility with old single annotation parameter
    if functional_annotations and os.path.exists(functional_annotations):
        logger.info(f"Loading legacy functional annotations from: {functional_annotations}")
        logger.info(f"Legacy annotation type: {annotation_type}")
        
        try:
            legacy_data = load_functional_annotations_enhanced(
                functional_annotations, annotation_type, project_output, tree_tips
            )
            
            if legacy_data:
                # Determine which type this is based on annotation_type
                if annotation_type in ["dbcan", "auto"]:
                    functional_data['cazyme'] = legacy_data
                    logger.info(f"Loaded legacy CAZyme data for {len(legacy_data)} MAGs")
                elif annotation_type == "eggnog":
                    functional_data['eggnog'] = legacy_data
                    logger.info(f"Loaded legacy EggNOG data for {len(legacy_data)} MAGs")
                elif annotation_type == "kegg":
                    functional_data['kegg'] = legacy_data
                    logger.info(f"Loaded legacy KEGG data for {len(legacy_data)} MAGs")
        except Exception as e:
            logger.error(f"Error loading legacy functional annotations: {e}")
    
    # Load specific annotation types (NEW)
    if cazyme_annotations and os.path.exists(cazyme_annotations):
        logger.info(f"Loading CAZyme annotations from: {cazyme_annotations}")
        cazyme_data = load_cazyme_annotations(cazyme_annotations, tree_tips)
        functional_data['cazyme'] = cazyme_data
        logger.info(f"Loaded CAZyme data for {len(cazyme_data)} MAGs")
    
    if cog_annotations and os.path.exists(cog_annotations):
        logger.info(f"Loading COG annotations from: {cog_annotations}")
        cog_data = load_cog_annotations(cog_annotations, tree_tips)
        functional_data['cog'] = cog_data
        logger.info(f"Loaded COG data for {len(cog_data)} MAGs")
    
    if kegg_annotations and os.path.exists(kegg_annotations):
        logger.info(f"Loading KEGG annotations from: {kegg_annotations}")
        kegg_data = load_kegg_annotations(kegg_annotations, tree_tips)
        functional_data['kegg'] = kegg_data
        logger.info(f"Loaded KEGG data for {len(kegg_data)} MAGs")
    
    # Find or create taxonomy metadata
    metadata_file = None
    
    if taxonomy_source and os.path.exists(taxonomy_source):
        logger.info(f"Using provided taxonomy source: {taxonomy_source}")
        try:
            metadata_file = process_metadata_file(taxonomy_source, 
                                                os.path.join(vis_output_dir, "Metadata_For_Tree.csv"))
            if metadata_file:
                logger.info(f"Successfully processed taxonomy metadata: {metadata_file}")
        except Exception as e:
            logger.error(f"Error processing taxonomy source: {e}")
    
    if not metadata_file:
        logger.info("Searching for taxonomy files in project output...")
        try:
            metadata_file = find_taxonomy_files(project_output, vis_output_dir)
            if metadata_file:
                logger.info(f"Found and processed taxonomy file: {metadata_file}")
        except Exception as e:
            logger.warning(f"Error searching for taxonomy files: {e}")
    
    if not metadata_file:
        logger.warning("No metadata file found. Creating basic metadata from tree file...")
        try:
            metadata_file = create_metadata_from_tree(tree_files[0], vis_output_dir)
            if metadata_file:
                logger.info(f"Created basic metadata: {metadata_file}")
            else:
                logger.error("Failed to create metadata from tree file")
                return False
        except Exception as e:
            logger.error(f"Error creating metadata from tree: {e}")
            return False
    
    # Process each tree file
    all_success = True
    total_files = len(tree_files)
    
    logger.info(f"Processing {total_files} tree file(s) for visualization...")
    logger.info("Will create the following visualization types:")
    logger.info("  1. Simplified taxonomic (phylum & genus rings only)")
    logger.info("  2. Rectangular tree with colored strips")
    if functional_data:
        logger.info("  3. Functional visualizations (based on available annotations)")
    
    for i, tree_file in enumerate(tree_files, 1):
        tree_basename = os.path.basename(tree_file)
        output_prefix = os.path.join(vis_output_dir, os.path.splitext(tree_basename)[0])
        
        logger.info(f"\nProcessing tree {i}/{total_files}: {tree_basename}")
        logger.info("-" * 40)
        
        # ALWAYS create simplified taxonomic visualization
        try:
            logger.info("Creating simplified taxonomic visualization...")
            success = create_simplified_taxonomic_visualization(
                tree_file, metadata_file, output_prefix, novel_mags
            )
            if success:
                logger.info("? Successfully created simplified taxonomic visualization")
            else:
                logger.warning("? Failed to create simplified taxonomic visualization")
                all_success = False
        except Exception as e:
            logger.error(f"? Error creating simplified taxonomic visualization: {e}")
            all_success = False
        
        # ALWAYS create rectangular tree visualization
        try:
            logger.info("Creating rectangular tree visualization...")
            success = create_rectangular_tree_visualization(
                tree_file, metadata_file, output_prefix, novel_mags
            )
            if success:
                logger.info("? Successfully created rectangular tree visualization")
            else:
                logger.warning("? Failed to create rectangular tree visualization")
                all_success = False
        except Exception as e:
            logger.error(f"? Error creating rectangular tree visualization: {e}")
            all_success = False
        
        # Create functional visualizations if data is available
        if functional_data:
            try:
                if len(functional_data) > 1:
                    # Multi-functional visualization
                    logger.info("Creating multi-functional visualization...")
                    success = create_multi_functional_visualization(
                        tree_file, metadata_file, output_prefix, 
                        novel_mags, functional_data
                    )
                    if success:
                        logger.info("? Successfully created multi-functional visualization")
                    else:
                        logger.warning("? Failed to create multi-functional visualization")
                        all_success = False
                else:
                    # Single functional visualization
                    logger.info("Creating single functional visualization...")
                    legacy_functional_data = {}
                    if 'cazyme' in functional_data:
                        legacy_functional_data = functional_data['cazyme']
                    elif functional_data:
                        legacy_functional_data = list(functional_data.values())[0]
                    
                    success = create_publication_quality_visualization(
                        tree_file, metadata_file, output_prefix, 
                        novel_mags, legacy_functional_data, annotation_type
                    )
                    if success:
                        logger.info("? Successfully created single functional visualization")
                    else:
                        logger.warning("? Failed to create single functional visualization")
                        all_success = False
                        
            except Exception as e:
                logger.error(f"? Error creating functional visualization: {e}")
                all_success = False
        
        # List all created files for this tree
        logger.info(f"\nCreated files for {tree_basename}:")
        potential_files = [
            f"{output_prefix}_simplified_taxonomic.pdf",
            f"{output_prefix}_rectangular_tree.pdf",
            f"{output_prefix}_multi_functional.pdf",
            f"{output_prefix}_enhanced_beautiful.pdf",
            f"{output_prefix}_enhanced_publication.pdf"
        ]
        
        for potential_file in potential_files:
            if os.path.exists(potential_file):
                file_size = os.path.getsize(potential_file) / 1024  # KB
                logger.info(f"  ? {os.path.basename(potential_file)} ({file_size:.1f} KB)")
    
    # Final summary
    end_time = time.time()
    duration = end_time - start_time
    
    logger.info("=" * 60)
    if all_success:
        logger.info("? ENHANCED TREE VISUALIZATION COMPLETED SUCCESSFULLY")
    else:
        logger.warning("? ENHANCED TREE VISUALIZATION COMPLETED WITH SOME ERRORS")
    
    logger.info(f"Total processing time: {duration:.2f} seconds")
    logger.info(f"Output directory: {vis_output_dir}")
    
    # Summary statistics
    logger.info("\nVisualization Summary:")
    logger.info(f"  Trees processed: {total_files}")
    logger.info(f"  Visualization types created: simplified_taxonomic, rectangular")
    if functional_data:
        logger.info(f"  Functional visualizations: {'multi-functional' if len(functional_data) > 1 else 'single-functional'}")
    
    if novel_mags:
        logger.info(f"  Novel MAGs highlighted: {len(novel_mags)}")
    for func_type, data in functional_data.items():
        logger.info(f"  MAGs with {func_type.upper()} data: {len(data)}")
    
    logger.info("=" * 60)
    
    return all_success

def find_tree_files(trees_dir, project_output):
    """Find phylogenetic tree files - prioritize standardized versions."""
    standardized_tree_files = []
    regular_tree_files = []
    
    def process_directory(directory):
        """Process a directory to find tree files, separating standardized from regular."""
        standardized = []
        regular = []
        
        if os.path.exists(directory):
            for file in os.listdir(directory):
                # Check if it's a tree file
                if file.endswith(".tree") or file.endswith(".nwk") or "tree" in file.lower():
                    # Skip log files
                    if ".log" in file:
                        continue
                    
                    full_path = os.path.join(directory, file)
                    
                    # Check if it's a standardized tree file
                    if "_standardized" in file:
                        # For standardized files, we don't exclude files with hyphens
                        standardized.append(full_path)
                        logger.info(f"Found standardized tree file: {file}")
                    else:
                        # For regular files, exclude those with hyphens (except for standardized)
                        if "-" not in file:
                            regular.append(full_path)
                            logger.info(f"Found regular tree file: {file}")
        
        return standardized, regular
    
    # First check the primary trees directory
    if trees_dir:
        std, reg = process_directory(trees_dir)
        standardized_tree_files.extend(std)
        regular_tree_files.extend(reg)
    
    # If no standardized trees found, check alternative locations
    if not standardized_tree_files:
        alternative_locations = [
            os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy2", "trees"),
            os.path.join(project_output, "Phylogeny", "gtdbtk"),
            os.path.join(project_output, "Phylogeny"),
            os.path.join(project_output, "Novel_Mags", "gtdbtk"),
            os.path.join(project_output, "gtdbtk"),
            project_output
        ]
        
        for location in alternative_locations:
            std, reg = process_directory(location)
            standardized_tree_files.extend(std)
            regular_tree_files.extend(reg)
            
            # Stop searching if we found standardized trees
            if standardized_tree_files:
                break
    
    # Prioritize standardized trees
    if standardized_tree_files:
        logger.info(f"Using {len(standardized_tree_files)} standardized tree file(s)")
        # Remove duplicates while preserving order
        seen = set()
        unique_standardized = []
        for f in standardized_tree_files:
            if f not in seen:
                seen.add(f)
                unique_standardized.append(f)
        return unique_standardized
    elif regular_tree_files:
        logger.warning("No standardized tree files found, falling back to regular tree files")
        logger.info(f"Using {len(regular_tree_files)} regular tree file(s)")
        # Remove duplicates while preserving order
        seen = set()
        unique_regular = []
        for f in regular_tree_files:
            if f not in seen:
                seen.add(f)
                unique_regular.append(f)
        return unique_regular
    else:
        logger.error("No tree files found (neither standardized nor regular)")
        return []

def find_specific_standardized_tree(project_output, tree_type="bac120"):
    """
    Find a specific standardized tree file (e.g., gtdbtk.bac120_standardized.tree).
    
    Args:
        project_output: Project output directory
        tree_type: Type of tree to find (e.g., "bac120", "ar122")
    
    Returns:
        Path to the standardized tree file if found, None otherwise
    """
    tree_filename = f"gtdbtk.{tree_type}_standardized.tree"
    alternative_tree_filename = f"{tree_type}_standardized.tree"
    
    search_locations = [
        os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy2", "trees"),
        os.path.join(project_output, "Phylogeny", "gtdbtk"),
        os.path.join(project_output, "Novel_Mags", "gtdbtk"),
        os.path.join(project_output, "gtdbtk"),
        project_output
    ]
    
    for location in search_locations:
        if os.path.exists(location):
            # Check for the standard filename
            tree_path = os.path.join(location, tree_filename)
            if os.path.exists(tree_path):
                logger.info(f"Found standardized {tree_type} tree: {tree_path}")
                return tree_path
            
            # Check for alternative filename
            alt_tree_path = os.path.join(location, alternative_tree_filename)
            if os.path.exists(alt_tree_path):
                logger.info(f"Found standardized {tree_type} tree: {alt_tree_path}")
                return alt_tree_path
    
    logger.warning(f"Standardized {tree_type} tree not found")
    return None

def find_taxonomy_files(project_output, output_dir):
    """Find GTDB-Tk taxonomy files."""
    potential_locations = [
        os.path.join(project_output, "Novel_Mags", "gtdbtk"),
        os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy"),
        os.path.join(project_output, "gtdbtk"),
        project_output
    ]
    
    taxonomy_patterns = [
        "gtdbtk.bac120.summary.tsv",
        "gtdbtk.ar122.summary.tsv",
        "*summary.tsv",
        "*.classifications.txt"
    ]
    
    for location in potential_locations:
        if os.path.exists(location):
            for pattern in taxonomy_patterns:
                pattern_path = os.path.join(location, pattern)
                matching_files = glob.glob(pattern_path)
                
                if matching_files:
                    taxonomy_file = matching_files[0]
                    return process_metadata_file(taxonomy_file, os.path.join(output_dir, "Metadata_For_Tree.csv"))
    
    return None

def process_metadata_file(taxonomy_file, output_path):
    """Process GTDB-Tk taxonomy file."""
    try:
        if taxonomy_file.endswith('.tsv'):
            df = pd.read_csv(taxonomy_file, sep='\t')
        else:
            for sep in ['\t', ',']:
                try:
                    df = pd.read_csv(taxonomy_file, sep=sep)
                    break
                except:
                    continue
            else:
                return None
        
        metadata_df = pd.DataFrame()
        
        id_column = None
        for potential_id in ['user_genome', 'genome', 'bin_id', 'name', 'Genome']:
            if potential_id in df.columns:
                id_column = potential_id
                break
        
        if not id_column:
            id_column = df.columns[0]
        
        metadata_df['user_genome'] = df[id_column]
        
        if 'classification' in df.columns:
            taxonomic_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            for level in taxonomic_levels:
                metadata_df[level] = 'Unknown'
            
            for idx, row in df.iterrows():
                if pd.isna(row['classification']):
                    continue
                
                taxa = row['classification'].split(';')
                for taxon in taxa:
                    if taxon.startswith('d__'):
                        metadata_df.at[idx, 'Domain'] = taxon[3:]
                    elif taxon.startswith('p__'):
                        metadata_df.at[idx, 'Phylum'] = taxon[3:]
                    elif taxon.startswith('c__'):
                        metadata_df.at[idx, 'Class'] = taxon[3:]
                    elif taxon.startswith('o__'):
                        metadata_df.at[idx, 'Order'] = taxon[3:]
                    elif taxon.startswith('f__'):
                        metadata_df.at[idx, 'Family'] = taxon[3:]
                    elif taxon.startswith('g__'):
                        metadata_df.at[idx, 'Genus'] = taxon[3:]
                    elif taxon.startswith('s__'):
                        metadata_df.at[idx, 'Species'] = taxon[3:]
        
        for col in metadata_df.columns:
            if col != 'user_genome':
                metadata_df[col] = metadata_df[col].astype(str).apply(
                    lambda x: x[3:] if isinstance(x, str) and len(x) > 3 and x[1:3] == '__' else x
                )
                metadata_df[col] = metadata_df[col].replace(['', 'nan', 'None', 'none'], 'Unknown')
                metadata_df[col] = metadata_df[col].fillna('Unknown')
        
        metadata_df.to_csv(output_path, index=False)
        return output_path
    
    except Exception as e:
        logger.error(f"Error processing metadata file: {e}")
        return None

def create_metadata_from_tree(tree_file, output_dir):
    """Create basic metadata from tree file."""
    output_path = os.path.join(output_dir, "Metadata_From_Tree.csv")
    
    try:
        with open(tree_file, 'r') as f:
            tree_content = f.read()
        
        genome_ids = re.findall(r'([a-zA-Z0-9_.]+)', tree_content)
        genome_ids = [gid for gid in genome_ids if not re.match(r'^\d+(\.\d+)?$', gid)]
        
        unique_ids = []
        for gid in genome_ids:
            if gid not in unique_ids and len(gid) > 1:
                unique_ids.append(gid)
        
        data = {
            'user_genome': unique_ids,
            'Domain': ['Bacteria'] * len(unique_ids),
            'Phylum': ['Unknown'] * len(unique_ids),
            'Class': ['Unknown'] * len(unique_ids),
            'Order': ['Unknown'] * len(unique_ids),
            'Family': ['Unknown'] * len(unique_ids),
            'Genus': ['Unknown'] * len(unique_ids),
            'Species': ['Unknown'] * len(unique_ids)
        }
        
        metadata_df = pd.DataFrame(data)
        metadata_df.to_csv(output_path, index=False)
        return output_path
    
    except Exception as e:
        logger.error(f"Error creating metadata from tree: {e}")
        return None

def run_r_command_with_conda(script_path):
    """Run R script with proper conda environment activation."""
    try:
        r_config = get_r_config()
        
        if not r_config:
            logger.error("No R configuration found")
            return False
        
        # Try different methods to run R
        if r_config.get('r_bin') and os.path.exists(r_config['r_bin']):
            cmd = f"{r_config['r_bin']}script {script_path}"
            logger.info(f"Using direct R path: {r_config['r_bin']}")
        elif r_config['env_name'] != 'system':
            cmd = f"conda run -n {r_config['env_name']} Rscript {script_path}"
            logger.info(f"Using conda run with environment: {r_config['env_name']}")
        else:
            cmd = f"Rscript {script_path}"
            logger.info("Using system Rscript")
        
        logger.info(f"Executing R command: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=False, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                               universal_newlines=True)
        
        for line in result.stdout.split('\n'):
            if line.strip():
                logger.info(f"R OUTPUT: {line}")
        
        for line in result.stderr.split('\n'):
            if line.strip():
                logger.error(f"R ERROR: {line}")
        
        return result.returncode == 0
        
    except Exception as e:
        logger.error(f"Error running R command: {e}")
        return False

def get_tree_tip_names(tree_file):
    """Extract tree tip names to help with MAG name matching."""
    try:
        with open(tree_file, 'r') as f:
            tree_content = f.read()
        
        # Extract potential tip names using patterns for BOTH old and standardized formats
        tip_names = []
        
        # Patterns for OLD format (SRR names with bin)
        old_patterns = [
            r'([A-Z]+\d+_[a-z0-9]+_bin[._]\d+(?:_[a-z]+)?)',
            r'([a-z]+_bin[._]\d+(?:_[a-z]+)?)',
            r'([A-Z]+\d+_[a-z0-9]+_bin[._]\d+)',
        ]
        
        # Pattern for STANDARDIZED format (MAG####)
        standardized_pattern = r'(MAG\d{4})'
        
        # Try standardized pattern first
        import re
        standardized_matches = re.findall(standardized_pattern, tree_content)
        tip_names.extend(standardized_matches)
        
        # Also try old patterns in case of mixed content
        for pattern in old_patterns:
            matches = re.findall(pattern, tree_content)
            tip_names.extend(matches)
        
        # Clean and deduplicate
        clean_tips = []
        seen = set()
        for tip in tip_names:
            if tip not in seen and len(tip) > 3:
                clean_tips.append(tip)
                seen.add(tip)
        
        logger.info(f"Extracted {len(clean_tips)} tree tip names")
        if clean_tips:
            logger.info(f"Sample tips: {clean_tips[:5]}")
        
        return clean_tips
        
    except Exception as e:
        logger.warning(f"Error extracting tree tip names: {e}")
        return []

def load_functional_annotations_enhanced(functional_annotations, annotation_type, project_output, tree_tips):
    """Enhanced functional annotation loading with FIXED MAG name matching."""
    if annotation_type == "auto":
        annotation_type = detect_annotation_type(functional_annotations, project_output)
    
    if annotation_type == "dbcan":
        return load_dbcan_annotations_with_matching(functional_annotations, tree_tips)
    elif annotation_type == "eggnog":
        return load_eggnog_annotations(functional_annotations, tree_tips)
    elif annotation_type == "kegg":
        return load_kegg_annotations(functional_annotations, tree_tips)
    else:
        logger.warning(f"Unknown annotation type: {annotation_type}")
        return {}

def load_eggnog_annotations(eggnog_file, tree_tips):
    """Load EggNOG annotations - placeholder implementation."""
    # This would need to be implemented based on your EggNOG file format
    logger.warning("EggNOG annotation loading not yet implemented")
    return {}

def load_cazyme_annotations(cazyme_file, tree_tips):
    """Load CAZyme annotations from Excel file - FIXED for standard format."""
    try:
        logger.info(f"Reading CAZyme Excel file: {cazyme_file}")
        df = pd.read_excel(cazyme_file)
        logger.info(f"CAZyme file shape: {df.shape}")
        logger.info(f"CAZyme columns: {list(df.columns)}")
        
        # YOUR FILE IS IN STANDARD FORMAT: MAGs in rows, categories in columns
        # Column 0: MAG names, Columns 1+: GT, CE, GH, CBM, PL, AA counts
        
        functional_data = {}
        
        # Get MAG names from first column (Unnamed: 0)
        mag_column = df.columns[0]  # This is 'Unnamed: 0'
        category_columns = ['GT', 'CE', 'GH', 'CBM', 'PL', 'AA']  # These are your data columns
        
        logger.info(f"Processing {len(df)} MAGs from standard format")
        logger.info(f"MAG name column: {mag_column}")
        logger.info(f"Category columns: {category_columns}")
        
        # Process each row (each MAG)
        for idx, row in df.iterrows():
            try:
                # Get and clean MAG name
                raw_mag_name = str(row[mag_column])
                mag_name = clean_mag_column_name(raw_mag_name)
                
                if not mag_name:
                    logger.debug(f"Skipping invalid MAG name: {raw_mag_name}")
                    continue
                
                # Extract category counts for this MAG
                categories = {}
                total_genes = 0
                
                for category in category_columns:
                    if category in row:
                        try:
                            count = int(float(row[category])) if not pd.isna(row[category]) else 0
                            categories[category] = count
                            total_genes += count
                        except (ValueError, TypeError):
                            categories[category] = 0
                    else:
                        categories[category] = 0
                
                # Only add if we have some data
                if total_genes > 0:
                    functional_data[mag_name] = {
                        'categories': categories,  # Direct category counts
                        'total_genes': total_genes,
                        'total_families': len([k for k, v in categories.items() if v > 0])
                    }
                    logger.debug(f"Added {mag_name}: {total_genes} total genes across {len([k for k, v in categories.items() if v > 0])} categories")
                
            except Exception as e:
                logger.warning(f"Error processing row {idx}: {e}")
                continue
        
        logger.info(f"Processed {len(functional_data)} MAGs with CAZyme data")
        
        # Match to tree tips
        matched_data = match_functional_data_to_tree_tips(functional_data, tree_tips)
        return matched_data
        
    except Exception as e:
        logger.error(f"Error loading CAZyme annotations: {e}")
        return {}

def load_cog_annotations(cog_file, tree_tips):
    """Load COG annotations from Excel file."""
    try:
        logger.info(f"Reading COG Excel file: {cog_file}")
        df = pd.read_excel(cog_file)
        logger.info(f"COG file shape: {df.shape}")
        
        # COG data has Main_Category, COG, Description, then MAG columns
        mag_columns = df.columns[3:].tolist()  # Skip first 3 columns
        
        functional_data = {}
        
        for mag_col in mag_columns:
            # Clean MAG name
            mag_name = clean_mag_column_name(mag_col)
            if not mag_name:
                continue
            
            # Extract COG category counts
            cog_counts = {}
            
            for idx, row in df.iterrows():
                try:
                    cog_code = str(row['COG']).strip()
                    count = row[mag_col]
                    
                    if pd.isna(count):
                        count = 0
                    else:
                        count = int(float(count))
                    
                    if count > 0 and cog_code and cog_code != 'nan':
                        cog_counts[cog_code] = count
                        
                except (ValueError, TypeError, KeyError):
                    continue
            
            # Add to functional data
            if cog_counts:
                functional_data[mag_name] = {
                    'categories': cog_counts,
                    'total_genes': sum(cog_counts.values()),
                    'total_categories': len(cog_counts)
                }
        
        # Match to tree tips
        matched_data = match_functional_data_to_tree_tips(functional_data, tree_tips)
        return matched_data
        
    except Exception as e:
        logger.error(f"Error loading COG annotations: {e}")
        return {}

def load_kegg_annotations(kegg_file, tree_tips):
    """Load KEGG annotations from Excel file."""
    try:
        logger.info(f"Reading KEGG Excel file: {kegg_file}")
        df = pd.read_excel(kegg_file)
        logger.info(f"KEGG file shape: {df.shape}")
        
        # KEGG data has Level_1, then MAG columns
        mag_columns = df.columns[1:].tolist()  # Skip first column
        
        functional_data = {}
        
        for mag_col in mag_columns:
            # Clean MAG name
            mag_name = clean_mag_column_name(mag_col)
            if not mag_name:
                continue
            
            # Extract KEGG pathway counts
            kegg_counts = {}
            
            for idx, row in df.iterrows():
                try:
                    pathway = str(row['Level_1']).strip()
                    count = row[mag_col]
                    
                    if pd.isna(count):
                        count = 0
                    else:
                        count = int(float(count))
                    
                    if count > 0 and pathway and pathway != 'nan':
                        # Simplify pathway names for visualization
                        pathway_short = simplify_kegg_pathway_name(pathway)
                        kegg_counts[pathway_short] = count
                        
                except (ValueError, TypeError, KeyError):
                    continue
            
            # Add to functional data
            if kegg_counts:
                functional_data[mag_name] = {
                    'categories': kegg_counts,
                    'total_genes': sum(kegg_counts.values()),
                    'total_pathways': len(kegg_counts)
                }
        
        # Match to tree tips
        matched_data = match_functional_data_to_tree_tips(functional_data, tree_tips)
        return matched_data
        
    except Exception as e:
        logger.error(f"Error loading KEGG annotations: {e}")
        return {}

def simplify_kegg_pathway_name(pathway):
    """Simplify KEGG pathway names for visualization."""
    # Remove numeric codes and simplify
    simplified = pathway.replace('09100', '').replace('09120', '').replace('09130', '')
    simplified = simplified.replace('09140', '').replace('09150', '').replace('09160', '')
    simplified = simplified.replace('09180', '').replace('09190', '').strip()
    
    # Map to shorter names
    name_map = {
        'Metabolism': 'Metabolism',
        'Genetic Information Processing': 'Genetic Info',
        'Environmental Information Processing': 'Environ Info',
        'Cellular Processes': 'Cellular',
        'Organismal Systems': 'Organismal',
        'Human Diseases': 'Diseases',
        'Brite Hierarchies': 'Brite',
        'Not Included in Pathway or Brite': 'Other'
    }
    
    return name_map.get(simplified, simplified)

def clean_mag_column_name(mag_col):
    """Enhanced MAG column name cleaning for your specific naming patterns."""
    if not mag_col or mag_col in ['nan', 'NaN', '']:
        return None
    
    mag_name = str(mag_col)
    original_name = mag_name
    
    # CRITICAL FIX: Handle your specific naming patterns
    # Remove multiple suffix combinations that appear in your data
    
    # Pattern 1: Remove _fa_faa (COG/KEGG pattern: SRR15883085_metabat_bin_3_fa_faa)
    if mag_name.endswith('_fa_faa'):
        mag_name = mag_name[:-7]  # Remove '_fa_faa'
    
    # Pattern 2: Remove _fa.counts (CAZyme pattern: SRR15883098_maxbin2_bin_3_fa.counts)  
    if mag_name.endswith('_fa.counts'):
        mag_name = mag_name[:-10]  # Remove '_fa.counts'
    
    # Pattern 3: Remove .counts
    if mag_name.endswith('.counts'):
        mag_name = mag_name[:-7]  # Remove '.counts'
    
    # Pattern 4: Remove _fa suffix
    if mag_name.endswith('_fa'):
        mag_name = mag_name[:-3]  # Remove '_fa'
    
    # Pattern 5: Remove other common suffixes
    other_suffixes = ['.fa', '.fasta', '.fna', '_faa', '.faa', '_contigs', '_scaffolds']
    for suffix in other_suffixes:
        if mag_name.endswith(suffix):
            mag_name = mag_name[:-len(suffix)]
            break  # Only remove one suffix per iteration
    
    # CRITICAL FIX: Convert bin_NUMBER to bin.NUMBER (your tree format)
    # This converts: SRR15883085_metabat_bin_3 -> SRR15883085_metabat_bin.3
    mag_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
    
    clean_name = mag_name.strip()
    
    # Debug logging to see conversions
    if original_name != clean_name:
        logger.debug(f"Name conversion: {original_name} -> {clean_name}")
    
    return clean_name if clean_name else None

def standardize_cazy_family_name(family_name):
    """Standardize CAZy family names."""
    if not family_name or family_name in ['nan', 'NaN', '']:
        return None
    
    family_name = str(family_name).upper().strip()
    family_name = family_name.replace('_', '').replace('-', '').replace(' ', '')
    
    pattern = r'^(GH|GT|PL|CE|AA|CBM)(\d+)$'
    match = re.match(pattern, family_name)
    
    if match:
        family_type = match.group(1)
        family_number = match.group(2)
        return f"{family_type}{family_number}"
    else:
        for prefix in ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM']:
            if family_name.startswith(prefix):
                number_part = family_name[len(prefix):]
                if number_part.isdigit():
                    return f"{prefix}{number_part}"
        return None

def get_cazy_family_type(family_name):
    """Extract CAZy family type."""
    if not family_name:
        return None
    
    for prefix in ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM']:
        if family_name.startswith(prefix):
            return prefix
    
    return None

def match_functional_data_to_tree_tips(raw_functional_data, tree_tips):
    """Match functional data to tree tips with comprehensive name conversion."""
    matched_data = {}
    
    logger.info(f"Matching {len(raw_functional_data)} functional entries to {len(tree_tips)} tree tips")
    
    # Strategy 1: Direct exact matching
    exact_matches = 0
    for mag_name, data in raw_functional_data.items():
        if mag_name in tree_tips:
            matched_data[mag_name] = data
            exact_matches += 1
    
    logger.info(f"Exact matches: {exact_matches}")
    
    # Strategy 2: Convert underscores to dots for bin numbers
    conversion_matches = 0
    for mag_name, data in raw_functional_data.items():
        if any(mag_name == k.replace('.', '_').replace('_sub', '') for k in matched_data.keys()):
            continue
            
        converted_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
        
        if converted_name in tree_tips:
            matched_data[converted_name] = data
            conversion_matches += 1
            continue
            
        converted_with_sub = converted_name + '_sub'
        if converted_with_sub in tree_tips:
            matched_data[converted_with_sub] = data
            conversion_matches += 1
            continue
            
        for tree_tip in tree_tips:
            if tree_tip in matched_data:
                continue
            if tree_tip.startswith(converted_name + '_'):
                matched_data[tree_tip] = data
                conversion_matches += 1
                break
    
    logger.info(f"Conversion matches: {conversion_matches}")
    
    # Strategy 3: Fuzzy matching for remaining
    fuzzy_matches = 0
    remaining_func = {}
    for mag_name, data in raw_functional_data.items():
        already_matched = False
        for matched_tip in matched_data.keys():
            base_matched = matched_tip.replace('.', '_').replace('_sub', '').replace('_contigs', '')
            base_func = mag_name.replace('_sub', '').replace('_contigs', '')
            if base_matched == base_func:
                already_matched = True
                break
        if not already_matched:
            remaining_func[mag_name] = data
    
    remaining_tips = [tip for tip in tree_tips if tip not in matched_data]
    
    for mag_name, data in remaining_func.items():
        normalized_mag = normalize_mag_name(mag_name)
        
        for tree_tip in remaining_tips[:]:
            if tree_tip in matched_data:
                continue
                
            normalized_tip = normalize_mag_name(tree_tip)
            
            if normalized_mag == normalized_tip or \
               (len(normalized_mag) > 6 and len(normalized_tip) > 6 and 
                (normalized_mag in normalized_tip or normalized_tip in normalized_mag)):
                matched_data[tree_tip] = data
                remaining_tips.remove(tree_tip)
                fuzzy_matches += 1
                break
    
    logger.info(f"Fuzzy matches: {fuzzy_matches}")
    total_matches = exact_matches + conversion_matches + fuzzy_matches
    logger.info(f"Total matched: {total_matches} out of {len(raw_functional_data)}")
    
    return matched_data

def normalize_mag_name(name):
    """Normalize MAG name for matching."""
    if not name:
        return ""
    
    norm = str(name).lower()
    
    # Remove common suffixes
    removals = ['_fa', '_fasta', '.fa', '.fasta', '.counts', '_contigs', '_scaffolds']
    for removal in removals:
        norm = norm.replace(removal, '')
    
    # Standardize separators
    norm = norm.replace('.', '_').replace('-', '_')
    
    # Remove multiple underscores
    while '__' in norm:
        norm = norm.replace('__', '_')
    
    return norm.strip('_')

def load_dbcan_annotations_with_matching(functional_annotations, tree_tips):
    """Load dbCAN annotations with FIXED MAG name matching."""
    functional_data = {}
    
    # Handle directory vs file input
    if os.path.isdir(functional_annotations):
        excel_files = []
        for file in os.listdir(functional_annotations):
            if file.endswith('.xlsx') and any(keyword in file.lower() for keyword in ['cazyme', 'cazy', 'dbcan']):
                excel_files.append(os.path.join(functional_annotations, file))
        
        if not excel_files:
            logger.warning(f"No CAZyme Excel files found in {functional_annotations}")
            return {}
        excel_file = excel_files[0]
    elif functional_annotations.endswith('.xlsx'):
        excel_file = functional_annotations
    else:
        logger.warning(f"Expected directory or Excel file, got: {functional_annotations}")
        return {}
    
    if not os.path.exists(excel_file):
        logger.error(f"Excel file not found: {excel_file}")
        return {}
    
    try:
        logger.info(f"Reading Excel file: {excel_file}")
        df = pd.read_excel(excel_file)
        logger.info(f"Successfully read Excel file with {df.shape[0]} rows, {df.shape[1]} columns")
        
        # Process the Excel data
        raw_functional_data = process_excel_cazy_data_enhanced(df)
        
        if not raw_functional_data:
            logger.warning("No functional data extracted from Excel file")
            return {}
        
        logger.info(f"Extracted functional data for {len(raw_functional_data)} MAGs from Excel")
        
        # FIXED: Better name matching strategy
        functional_data = match_functional_data_to_tree_tips_fixed(raw_functional_data, tree_tips)
        
        logger.info(f"Successfully matched {len(functional_data)} functional data entries to tree tips")
        
        return functional_data
        
    except Exception as e:
        logger.error(f"Error reading Excel file {excel_file}: {e}")
        return {}

def match_functional_data_to_tree_tips_fixed(raw_functional_data, tree_tips):
    """FIXED matching function with comprehensive name conversion strategies."""
    matched_data = {}
    
    logger.info("=== FIXED MATCHING: Functional data to tree tips ===")
    logger.info(f"Raw functional MAGs: {list(raw_functional_data.keys())[:5]}")
    logger.info(f"Tree tips: {tree_tips[:5]}")
    
    # Strategy 1: Direct exact matching
    exact_matches = 0
    for mag_name, data in raw_functional_data.items():
        if mag_name in tree_tips:
            matched_data[mag_name] = data
            exact_matches += 1
            logger.info(f"Exact match: {mag_name}")
    
    logger.info(f"Strategy 1 - Exact matches: {exact_matches}")
    
    # Strategy 2: Convert underscores to dots for bin numbers (CRITICAL FIX)
    conversion_matches = 0
    for mag_name, data in raw_functional_data.items():
        if any(mag_name == k.replace('.', '_').replace('_sub', '') for k in matched_data.keys()):
            continue  # Already matched (check by base name)
            
        # Convert patterns: bin_NUMBER -> bin.NUMBER
        converted_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
        
        # Try direct match with converted name
        if converted_name in tree_tips:
            matched_data[converted_name] = data
            conversion_matches += 1
            logger.info(f"Conversion match: {mag_name} -> {converted_name}")
            continue
            
        # Try match with _sub suffix
        converted_with_sub = converted_name + '_sub'
        if converted_with_sub in tree_tips:
            matched_data[converted_with_sub] = data
            conversion_matches += 1
            logger.info(f"Conversion+sub match: {mag_name} -> {converted_with_sub}")
            continue
            
        # Try other suffix patterns
        for tree_tip in tree_tips:
            if tree_tip in matched_data:
                continue
            # Check if tree tip starts with converted name + underscore (for other suffixes)
            if tree_tip.startswith(converted_name + '_'):
                matched_data[tree_tip] = data
                conversion_matches += 1
                logger.info(f"Conversion+suffix match: {mag_name} -> {tree_tip}")
                break
    
    logger.info(f"Strategy 2 - Conversion matches: {conversion_matches}")
    
    # Strategy 3: Handle coassembly bins (no conversion needed)
    coassembly_matches = 0
    for mag_name, data in raw_functional_data.items():
        if any(mag_name == k.replace('_sub', '') for k in matched_data.keys()):
            continue  # Already matched
            
        if mag_name.startswith('coassembly_bin'):
            # Try direct match
            if mag_name in tree_tips:
                matched_data[mag_name] = data
                coassembly_matches += 1
                logger.info(f"Coassembly exact match: {mag_name}")
                continue
                
            # Try with dot conversion (just in case)
            converted_coassembly = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
            if converted_coassembly in tree_tips:
                matched_data[converted_coassembly] = data
                coassembly_matches += 1
                logger.info(f"Coassembly conversion match: {mag_name} -> {converted_coassembly}")
                continue
    
    logger.info(f"Strategy 3 - Coassembly matches: {coassembly_matches}")
    
    # Strategy 4: Fuzzy matching for remaining unmatched
    fuzzy_matches = 0
    remaining_func = {}
    for mag_name, data in raw_functional_data.items():
        already_matched = False
        for matched_tip in matched_data.keys():
            # Check if this functional data is already matched
            base_matched = matched_tip.replace('.', '_').replace('_sub', '').replace('_contigs', '')
            base_func = mag_name.replace('_sub', '').replace('_contigs', '')
            if base_matched == base_func:
                already_matched = True
                break
        if not already_matched:
            remaining_func[mag_name] = data
    
    remaining_tips = [tip for tip in tree_tips if tip not in matched_data]
    
    logger.info(f"Remaining unmatched functional: {len(remaining_func)} -> {list(remaining_func.keys())[:3]}")
    logger.info(f"Remaining unmatched tree tips: {len(remaining_tips)} -> {remaining_tips[:3]}")
    
    for mag_name, data in remaining_func.items():
        normalized_mag = normalize_mag_name_fixed(mag_name)
        
        for tree_tip in remaining_tips[:]:  # Use slice to avoid modification during iteration
            if tree_tip in matched_data:
                continue
                
            normalized_tip = normalize_mag_name_fixed(tree_tip)
            
            # More flexible similarity check
            if normalized_mag == normalized_tip or \
               (len(normalized_mag) > 6 and len(normalized_tip) > 6 and 
                (normalized_mag in normalized_tip or normalized_tip in normalized_mag)):
                matched_data[tree_tip] = data
                remaining_tips.remove(tree_tip)
                fuzzy_matches += 1
                logger.info(f"Fuzzy match: {mag_name} ({normalized_mag}) -> {tree_tip} ({normalized_tip})")
                break
    
    logger.info(f"Strategy 4 - Fuzzy matches: {fuzzy_matches}")
    
    total_matches = exact_matches + conversion_matches + coassembly_matches + fuzzy_matches
    logger.info(f"TOTAL MATCHED: {total_matches} out of {len(raw_functional_data)} functional entries")
    
    # Debug: Show final matched pairs
    logger.info("=== FINAL MATCHED PAIRS ===")
    for tree_tip, data in list(matched_data.items())[:5]:
        families_count = len([k for k in data.get('cazy_families', {}).keys() if not k.endswith('_total')])
        logger.info(f"  {tree_tip}: {families_count} CAZy families")
    
    return matched_data

def normalize_mag_name_fixed(name):
    """Improved MAG name normalization for matching."""
    if not name:
        return ""
    
    # Convert to lowercase
    norm = str(name).lower()
    
    # Remove common suffixes
    removals = ['_fa', '_fasta', '.fa', '.fasta', '.counts', '_contigs', '_scaffolds']
    for removal in removals:
        norm = norm.replace(removal, '')
    
    # Standardize separators - convert dots to underscores for comparison
    norm = norm.replace('.', '_').replace('-', '_')
    
    # Remove multiple underscores
    while '__' in norm:
        norm = norm.replace('__', '_')
    
    # Remove leading/trailing underscores
    norm = norm.strip('_')
    
    return norm

def process_excel_cazy_data_enhanced(df):
    """Enhanced processing of CAZy data from Excel with better debugging."""
    functional_data = {}
    
    try:
        logger.info(f"Processing Excel data with shape: {df.shape}")
        logger.info(f"First few column names: {list(df.columns)[:10]}")
        
        # Check if data is transposed
        first_col_name = df.columns[0]
        logger.info(f"First column name: '{first_col_name}'")
        
        # Sample first column values
        sample_first_col_values = df[first_col_name].dropna().head(10).tolist()
        logger.info(f"Sample values from first column: {sample_first_col_values}")
        
        # Count CAZy-like values
        cazy_like_count = 0
        for val in sample_first_col_values:
            val_str = str(val).upper()
            if any(val_str.startswith(prefix) for prefix in ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM']):
                cazy_like_count += 1
        
        logger.info(f"Found {cazy_like_count} CAZy-like values in first column out of {len(sample_first_col_values)} samples")
        
        # Determine if transposed
        is_transposed = cazy_like_count > len(sample_first_col_values) * 0.5
        
        if is_transposed:
            logger.info("Detected TRANSPOSED format: CAZy families in rows, MAGs in columns")
            functional_data = process_transposed_cazy_data_enhanced(df)
        else:
            logger.info("Detected STANDARD format: MAGs in rows, CAZy families in columns")
            functional_data = process_standard_cazy_data_enhanced(df)
        
        # Log results
        if functional_data:
            logger.info(f"Successfully processed {len(functional_data)} MAGs from Excel file")
            
            # Show detailed statistics
            total_families = 0
            total_counts = 0
            for mag, data in functional_data.items():
                if 'cazy_families' in data:
                    non_total_families = [k for k in data['cazy_families'].keys() if not k.endswith('_total')]
                    mag_families = len(non_total_families)
                    mag_counts = sum(data['cazy_families'][k] for k in non_total_families)
                    total_families += mag_families
                    total_counts += mag_counts
                    
                    logger.debug(f"  {mag}: {mag_families} families, {mag_counts} total counts")
            
            if len(functional_data) > 0:
                avg_families = total_families / len(functional_data)
                logger.info(f"Average CAZy families per MAG: {avg_families:.1f}")
        else:
            logger.warning("No functional data processed from Excel file")
        
        return functional_data
        
    except Exception as e:
        logger.error(f"Error processing Excel CAZy data: {e}")
        return {}

def process_transposed_cazy_data_enhanced(df):
    """Enhanced processing of transposed CAZy data with better error handling."""
    functional_data = {}
    
    try:
        # Get CAZy family names from first column
        cazy_col = df.columns[0]
        cazy_families = df[cazy_col].dropna().tolist()
        
        logger.info(f"Processing {len(cazy_families)} CAZy families")
        
        # Get MAG columns (all columns except first)
        mag_columns = df.columns[1:].tolist()
        logger.info(f"Found {len(mag_columns)} MAG columns")
        logger.info(f"Sample MAG columns: {mag_columns[:5]}")
        
        # Process each MAG column
        for mag_col in mag_columns:
            try:
                # Clean MAG name
                mag_name = clean_mag_column_name_enhanced(mag_col)
                if not mag_name:
                    logger.debug(f"Skipping invalid MAG column: {mag_col}")
                    continue
                
                logger.debug(f"Processing MAG: {mag_col} -> {mag_name}")
                
                # Extract counts for this MAG
                cazy_counts = {}
                family_type_totals = {}
                
                for idx, family_name in enumerate(cazy_families):
                    try:
                        # Get the count value
                        if idx < len(df):
                            count = df.iloc[idx][mag_col]
                        else:
                            continue
                        
                        # Handle different data types
                        if pd.isna(count):
                            count = 0
                        else:
                            try:
                                count = int(float(count))
                            except (ValueError, TypeError):
                                count = 0
                        
                        if count > 0:
                            # Standardize family name
                            clean_family = standardize_cazy_family_name(str(family_name))
                            
                            if clean_family:
                                cazy_counts[clean_family] = count
                                
                                # Count family type totals
                                family_type = get_cazy_family_type(clean_family)
                                if family_type:
                                    total_key = f"{family_type}_total"
                                    if total_key not in family_type_totals:
                                        family_type_totals[total_key] = 0
                                    family_type_totals[total_key] += count
                                        
                    except (ValueError, TypeError, IndexError) as e:
                        logger.debug(f"Skipping invalid count for {family_name} in {mag_name}: {e}")
                        continue
                
                # Add to functional data if we found any families
                if cazy_counts:
                    # Combine individual families with totals
                    all_counts = {**cazy_counts, **family_type_totals}
                    
                    functional_data[mag_name] = {
                        'cazy_families': all_counts,
                        'total_genes': sum(cazy_counts.values()),
                        'total_cazy_families': len(cazy_counts)
                    }
                    
                    logger.debug(f"  {mag_name}: {len(cazy_counts)} families, {sum(cazy_counts.values())} total counts")
                else:
                    logger.debug(f"No CAZy families found for {mag_name}")
                    
            except Exception as e:
                logger.warning(f"Error processing MAG column {mag_col}: {e}")
                continue
        
        return functional_data
        
    except Exception as e:
        logger.error(f"Error in process_transposed_cazy_data_enhanced: {e}")
        return {}

def process_standard_cazy_data_enhanced(df):
    """Processing of standard CAZy data format."""
    functional_data = {}
    
    try:
        # In standard format, MAGs are in rows and categories are in columns
        mag_column = df.columns[0]
        category_columns = ['GT', 'CE', 'GH', 'CBM', 'PL', 'AA']
        
        logger.info(f"Processing standard format with MAGs in column: {mag_column}")
        
        for idx, row in df.iterrows():
            try:
                mag_name = clean_mag_column_name(str(row[mag_column]))
                if not mag_name:
                    continue
                
                categories = {}
                total_genes = 0
                
                for category in category_columns:
                    if category in row:
                        try:
                            count = int(float(row[category])) if not pd.isna(row[category]) else 0
                            categories[category] = count
                            total_genes += count
                        except (ValueError, TypeError):
                            categories[category] = 0
                
                if total_genes > 0:
                    # Add category totals
                    for cat in categories:
                        categories[f"{cat}_total"] = categories[cat]
                    
                    functional_data[mag_name] = {
                        'cazy_families': categories,
                        'total_genes': total_genes,
                        'total_cazy_families': len([k for k in categories.keys() if not k.endswith('_total')])
                    }
                    
            except Exception as e:
                logger.warning(f"Error processing row {idx}: {e}")
                continue
        
        return functional_data
        
    except Exception as e:
        logger.error(f"Error in process_standard_cazy_data_enhanced: {e}")
        return {}

def clean_mag_column_name_enhanced(mag_col):
    """Enhanced MAG column name cleaning with multiple patterns."""
    if not mag_col or mag_col in ['nan', 'NaN', '']:
        return None
    
    mag_name = str(mag_col)
    original_name = mag_name
    
    # Remove .counts suffix
    if mag_name.endswith('.counts'):
        mag_name = mag_name[:-7]
    
    # Remove _fa suffix  
    if mag_name.endswith('_fa'):
        mag_name = mag_name[:-3]
    
    # Remove file extensions
    for ext in ['.fa', '.fasta', '.fna', '.faa']:
        if mag_name.endswith(ext):
            mag_name = mag_name[:-len(ext)]
    
    # Remove other common suffixes
    for pattern in ['_contigs', '_scaffolds']:
        if mag_name.endswith(pattern):
            mag_name = mag_name[:-len(pattern)]
    
    # Clean up the name
    clean_name = mag_name.strip()
    
    logger.debug(f"Cleaned MAG name: {original_name} -> {clean_name}")
    
    return clean_name if clean_name else None

def load_novel_mags_list_enhanced(novel_mags_file):
    """Enhanced novel MAGs list loading with better validation and name conversion."""
    novel_mags = []
    
    if not os.path.exists(novel_mags_file):
        logger.error(f"Novel MAGs file not found: {novel_mags_file}")
        return []
    
    try:
        with open(novel_mags_file, 'r') as f:
            line_count = 0
            for line_num, line in enumerate(f, 1):
                line_count += 1
                mag_name = line.strip()
                
                # Skip empty lines and comments
                if not mag_name or mag_name.startswith('#'):
                    continue
                
                # Convert to tree format to match tree tip names (bin_NUMBER -> bin.NUMBER)
                tree_format_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
                if tree_format_name:
                    novel_mags.append(tree_format_name)
                    logger.debug(f"Converted novel MAG: {mag_name} -> {tree_format_name}")
        
        logger.info(f"Processed {line_count} lines from novel MAGs file")
        logger.info(f"Loaded {len(novel_mags)} valid novel MAG names")
        
        return novel_mags
        
    except Exception as e:
        logger.error(f"Error reading novel MAGs file: {e}")
        return []

def detect_annotation_type(functional_annotations, project_output):
    """Auto-detect the type of functional annotation."""
    # Check for dbCAN files first (most likely based on your data)
    if functional_annotations.endswith('.xlsx') or 'cazyme' in functional_annotations.lower():
        return "dbcan"
        
    # Check for EggNOG files
    eggnog_paths = [
        os.path.join(project_output, "Functional_Annotation", "EggNOG"),
        os.path.join(project_output, "eggnog"),
        functional_annotations
    ]
    
    for path in eggnog_paths:
        if os.path.exists(path):
            if any("emapper" in f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))):
                return "eggnog"
    
    # Check for KEGG files
    if functional_annotations and ("kegg" in functional_annotations.lower() or "ko" in functional_annotations.lower()):
        return "kegg"
    
    return "dbcan"  # Default to dbcan for Excel files

def create_simplified_taxonomic_visualization(tree_file, metadata_file, output_prefix, novel_mags):
    """Create simplified circular tree with only phylum and genus rings."""
    try:
        script_path = f"{output_prefix}_simplified_taxonomic_vis.R"
        
        logger.info(f"Creating simplified taxonomic visualization")
        logger.info(f"  - {len(novel_mags)} novel MAG entries")
        
        # Prepare data for R script
        novel_mags_json = json.dumps(novel_mags)
        
        r_script = f'''#!/usr/bin/env Rscript

# Simplified Taxonomic Tree Visualization - Only Phylum and Genus Rings

# Load required packages
tryCatch({{
  library(ape)
  library(RColorBrewer)
  library(jsonlite)
  cat("Base packages loaded successfully\\n")
}}, error = function(e) {{
  cat("Error loading packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- fromJSON('{novel_mags_json}')

set.seed(42)

cat("=== Simplified Taxonomic Tree Visualization ===\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")

# Read and process tree
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded successfully with", length(tree$tip.label), "tips\\n")
  
  if (!is.rooted(tree)) {{
    if (length(tree$tip.label) > 2) {{
      outgroup_tip <- tree$tip.label[1]
      tryCatch({{
        tree <- root(tree, outgroup = outgroup_tip)
        cat("Applied root using tip:", outgroup_tip, "\\n")
      }}, error = function(e) {{
        cat("Rooting failed, using tree as-is\\n")
      }})
    }}
  }}
}}, error = function(e) {{
  cat("Error loading tree:", e$message, "\\n")
  quit(status=1)
}})

# Read metadata
tryCatch({{
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns\\n")
  id_col <- colnames(metadata)[1]
}}, error = function(e) {{
  cat("Error loading metadata:", e$message, "\\n")
  quit(status=1)
}})

# Process taxonomy data
taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]])

for (level in taxonomic_levels) {{
  if (level %in% colnames(metadata)) {{
    taxa_data[[level]] <- as.character(metadata[[level]])
    taxa_data[[level]][taxa_data[[level]] %in% c("Unknown", "", "NA", "na", "None")] <- "Unknown"
  }} else {{
    taxa_data[[level]] <- "Unknown"
  }}
}}

# Parse classification column if exists
if ("classification" %in% colnames(metadata)) {{
  cat("Extracting taxonomy from classification column\\n")
  for (i in 1:nrow(metadata)) {{
    if (!is.na(metadata$classification[i])) {{
      taxa <- strsplit(metadata$classification[i], ";")[[1]]
      for (taxon in taxa) {{
        if (startsWith(taxon, "d__")) {{
          taxa_data$Domain[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "p__")) {{
          taxa_data$Phylum[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "c__")) {{
          taxa_data$Class[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "o__")) {{
          taxa_data$Order[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "f__")) {{
          taxa_data$Family[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "g__")) {{
          taxa_data$Genus[i] <- substring(taxon, 4)  
        }} else if (startsWith(taxon, "s__")) {{
          taxa_data$Species[i] <- substring(taxon, 4)
        }}
      }}
    }}
  }}
}}

# Match tree tips with taxa data
match_tip_to_taxa <- function(tip, taxa_data) {{
  match_idx <- which(taxa_data$ID == tip)
  if (length(match_idx) == 0) {{
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  if (length(match_idx) == 0) {{
    converted_tip <- gsub("\\\\.", "_", tip)
    match_idx <- which(taxa_data$ID == converted_tip)
  }}
  if (length(match_idx) == 0) {{
    for (i in 1:nrow(taxa_data)) {{
      base_taxa <- gsub("_sub$", "", taxa_data$ID[i])
      base_tip <- gsub("_sub$", "", tip)
      if (base_taxa == base_tip || 
          grepl(base_taxa, base_tip, fixed=TRUE) || 
          grepl(base_tip, base_taxa, fixed=TRUE)) {{
        match_idx <- i
        break
      }}
    }}
  }}
  
  if (length(match_idx) > 0) {{
    return(taxa_data[match_idx[1], ])
  }} else {{
    unknown_row <- taxa_data[1, ]
    unknown_row[] <- c(tip, rep("Unknown", length(taxonomic_levels)))
    return(unknown_row)
  }}
}}

tip_data <- do.call(rbind, lapply(tree$tip.label, function(tip) {{
  return(match_tip_to_taxa(tip, taxa_data))
}}))

# Add novel MAG information
tip_data$is_novel <- FALSE

for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  if (tip %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    next
  }}
  
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    next
  }}
  
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
  }}
}}

novel_count <- sum(tip_data$is_novel)
cat("Novel MAGs matched:", novel_count, "\\n")

# Create color palettes
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

if (length(phyla) > 0) {{
  if (length(phyla) <= 12) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Paired")[1:length(phyla)]
  }} else {{
    phylum_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(phyla))
  }}
  names(phylum_colors) <- phyla
  phylum_colors <- c(phylum_colors, "Unknown" = "#E0E0E0")
}} else {{
  phylum_colors <- c("Unknown" = "#E0E0E0")
}}

if (length(genera) > 0) {{
  if (length(genera) <= 8) {{
    genus_colors <- brewer.pal(max(3, length(genera)), "Dark2")[1:length(genera)]
  }} else {{
    genus_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "#F0F0F0")
}} else {{
  genus_colors <- c("Unknown" = "#F0F0F0")
}}

# Create maps
tip_to_phylum <- setNames(tip_data$Phylum, tip_data$ID)
tip_to_genus <- setNames(tip_data$Genus, tip_data$ID)
tip_to_novel <- setNames(tip_data$is_novel, tip_data$ID)

# Helper functions
get_taxon <- function(tip_label, taxon_map) {{
  if (tip_label %in% names(taxon_map)) {{
    return(taxon_map[tip_label])
  }} else {{
    return("Unknown")
  }}
}}

get_color <- function(value, palette) {{
  if (is.null(value) || is.na(value) || value == "" || value == "Unknown") {{
    return(palette["Unknown"])
  }} else if (value %in% names(palette)) {{
    return(palette[value])
  }} else {{
    return(palette["Unknown"])
  }}
}}

draw_arc_segment <- function(center_x, center_y, radius_inner, radius_outer, 
                           start_angle, end_angle, color) {{
  if (start_angle > end_angle) {{
    tmp <- start_angle
    start_angle <- end_angle
    end_angle <- tmp
  }}
  
  if ((end_angle - start_angle) < 0.01) {{
    end_angle <- start_angle + 0.01
  }}
  
  angles <- seq(start_angle, end_angle, length.out = max(100, ceiling((end_angle - start_angle) * 200/pi)))
  
  x_inner <- center_x + radius_inner * cos(angles)
  y_inner <- center_y + radius_inner * sin(angles)
  x_outer <- center_x + radius_outer * cos(rev(angles))
  y_outer <- center_y + radius_outer * sin(rev(angles))
  
  polygon(c(x_inner, x_outer), c(y_inner, y_outer), 
          col = color, border = "white", lwd = 0.5)
}}

##############################################
# CREATE SIMPLIFIED TAXONOMIC VISUALIZATION
##############################################

simplified_file <- paste0(output_prefix, "_simplified_taxonomic.pdf")
pdf(simplified_file, width = 16, height = 14)

# Clean layout with minimal elements
layout(matrix(c(1,2), nrow=1), widths=c(3,1))
par(mar=c(1,1,3,1), oma=c(0,0,0,0))

# Draw tree
plot(tree, type="fan", show.tip.label=FALSE, 
     main="Phylogenetic Tree - Taxonomic Classification", 
     edge.width=2.5, cex.main=1.8, edge.color="#333333", 
     no.margin=FALSE, use.edge.length=FALSE,
     x.lim=c(-1.5, 1.5), y.lim=c(-1.5, 1.5))

# Get coordinates
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords_x <- last_plot$xx[1:length(tree$tip.label)]
tip_coords_y <- last_plot$yy[1:length(tree$tip.label)]

center_x <- mean(range(tip_coords_x))
center_y <- mean(range(tip_coords_y))
max_radius <- max(sqrt((tip_coords_x - center_x)^2 + (tip_coords_y - center_y)^2))

# Calculate angles
tip_angles <- atan2(tip_coords_y - center_y, tip_coords_x - center_x)
tip_angles[tip_angles < 0] <- tip_angles[tip_angles < 0] + 2 * pi

# Sort tips by angle
sorted_indices <- order(tip_angles)
sorted_tips <- tree$tip.label[sorted_indices]
sorted_angles <- tip_angles[sorted_indices]
sorted_tips <- c(sorted_tips, sorted_tips[1])
sorted_angles <- c(sorted_angles, sorted_angles[1] + 2*pi)

# Define ring positions - only two rings
ring_gap <- 0.12 * max_radius
ring_width <- 0.12 * max_radius
ring_separation <- 0.06 * max_radius

phylum_inner <- max_radius + ring_gap
phylum_outer <- phylum_inner + ring_width

genus_inner <- phylum_outer + ring_separation
genus_outer <- genus_inner + ring_width

# Draw segments
for (i in 1:(length(sorted_tips)-1)) {{
  tryCatch({{
    tip <- sorted_tips[i]
    start_angle <- sorted_angles[i]
    end_angle <- sorted_angles[i+1]
    
    # Phylum ring
    phylum <- get_taxon(tip, tip_to_phylum)
    phylum_color <- get_color(phylum, phylum_colors)
    draw_arc_segment(center_x, center_y, phylum_inner, phylum_outer, 
                    start_angle, end_angle, phylum_color)
    
    # Genus ring
    genus <- get_taxon(tip, tip_to_genus)
    genus_color <- get_color(genus, genus_colors)
    draw_arc_segment(center_x, center_y, genus_inner, genus_outer, 
                    start_angle, end_angle, genus_color)
    
  }}, error = function(e) {{
    cat("Error drawing segment", i, ":", e$message, "\\n")
  }})
}}

# Add tip points
for (i in 1:length(tree$tip.label)) {{
  tip <- tree$tip.label[i]
  phylum <- get_taxon(tip, tip_to_phylum)
  phylum_color <- get_color(phylum, phylum_colors)
  
  # Base point colored by phylum
  points(tip_coords_x[i], tip_coords_y[i], pch=19, col=phylum_color, cex=1.8)
  points(tip_coords_x[i], tip_coords_y[i], pch=1, col="#333333", cex=1.8, lwd=1.2)
  
  # Highlight novel MAGs
  if (get_taxon(tip, tip_to_novel)) {{
    points(tip_coords_x[i], tip_coords_y[i], pch=1, col="#FF0000", cex=4, lwd=2.5)
    points(tip_coords_x[i], tip_coords_y[i], pch=8, col="#FF0000", cex=3, lwd=2)
  }}
}}

# Add connecting lines
for (i in 1:length(sorted_tips) - 1) {{
  tip_idx <- which(tree$tip.label == sorted_tips[i])
  if (length(tip_idx) > 0) {{
    angle <- sorted_angles[i]
    x1 <- tip_coords_x[tip_idx]
    y1 <- tip_coords_y[tip_idx]
    x2 <- center_x + phylum_inner * cos(angle)
    y2 <- center_y + phylum_inner * sin(angle)
    
    lines(c(x1, x2), c(y1, y2), col="#E0E0E0", lwd=1, lty=1)
  }}
}}

# Legend
par(mar=c(2,1,3,1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Title
text(0.5, 0.98, "Legend", font=2, cex=1.8, adj=0.5)

# Phylum legend
text(0.05, 0.90, "Phylum", font=2, pos=4, cex=1.5)
legend_y <- 0.84
for (i in 1:min(12, length(phyla))) {{
  rect(0.05, legend_y - 0.018, 0.13, legend_y + 0.018, 
       col=phylum_colors[phyla[i]], border="black", lwd=1)
  text(0.15, legend_y, phyla[i], pos=4, cex=1.0)
  legend_y <- legend_y - 0.038
}}

# Genus legend
legend_y <- legend_y - 0.05
text(0.05, legend_y, "Genus", font=2, pos=4, cex=1.5)
legend_y <- legend_y - 0.04
for (i in 1:min(8, length(genera))) {{
  rect(0.05, legend_y - 0.018, 0.13, legend_y + 0.018, 
       col=genus_colors[genera[i]], border="black", lwd=1)
  text(0.15, legend_y, genera[i], pos=4, cex=1.0)
  legend_y <- legend_y - 0.038
}}

# Novel MAG indicator
legend_y <- legend_y - 0.05
text(0.05, legend_y, "Novel MAGs", font=2, pos=4, cex=1.5)
legend_y <- legend_y - 0.04
rect(0.05, legend_y - 0.018, 0.13, legend_y + 0.018, 
     col="#FF0000", border="black", lwd=1)
text(0.15, legend_y, paste("Novel (n=", novel_count, ")", sep=""), pos=4, cex=1.0)

dev.off()

cat("\\n=== Simplified Taxonomic Visualization Complete ===\\n")
cat("Created:", simplified_file, "\\n")
cat(sprintf("Summary: %d total MAGs, %d novel MAGs\\n", nrow(tip_data), novel_count))

quit(status=0)
'''
        
        # Write and execute R script
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        os.chmod(script_path, 0o755)
        
        # Run the R script
        success = run_r_command_with_conda(script_path)
        
        return success
        
    except Exception as e:
        logger.error(f"Error creating simplified taxonomic visualization: {e}")
        return False


def create_rectangular_tree_visualization(tree_file, metadata_file, output_prefix, novel_mags):
    """Create rectangular tree with colored strips for phylum and genus."""
    try:
        script_path = f"{output_prefix}_rectangular_vis.R"
        
        logger.info(f"Creating rectangular tree visualization")
        logger.info(f"  - {len(novel_mags)} novel MAG entries")
        
        # Prepare data for R script
        novel_mags_json = json.dumps(novel_mags)
        
        r_script = f'''#!/usr/bin/env Rscript

# Rectangular Tree Visualization with Taxonomic Strips

# Load required packages
tryCatch({{
  library(ape)
  library(RColorBrewer)
  library(jsonlite)
  cat("Base packages loaded successfully\\n")
}}, error = function(e) {{
  cat("Error loading packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- fromJSON('{novel_mags_json}')

set.seed(42)

cat("=== Rectangular Tree Visualization ===\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")

# Read and process tree
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded successfully with", length(tree$tip.label), "tips\\n")
  
  # Ladderize the tree for better visualization
  tree <- ladderize(tree, right = FALSE)
  
  if (!is.rooted(tree)) {{
    if (length(tree$tip.label) > 2) {{
      outgroup_tip <- tree$tip.label[1]
      tryCatch({{
        tree <- root(tree, outgroup = outgroup_tip)
        tree <- ladderize(tree, right = FALSE)
        cat("Applied root using tip:", outgroup_tip, "\\n")
      }}, error = function(e) {{
        cat("Rooting failed, using tree as-is\\n")
      }})
    }}
  }}
}}, error = function(e) {{
  cat("Error loading tree:", e$message, "\\n")
  quit(status=1)
}})

# Read metadata
tryCatch({{
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns\\n")
  id_col <- colnames(metadata)[1]
}}, error = function(e) {{
  cat("Error loading metadata:", e$message, "\\n")
  quit(status=1)
}})

# Process taxonomy data
taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]])

for (level in taxonomic_levels) {{
  if (level %in% colnames(metadata)) {{
    taxa_data[[level]] <- as.character(metadata[[level]])
    taxa_data[[level]][taxa_data[[level]] %in% c("Unknown", "", "NA", "na", "None")] <- "Unknown"
  }} else {{
    taxa_data[[level]] <- "Unknown"
  }}
}}

# Parse classification column if exists
if ("classification" %in% colnames(metadata)) {{
  cat("Extracting taxonomy from classification column\\n")
  for (i in 1:nrow(metadata)) {{
    if (!is.na(metadata$classification[i])) {{
      taxa <- strsplit(metadata$classification[i], ";")[[1]]
      for (taxon in taxa) {{
        if (startsWith(taxon, "d__")) {{
          taxa_data$Domain[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "p__")) {{
          taxa_data$Phylum[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "c__")) {{
          taxa_data$Class[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "o__")) {{
          taxa_data$Order[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "f__")) {{
          taxa_data$Family[i] <- substring(taxon, 4)
        }} else if (startsWith(taxon, "g__")) {{
          taxa_data$Genus[i] <- substring(taxon, 4)  
        }} else if (startsWith(taxon, "s__")) {{
          taxa_data$Species[i] <- substring(taxon, 4)
        }}
      }}
    }}
  }}
}}

# Match tree tips with taxa data
match_tip_to_taxa <- function(tip, taxa_data) {{
  match_idx <- which(taxa_data$ID == tip)
  if (length(match_idx) == 0) {{
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  if (length(match_idx) == 0) {{
    converted_tip <- gsub("\\\\.", "_", tip)
    match_idx <- which(taxa_data$ID == converted_tip)
  }}
  if (length(match_idx) == 0) {{
    for (i in 1:nrow(taxa_data)) {{
      base_taxa <- gsub("_sub$", "", taxa_data$ID[i])
      base_tip <- gsub("_sub$", "", tip)
      if (base_taxa == base_tip || 
          grepl(base_taxa, base_tip, fixed=TRUE) || 
          grepl(base_tip, base_taxa, fixed=TRUE)) {{
        match_idx <- i
        break
      }}
    }}
  }}
  
  if (length(match_idx) > 0) {{
    return(taxa_data[match_idx[1], ])
  }} else {{
    unknown_row <- taxa_data[1, ]
    unknown_row[] <- c(tip, rep("Unknown", length(taxonomic_levels)))
    return(unknown_row)
  }}
}}

# Create ordered tip data based on tree tip order
tip_data <- do.call(rbind, lapply(tree$tip.label, function(tip) {{
  return(match_tip_to_taxa(tip, taxa_data))
}}))

# Add novel MAG information
tip_data$is_novel <- FALSE

for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  if (tip %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    next
  }}
  
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    next
  }}
  
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
  }}
}}

novel_count <- sum(tip_data$is_novel)
cat("Novel MAGs matched:", novel_count, "\\n")

# Create color palettes
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

# Use distinct color palettes
if (length(phyla) > 0) {{
  if (length(phyla) <= 8) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Set1")[1:length(phyla)]
  }} else if (length(phyla) <= 12) {{
    phylum_colors <- c(brewer.pal(8, "Set1"), brewer.pal(4, "Dark2")[1:(length(phyla)-8)])
  }} else {{
    phylum_colors <- colorRampPalette(c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2")))(length(phyla))
  }}
  names(phylum_colors) <- phyla
  phylum_colors <- c(phylum_colors, "Unknown" = "#F5F5F5")
}} else {{
  phylum_colors <- c("Unknown" = "#F5F5F5")
}}

if (length(genera) > 0) {{
  if (length(genera) <= 9) {{
    genus_colors <- brewer.pal(max(3, length(genera)), "Pastel1")[1:length(genera)]
  }} else {{
    genus_colors <- colorRampPalette(brewer.pal(9, "Pastel1"))(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "#FAFAFA")
}} else {{
  genus_colors <- c("Unknown" = "#FAFAFA")
}}

##############################################
# CREATE RECTANGULAR TREE VISUALIZATION
##############################################

rectangular_file <- paste0(output_prefix, "_rectangular_tree.pdf")
pdf(rectangular_file, width = 18, height = max(12, length(tree$tip.label) * 0.25))

# Layout for tree and legend
layout(matrix(c(1,2), nrow=1), widths=c(4,1))

# Main tree plot
par(mar=c(3,1,3,10), xpd=TRUE)

# Plot rectangular tree without tip labels
plot(tree, type="phylogram", show.tip.label=FALSE, 
     main="Rectangular Phylogenetic Tree with Taxonomic Classification", 
     edge.width=2, cex.main=1.6, edge.color="#2C3E50",
     x.lim=c(0, max(node.depth.edgelength(tree)) * 1.3))

# Get tree layout information
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
n_tips <- length(tree$tip.label)

# Calculate positions for colored strips
x_tree_max <- max(last_plot$xx[1:n_tips])
strip_start_x <- x_tree_max + 0.02 * max(last_plot$xx)
strip_width <- 0.025 * max(last_plot$xx)
strip_gap <- 0.005 * max(last_plot$xx)

# Draw colored strips for each tip
for (i in 1:n_tips) {{
  tip <- tree$tip.label[i]
  y_pos <- last_plot$yy[i]
  
  # Get taxonomy for this tip
  tip_idx <- which(tip_data$ID == tip)
  if (length(tip_idx) > 0) {{
    phylum <- tip_data$Phylum[tip_idx[1]]
    genus <- tip_data$Genus[tip_idx[1]]
    is_novel <- tip_data$is_novel[tip_idx[1]]
    
    # Phylum strip
    phylum_color <- if (phylum != "Unknown" && phylum %in% names(phylum_colors)) {{
      phylum_colors[phylum]
    }} else {{
      phylum_colors["Unknown"]
    }}
    
    rect(strip_start_x, y_pos - 0.4, 
         strip_start_x + strip_width, y_pos + 0.4,
         col = phylum_color, border = NA)
    
    # Genus strip
    genus_color <- if (genus != "Unknown" && genus %in% names(genus_colors)) {{
      genus_colors[genus]
    }} else {{
      genus_colors["Unknown"]
    }}
    
    rect(strip_start_x + strip_width + strip_gap, y_pos - 0.4, 
         strip_start_x + 2*strip_width + strip_gap, y_pos + 0.4,
         col = genus_color, border = NA)
    
    # Novel MAG indicator
    if (is_novel) {{
      # Add red border around both strips
      rect(strip_start_x - 0.001 * max(last_plot$xx), y_pos - 0.45, 
           strip_start_x + 2*strip_width + strip_gap + 0.001 * max(last_plot$xx), y_pos + 0.45,
           col = NA, border = "#FF0000", lwd = 2)
      
      # Add star symbol
      points(strip_start_x + 2*strip_width + strip_gap + 0.015 * max(last_plot$xx), 
             y_pos, pch = 8, col = "#FF0000", cex = 1.5, lwd = 2)
    }}
  }}
}}

# Add strip labels at the top
text(strip_start_x + strip_width/2, n_tips + 1.5, "Phylum", 
     font = 2, cex = 1.2, srt = 0)
text(strip_start_x + strip_width + strip_gap + strip_width/2, n_tips + 1.5, "Genus", 
     font = 2, cex = 1.2, srt = 0)

# Legend
par(mar=c(5,2,5,2), xpd=TRUE)
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Title
text(0.5, 0.98, "Legend", font=2, cex=1.6, adj=0.5)

# Phylum colors
text(0.05, 0.92, "Phylum", font=2, pos=4, cex=1.3)
legend_y <- 0.86
n_phyla_shown <- min(15, length(phyla))
for (i in 1:n_phyla_shown) {{
  rect(0.05, legend_y - 0.015, 0.15, legend_y + 0.015, 
       col=phylum_colors[phyla[i]], border="black", lwd=0.8)
  # Truncate long names
  phylum_name <- phyla[i]
  if (nchar(phylum_name) > 20) {{
    phylum_name <- paste0(substr(phylum_name, 1, 18), "...")
  }}
  text(0.17, legend_y, phylum_name, pos=4, cex=0.9)
  legend_y <- legend_y - 0.032
}}

if (length(phyla) > n_phyla_shown) {{
  text(0.17, legend_y, paste0("... and ", length(phyla) - n_phyla_shown, " more"), 
       pos=4, cex=0.8, font=3)
  legend_y <- legend_y - 0.032
}}

# Genus colors
legend_y <- legend_y - 0.04
text(0.05, legend_y, "Genus", font=2, pos=4, cex=1.3)
legend_y <- legend_y - 0.035
n_genera_shown <- min(10, length(genera))
for (i in 1:n_genera_shown) {{
  rect(0.05, legend_y - 0.015, 0.15, legend_y + 0.015, 
       col=genus_colors[genera[i]], border="black", lwd=0.8)
  # Truncate long names
  genus_name <- genera[i]
  if (nchar(genus_name) > 20) {{
    genus_name <- paste0(substr(genus_name, 1, 18), "...")
  }}
  text(0.17, legend_y, genus_name, pos=4, cex=0.9)
  legend_y <- legend_y - 0.032
}}

if (length(genera) > n_genera_shown) {{
  text(0.17, legend_y, paste0("... and ", length(genera) - n_genera_shown, " more"), 
       pos=4, cex=0.8, font=3)
  legend_y <- legend_y - 0.032
}}

# Novel MAG indicator
legend_y <- legend_y - 0.04
text(0.05, legend_y, "Novel MAGs", font=2, pos=4, cex=1.3)
legend_y <- legend_y - 0.035
rect(0.05, legend_y - 0.015, 0.15, legend_y + 0.015, 
     col=NA, border="#FF0000", lwd=2)
points(0.10, legend_y, pch=8, col="#FF0000", cex=1.5, lwd=2)
text(0.17, legend_y, paste0("Novel (n=", novel_count, ")"), pos=4, cex=0.9)

# Statistics
legend_y <- legend_y - 0.06
text(0.05, legend_y, "Statistics", font=2, pos=4, cex=1.3)
legend_y <- legend_y - 0.035
text(0.05, legend_y, paste0("Total MAGs: ", n_tips), pos=4, cex=0.9)
legend_y <- legend_y - 0.025
text(0.05, legend_y, paste0("Phyla: ", length(phyla)), pos=4, cex=0.9)
legend_y <- legend_y - 0.025
text(0.05, legend_y, paste0("Genera: ", length(genera)), pos=4, cex=0.9)

dev.off()

cat("\\n=== Rectangular Tree Visualization Complete ===\\n")
cat("Created:", rectangular_file, "\\n")
cat(sprintf("Summary: %d total MAGs, %d novel MAGs, %d phyla, %d genera\\n", 
           n_tips, novel_count, length(phyla), length(genera)))

quit(status=0)
'''
        
        # Write and execute R script
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        os.chmod(script_path, 0o755)
        
        # Run the R script
        success = run_r_command_with_conda(script_path)
        
        return success
        
    except Exception as e:
        logger.error(f"Error creating rectangular tree visualization: {e}")
        return False

def create_multi_functional_visualization(tree_file, metadata_file, output_prefix, 
                                        novel_mags, functional_data):
    """Create multi-functional visualization with CAZyme, COG, and KEGG pie charts."""
    if PIE_MODULE_AVAILABLE:
        logger.info("Delegating multi-functional visualization to pie module...")
        return pie_viz.create_multi_functional_visualization(
            tree_file, metadata_file, output_prefix, novel_mags, functional_data
        )
    else:
        logger.warning("Pie module not available, skipping multi-functional visualization")
        return False

def create_publication_quality_visualization(tree_file, metadata_file, output_prefix, 
                                           novel_mags, functional_data, annotation_type):
    """Create publication-quality visualization with single functional annotation."""
    if PIE_MODULE_AVAILABLE:
        logger.info("Delegating publication quality visualization to pie module...")
        return pie_viz.create_publication_quality_visualization(
            tree_file, metadata_file, output_prefix, novel_mags, functional_data, annotation_type
        )
    else:
        logger.warning("Pie module not available, skipping publication quality visualization")
        return False

if __name__ == "__main__":
    import sys
    import os
    import shutil
    
    if len(sys.argv) < 6:
        print("Usage: python enhanced_phylo_tree_vis.py <cpus> <memory> <time_limit> <project_input> <project_output> [taxonomy_source] [novel_mags_file] [functional_annotations] [annotation_type] [cazyme_annotations] [cog_annotations] [kegg_annotations]")
        sys.exit(1)
    
    cpus = int(sys.argv[1])
    memory = sys.argv[2]
    time_limit = sys.argv[3]
    project_input = sys.argv[4]
    project_output = sys.argv[5]  # Save the original
    taxonomy_source = sys.argv[6] if len(sys.argv) > 6 else None
    novel_mags_file = sys.argv[7] if len(sys.argv) > 7 else None
    functional_annotations = sys.argv[8] if len(sys.argv) > 8 else None
    annotation_type = sys.argv[9] if len(sys.argv) > 9 else "auto"
    cazyme_annotations = sys.argv[10] if len(sys.argv) > 10 else None
    cog_annotations = sys.argv[11] if len(sys.argv) > 11 else None
    kegg_annotations = sys.argv[12] if len(sys.argv) > 12 else None
    
    # Save the base output directory
    base_output = project_output
    
    # PHASE 1: Run with modified path that will create enhanced_visualizations_finalrec
    # Create a fake project_output that when combined with hardcoded path gives us what we want
    # The function adds: "Phylogeny/gtdbtk/taxonomy2/enhanced_visualizations"
    # So we create a temp directory structure
    
    import tempfile
    temp_base = tempfile.mkdtemp()
    
    # Create the expected structure for first run
    fake_output_1 = os.path.join(temp_base, "run1")
    os.makedirs(fake_output_1, exist_ok=True)
    
    logger.info("=" * 60)
    logger.info("PHASE 1: Running simplified and rectangular visualizations...")
    logger.info("=" * 60)
    
    success1 = run_enhanced_tree_visualization(
        cpus, memory, time_limit, project_input, fake_output_1,
        taxonomy_source, novel_mags_file, functional_annotations, annotation_type,
        cazyme_annotations, cog_annotations, kegg_annotations
    )
    
    # Move the results to the real location with new name
    source_dir_1 = os.path.join(fake_output_1, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations")
    dest_dir_1 = os.path.join(base_output, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations_finalrec")
    
    if os.path.exists(source_dir_1):
        os.makedirs(os.path.dirname(dest_dir_1), exist_ok=True)
        if os.path.exists(dest_dir_1):
            shutil.rmtree(dest_dir_1)
        shutil.move(source_dir_1, dest_dir_1)
        logger.info(f"Moved visualizations to: {dest_dir_1}")
    
    # PHASE 2: Run pie module
    if PIE_MODULE_AVAILABLE:
        fake_output_2 = os.path.join(temp_base, "run2")
        os.makedirs(fake_output_2, exist_ok=True)
        
        logger.info("")
        logger.info("=" * 60)
        logger.info("PHASE 2: Running pie chart visualizations...")
        logger.info("=" * 60)
        
        success2 = pie_viz.run_enhanced_tree_visualization(
            cpus, memory, time_limit, project_input, fake_output_2,
            taxonomy_source, novel_mags_file, functional_annotations, annotation_type,
            cazyme_annotations, cog_annotations, kegg_annotations
        )
        
        # Move the results to the real location with new name
        source_dir_2 = os.path.join(fake_output_2, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations")
        dest_dir_2 = os.path.join(base_output, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations_pie")
        
        if os.path.exists(source_dir_2):
            os.makedirs(os.path.dirname(dest_dir_2), exist_ok=True)
            if os.path.exists(dest_dir_2):
                shutil.rmtree(dest_dir_2)
            shutil.move(source_dir_2, dest_dir_2)
            logger.info(f"Moved pie visualizations to: {dest_dir_2}")
    else:
        logger.warning("Pie module not available")
        success2 = False
    
    # Clean up temp directory
    shutil.rmtree(temp_base, ignore_errors=True)
    
    # List all files in both directories
    logger.info("")
    logger.info("=" * 60)
    logger.info("ALL CREATED VISUALIZATION FILES:")
    logger.info("=" * 60)
    
    finalrec_dir = os.path.join(base_output, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations_finalrec")
    if os.path.exists(finalrec_dir):
        files = [f for f in os.listdir(finalrec_dir) if f.endswith('.pdf')]
        if files:
            logger.info(f"\nFinal_rec visualizations ({finalrec_dir}):")
            for file in sorted(files):
                file_size = os.path.getsize(os.path.join(finalrec_dir, file)) / 1024
                logger.info(f"   {file} ({file_size:.1f} KB)")
    
    pie_dir = os.path.join(base_output, "Phylogeny", "gtdbtk", "taxonomy2", "enhanced_visualizations_pie")
    if os.path.exists(pie_dir):
        files = [f for f in os.listdir(pie_dir) if f.endswith('.pdf')]
        if files:
            logger.info(f"\nPie visualizations ({pie_dir}):")
            for file in sorted(files):
                file_size = os.path.getsize(os.path.join(pie_dir, file)) / 1024
                logger.info(f"   {file} ({file_size:.1f} KB)")
    
    # Final summary
    logger.info("")
    logger.info("=" * 60)
    if PIE_MODULE_AVAILABLE:
        if success1 and success2:
            logger.info(" Both visualization modules completed successfully")
            sys.exit(0)
        elif success1:
            logger.info(" Only simplified/rectangular visualizations completed")
            sys.exit(1)
        else:
            logger.info(" Errors occurred in visualization")
            sys.exit(1)
    else:
        if success1:
            logger.info(" Simplified/rectangular visualizations completed (pie module not available)")
            sys.exit(0)
        else:
            logger.info(" Visualization failed")
            sys.exit(1)
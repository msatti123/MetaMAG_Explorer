"""
Simple Phylogenetic Tree Visualization with Base R

This module provides functions to create publication-quality tree visualizations 
using only base R packages (no ggtree or Bioconductor required).
"""

import os
import subprocess
import pandas as pd
import numpy as np
import re
import time
import shutil
import logging
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import setup_r_environment as config_setup_r, get_r_config, detect_r_environment

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_r_command_with_conda(script_path):
    """Run R script with proper conda environment activation."""
    try:
        from MetaMAG.config import get_r_config
        r_config = get_r_config()
        
        if not r_config:
            logger.error("No R configuration found")
            return False
        
        if r_config['env_name'] != 'system' and r_config.get('activate_cmd'):
            cmd = f"{r_config['activate_cmd']} && Rscript {script_path}"
            logger.info(f"Using conda R environment: {r_config['env_name']}")
        else:
            cmd = f"Rscript {script_path}"
            logger.info("Using system R")
        
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

def run_tree_visualization(cpus, memory, time_limit, project_input, project_output, taxonomy_source=None):
    """
    Runs the enhanced phylogenetic tree visualization on the generated tree files.
    
    Parameters:
    - cpus: Number of CPUs to use
    - memory: Memory allocation
    - time_limit: Time limit
    - project_input: Input directory path
    - project_output: Output directory path
    - taxonomy_source: Path to GTDB-Tk taxonomy file (optional, will search if not provided)
    
    Returns:
    - Boolean indicating success or failure
    """
    start_time = time.time()
    logger.info("Starting enhanced phylogenetic tree visualization...")
    
    # Set up custom R environment
    setup_r_environment()
    
    # Define paths
    phylogeny_dir = os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy2")
    trees_dir = os.path.join(phylogeny_dir, "trees")
    vis_output_dir = os.path.join(phylogeny_dir, "visualizations")
    
    # Ensure output directory exists
    ensure_directory_exists(vis_output_dir)
    
    # Find tree files
    tree_files = find_tree_files(trees_dir, project_output)
    
    if not tree_files:
        logger.error(f"No tree files found in {trees_dir} or {project_output}")
        return False
    
    logger.info(f"Found {len(tree_files)} tree file(s) for visualization")
    
    # Find or create taxonomy metadata
    if taxonomy_source and os.path.exists(taxonomy_source):
        logger.info(f"Using provided taxonomy source: {taxonomy_source}")
        metadata_file = process_metadata_file(taxonomy_source, os.path.join(vis_output_dir, "Metadata_For_Tree.csv"))
    else:
        # Search for taxonomy files in common locations
        metadata_file = find_taxonomy_files(project_output, vis_output_dir)
    
    if not metadata_file:
        logger.warning("No metadata file found. Creating from tree file...")
        metadata_file = create_metadata_from_tree(tree_files[0], vis_output_dir)
    
    # Process each tree file
    all_success = True
    for tree_file in tree_files:
        tree_basename = os.path.basename(tree_file)
        output_prefix = os.path.join(vis_output_dir, os.path.splitext(tree_basename)[0])
        
        # Create visualization
        success = create_publication_quality_visualization(tree_file, metadata_file, output_prefix)
        
        if not success:
            logger.warning(f"Failed to create visualization for {tree_basename}")
            all_success = False
        else:
            logger.info(f"Successfully created visualization for {tree_basename}")
    
    end_time = time.time()
    logger.info(f"Tree visualization completed in {end_time - start_time:.2f} seconds")
    return all_success

def setup_r_environment():
    """
    Set up the custom R environment using portable config detection.
    """
    from MetaMAG.config import setup_r_environment as config_setup_r, get_r_config
    
    logger.info("Setting up R environment using portable detection...")
    
    # Use the config system to set up R environment
    success = config_setup_r()
    
    if success:
        r_config = get_r_config()
        if r_config:
            logger.info(f"R environment configured: {r_config['env_name']}")
            if r_config.get('r_home'):
                logger.info(f"R_HOME: {r_config['r_home']}")
            if r_config.get('r_libs'):
                logger.info(f"R_LIBS_USER: {r_config['r_libs']}")
            return True
        else:
            logger.error("R configuration not available")
            return False
    else:
        logger.error("Failed to configure R environment")
        logger.error("Please ensure you have:")
        logger.error("  1. Conda installed")
        logger.error("  2. An R environment created (conda create -n R r-base)")
        logger.error("  3. Or system R installed")
        return False


def find_tree_files(trees_dir, project_output):
    """
    Find all phylogenetic tree files in the project directories.
    """
    tree_files = []
    
    # Search in the trees_dir first
    if os.path.exists(trees_dir):
        logger.info(f"Searching for tree files in {trees_dir}")
        for file in os.listdir(trees_dir):
            if file.endswith(".tree") or file.endswith(".nwk") or "tree" in file.lower():
                # Skip invalid tree files like .tree-table or .tree-taxonomy
                if "-" not in file and ".log" not in file:
                    tree_files.append(os.path.join(trees_dir, file))
    
    # If no tree files found, search in alternative locations
    if not tree_files:
        logger.info("No tree files found in primary directory, searching alternative locations")
        
        # Common locations for tree files
        alternative_locations = [
            os.path.join(project_output, "Phylogeny", "gtdbtk"),
            os.path.join(project_output, "Phylogeny"),
            os.path.join(project_output, "Novel_Mags", "gtdbtk"),
            os.path.join(project_output, "gtdbtk"),
            project_output
        ]
        
        for location in alternative_locations:
            if os.path.exists(location):
                for file in os.listdir(location):
                    if (file.endswith(".tree") or file.endswith(".nwk") or "tree" in file.lower()) and "-" not in file and ".log" not in file:
                        tree_files.append(os.path.join(location, file))
                
                # If we found tree files, no need to check other locations
                if tree_files:
                    logger.info(f"Found tree files in {location}")
                    break
    
    return tree_files

def find_taxonomy_files(project_output, output_dir):
    """
    Find GTDB-Tk taxonomy files in the project output.
    """
    # Common locations for taxonomy files
    potential_locations = [
        os.path.join(project_output, "Novel_Mags", "gtdbtk"),
        os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy"),
        os.path.join(project_output, "gtdbtk"),
        project_output
    ]
    
    # Search patterns for taxonomy files
    taxonomy_patterns = [
        "gtdbtk.bac120.summary.tsv",
        "gtdbtk.ar122.summary.tsv",
        "*summary.tsv",
        "*.classifications.txt"
    ]
    
    # Search for taxonomy files
    for location in potential_locations:
        if os.path.exists(location):
            for pattern in taxonomy_patterns:
                pattern_path = os.path.join(location, pattern)
                import glob
                matching_files = glob.glob(pattern_path)
                
                if matching_files:
                    taxonomy_file = matching_files[0]
                    logger.info(f"Found taxonomy file: {taxonomy_file}")
                    return process_metadata_file(taxonomy_file, os.path.join(output_dir, "Metadata_For_Tree.csv"))
    
    logger.warning("No taxonomy files found in common locations")
    return None

def process_metadata_file(taxonomy_file, output_path):
    """
    Process a GTDB-Tk taxonomy file into metadata suitable for tree visualization.
    """
    try:
        # Determine file type and read accordingly
        if taxonomy_file.endswith('.tsv'):
            df = pd.read_csv(taxonomy_file, sep='\t')
        else:
            # Try different separators
            for sep in ['\t', ',']:
                try:
                    df = pd.read_csv(taxonomy_file, sep=sep)
                    break
                except:
                    continue
            else:
                logger.error(f"Could not determine file format for {taxonomy_file}")
                return None
        
        logger.info(f"Read taxonomy file with {df.shape[0]} rows and {df.shape[1]} columns")
        
        # Process the dataframe to extract taxonomy information
        metadata_df = pd.DataFrame()
        
        # Get the ID column (user_genome, bin_id, etc.)
        id_column = None
        for potential_id in ['user_genome', 'genome', 'bin_id', 'name', 'Genome']:
            if potential_id in df.columns:
                id_column = potential_id
                break
        
        if not id_column:
            # Use first column as ID
            id_column = df.columns[0]
            logger.warning(f"No standard ID column found, using {id_column}")
        
        metadata_df['user_genome'] = df[id_column]
        
        # Extract taxonomy information
        if 'classification' in df.columns:
            # Extract from GTDB-style classification string
            logger.info("Extracting taxonomy from classification column")
            
            # Create columns for different taxonomic levels
            taxonomic_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            for level in taxonomic_levels:
                metadata_df[level] = 'Unknown'
            
            # Extract taxonomy information
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
        else:
            # Try to find direct taxonomic columns
            tax_map = {
                'phylum': 'Phylum',
                'genus': 'Genus',
                'family': 'Family',
                'domain': 'Domain',
                'class': 'Class',
                'order': 'Order',
                'species': 'Species'
            }
            
            # Check for each taxonomic level
            for orig, new in tax_map.items():
                # Find columns that match this taxonomic level
                cols = [col for col in df.columns if col.lower() == orig]
                if cols:
                    metadata_df[new] = df[cols[0]]
                else:
                    metadata_df[new] = 'Unknown'
            
        # Clean up values
        for col in metadata_df.columns:
            if col != 'user_genome':
                # Remove prefix from taxonomic names
                metadata_df[col] = metadata_df[col].astype(str).apply(
                    lambda x: x[3:] if isinstance(x, str) and len(x) > 3 and x[1:3] == '__' else x
                )
                
                # Replace empty or null values
                metadata_df[col] = metadata_df[col].replace(['', 'nan', 'None', 'none'], 'Unknown')
                metadata_df[col] = metadata_df[col].fillna('Unknown')
        
        # Save processed metadata
        metadata_df.to_csv(output_path, index=False)
        logger.info(f"Processed metadata saved to {output_path}")
        
        return output_path
    
    except Exception as e:
        logger.error(f"Error processing metadata file: {e}")
        return None

def create_metadata_from_tree(tree_file, output_dir):
    """
    Create a basic metadata file from the tree file.
    """
    output_path = os.path.join(output_dir, "Metadata_From_Tree.csv")
    
    try:
        # Read the tree file
        with open(tree_file, 'r') as f:
            tree_content = f.read()
        
        # Extract genome IDs using a more robust pattern
        # This pattern looks for strings that are not tree formatting characters
        genome_ids = re.findall(r'([a-zA-Z0-9_.]+)', tree_content)
        
        # Filter out branch lengths and other numerical values
        genome_ids = [gid for gid in genome_ids if not re.match(r'^\d+(\.\d+)?$', gid)]
        
        # Remove duplicate IDs
        unique_ids = []
        for gid in genome_ids:
            if gid not in unique_ids and len(gid) > 1:  # Avoid single-character IDs
                unique_ids.append(gid)
        
        logger.info(f"Extracted {len(unique_ids)} unique genome IDs from tree")
        
        # Create basic metadata with placeholder taxonomic information
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
        
        # Create dataframe and save
        metadata_df = pd.DataFrame(data)
        metadata_df.to_csv(output_path, index=False)
        logger.info(f"Created basic metadata from tree file: {output_path}")
        
        return output_path
    
    except Exception as e:
        logger.error(f"Error creating metadata from tree: {e}")
        return None

def create_publication_quality_visualization(tree_file, metadata_file, output_prefix):
    """
    Create publication-quality phylogenetic tree visualizations using only the most basic ggtree functions.
    No external dependencies beyond core ggtree, ape, ggplot2, and RColorBrewer.
    
    Parameters:
    - tree_file: Path to the Newick tree file
    - metadata_file: Path to the metadata CSV file
    - output_prefix: Output directory and prefix for saving files
    
    Returns:
    - Boolean indicating success or failure
    """
    try:
        # Create an R script using only the most basic ggtree functions
        script_path = f"{output_prefix}_minimal_ggtree.R"
        
        # Write R script content using absolutely minimal dependencies
        r_script = f"""#!/usr/bin/env Rscript

# Enable detailed error reporting
options(error = function() {{ 
  traceback(20)
  if (!interactive()) {{ quit("no", status = 1, runLast = FALSE) }}
}})

# Load only absolutely necessary packages
tryCatch({{
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(RColorBrewer)
  cat("Base packages loaded successfully\\n")
}}, error = function(e) {{
  cat("Error loading packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"

# Read tree
cat("Reading tree from", tree_file, "\\n")
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded with", length(tree$tip.label), "tips\\n")
}}, error = function(e) {{
  cat("Error reading tree file:", e$message, "\\n")
  quit(status = 1)
}})

# Read metadata
cat("Reading metadata from", metadata_file, "\\n")
tryCatch({{
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns\\n")
}}, error = function(e) {{
  cat("Error reading metadata file:", e$message, "\\n")
  quit(status = 1)
}})

# Print sample of the metadata to debug
cat("First few rows of metadata:\\n")
print(head(metadata))

# Identify ID column
id_col <- colnames(metadata)[1]
cat("Using", id_col, "as the ID column\\n")

# Handle missing values for visualization
if ("Phylum" %in% colnames(metadata)) {{
  metadata$Phylum <- as.character(metadata$Phylum)
  cat("Phylum column found\\n")
}} else {{
  metadata$Phylum <- NA
  cat("No Phylum column found, adding empty column\\n")
}}

if ("Genus" %in% colnames(metadata)) {{
  metadata$Genus <- as.character(metadata$Genus)
  cat("Genus column found\\n")
}} else {{
  metadata$Genus <- NA
  cat("No Genus column found, adding empty column\\n")
}}

# Replace 'Unknown' with NA for consistent handling
metadata$Phylum[metadata$Phylum == "Unknown" | metadata$Phylum == ""] <- NA
metadata$Genus[metadata$Genus == "Unknown" | metadata$Genus == ""] <- NA

# Make sure tree and metadata match
tree_tips <- tree$tip.label
metadata_ids <- metadata[[id_col]]
matching_ids <- intersect(tree_tips, metadata_ids)

cat("Tree has", length(tree_tips), "tips\\n")
cat("Metadata has", length(metadata_ids), "entries\\n")
cat("Found", length(matching_ids), "matching IDs\\n")

if (length(matching_ids) == 0) {{
  cat("WARNING: No matching IDs found between tree and metadata!\\n")
  cat("This may result in a tree without proper coloring.\\n")
  
  # Create a dummy metadata dataframe mapping tree tips to default values
  metadata <- data.frame(
    id = tree_tips,
    Phylum = rep(NA, length(tree_tips)),
    Genus = rep(NA, length(tree_tips)),
    stringsAsFactors = FALSE
  )
  colnames(metadata)[1] <- id_col
  cat("Created default metadata for tree tips\\n")
}}

# Prepare colors for visualization
# Get unique non-NA values for Phylum and Genus
phyla <- unique(metadata$Phylum[!is.na(metadata$Phylum)])
genera <- unique(metadata$Genus[!is.na(metadata$Genus)])

cat("Found", length(phyla), "unique phyla:\\n")
print(phyla)
cat("Found", length(genera), "unique genera:\\n")
print(head(genera))

# Create color palettes
# For Phylum - use Set1 which has bright distinguishable colors
phylum_colors <- brewer.pal(min(9, max(3, length(phyla))), "Set1")
if (length(phyla) > 9) {{
  phylum_colors <- colorRampPalette(phylum_colors)(length(phyla))
}}
names(phylum_colors) <- phyla

# For Genus - use Set3 for a different color scheme (more colors than Set2)
genus_colors <- brewer.pal(min(12, max(3, length(genera))), "Set3")
if (length(genera) > 12) {{
  genus_colors <- colorRampPalette(genus_colors)(length(genera))
}}
names(genus_colors) <- genera

# Print color assignments for debugging
cat("Phylum color assignments:\\n")
print(phylum_colors)
cat("Genus color assignments:\\n")
print(head(genus_colors))

# Create a data frame for the tree with taxonomy
tree_data <- data.frame(
  label = tree$tip.label,
  stringsAsFactors = FALSE
)

# Add taxonomy columns
tree_data$Phylum <- NA
tree_data$Genus <- NA

# Fill in taxonomic info where available
for (i in 1:nrow(tree_data)) {{
  tip <- tree_data$label[i]
  idx <- match(tip, metadata[[id_col]])
  if (!is.na(idx)) {{
    tree_data$Phylum[i] <- metadata$Phylum[idx]
    tree_data$Genus[i] <- metadata$Genus[idx]
  }}
}}

# Print sample data for debugging
cat("Sample of tree data with taxonomy:\\n")
print(head(tree_data))

# Function to create PDF output
create_circular_tree <- function() {{
  # Create circular tree visualization with colored tip points
  circular_file <- paste0(output_prefix, "_circular_tips.pdf")
  pdf(circular_file, width = 10, height = 10)
  
  # Create basic tree
  p <- ggtree(tree, layout = "circular", branch.length = "none") + 
       ggtitle("Phylogenetic Tree with Colored Tips") +
       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Add data
  p <- p %<+% tree_data
  
  # Color tips by phylum
  p <- p + geom_tippoint(aes(color = Phylum), size = 3, na.rm = TRUE) +
         scale_color_manual(values = phylum_colors, na.value = "gray70", 
                           name = "Phylum")
  
  # Improve legend
  p <- p + theme(legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10))
  
  # Print plot
  print(p)
  dev.off()
  
  # Return file path
  return(circular_file)
}}

# Function to create a rectangular tree with tip labels
create_rectangular_tree <- function() {{
  # Create rectangular tree with labels
  rect_file <- paste0(output_prefix, "_rectangular.pdf")
  pdf(rect_file, width = 12, height = 8)
  
  # Create rectangular tree
  p <- ggtree(tree, layout = "rectangular", branch.length = "none") + 
       ggtitle("Phylogenetic Tree (Rectangular)") +
       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Add data
  p <- p %<+% tree_data
  
  # Add colored tip points and labels
  p <- p + geom_tippoint(aes(color = Phylum), size = 3, na.rm = TRUE) +
         geom_tiplab(align = TRUE, size = 3, hjust = -0.05) +
         scale_color_manual(values = phylum_colors, na.value = "gray70", 
                           name = "Phylum")
  
  # Add legend
  p <- p + theme(legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10))
  
  # Print plot
  print(p)
  dev.off()
  
  # Return file path
  return(rect_file)
}}

# Function to create a heatmap visualization
create_heatmap_tree <- function() {{
  # Create tree with heatmap
  heatmap_file <- paste0(output_prefix, "_heatmap.pdf")
  pdf(heatmap_file, width = 12, height = 8)
  
  # Create rectangular tree
  p <- ggtree(tree, layout = "rectangular", branch.length = "none") + 
       ggtitle("Phylogenetic Tree with Taxonomy Heatmap") +
       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Add data
  p <- p %<+% tree_data
  
  # Create a simple heatmap matrix with just phylum
  phylum_mat <- matrix(0, nrow = length(tree$tip.label), ncol = length(phyla))
  rownames(phylum_mat) <- tree$tip.label
  colnames(phylum_mat) <- phyla
  
  for (i in 1:nrow(tree_data)) {{
    if (!is.na(tree_data$Phylum[i])) {{
      phy <- tree_data$Phylum[i]
      if (phy %in% phyla) {{
        phylum_mat[tree_data$label[i], phy] <- 1
      }}
    }}
  }}
  
  # Convert to data frame for gheatmap
  phylum_df <- as.data.frame(phylum_mat)
  
  # Add heatmap
  p2 <- gheatmap(p, phylum_df, offset = 0.05, width = 0.2, 
             colnames_position = "top", colnames_angle = 45, 
             colnames_offset_y = 0) +
       scale_fill_gradient(low = "white", high = "steelblue")
  
  # Print plot
  print(p2)
  dev.off()
  
  # Return file path
  return(heatmap_file)
}}

# Create PNG version of circular tree 
create_tree_png <- function() {{
  png_file <- paste0(output_prefix, "_circular_tips.png")
  png(png_file, width = 2000, height = 2000, res = 300)
  
  # Create tree
  p <- ggtree(tree, layout = "circular", branch.length = "none") + 
       ggtitle("Phylogenetic Tree with Colored Tips") +
       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Add data
  p <- p %<+% tree_data
  
  # Color tips by phylum
  p <- p + geom_tippoint(aes(color = Phylum), size = 3, na.rm = TRUE) +
         scale_color_manual(values = phylum_colors, na.value = "gray70", 
                           name = "Phylum")
  
  # Improve legend
  p <- p + theme(legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10))
  
  # Print plot
  print(p)
  dev.off()
  
  # Return file path
  return(png_file)
}}

# Create each type of visualization
circular_file <- ""
rect_file <- ""
heatmap_file <- ""
png_file <- ""

tryCatch({{
  cat("Creating circular tree...\\n")
  circular_file <- create_circular_tree()
  cat("Circular tree created at", circular_file, "\\n")
}}, error = function(e) {{
  cat("Error creating circular tree:", e$message, "\\n")
}})

tryCatch({{
  cat("Creating rectangular tree...\\n")
  rect_file <- create_rectangular_tree()
  cat("Rectangular tree created at", rect_file, "\\n")
}}, error = function(e) {{
  cat("Error creating rectangular tree:", e$message, "\\n")
}})

tryCatch({{
  cat("Creating heatmap tree...\\n")
  heatmap_file <- create_heatmap_tree()
  cat("Heatmap tree created at", heatmap_file, "\\n")
}}, error = function(e) {{
  cat("Error creating heatmap tree:", e$message, "\\n")
}})

tryCatch({{
  cat("Creating PNG version...\\n")
  png_file <- create_tree_png()
  cat("PNG version created at", png_file, "\\n")
}}, error = function(e) {{
  cat("Error creating PNG version:", e$message, "\\n")
}})

# Output completion message
cat("Visualizations completed!\\n")
cat("Created files:\\n")
if (circular_file != "") cat("  ", circular_file, "\\n")
if (rect_file != "") cat("  ", rect_file, "\\n")
if (heatmap_file != "") cat("  ", heatmap_file, "\\n")
if (png_file != "") cat("  ", png_file, "\\n")

quit(status=0)
"""
        
        # Write the R script to file
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        # Make the script executable
        import os
        os.chmod(script_path, 0o755)
        
        # Run the R script with better error handling
        import subprocess
        import logging
        logging.info(f"Running minimal ggtree script: {script_path}")
        
        # Define a function to run the command and capture all output
        def run_command(cmd):
            try:
                logging.info(f"Executing command: {cmd}")
                
                # Run with full output capture
                result = subprocess.run(
                    cmd, 
                    shell=True, 
                    check=False,  # Don't raise exception on non-zero exit
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                
                # Log all stdout
                for line in result.stdout.split('\n'):
                    if line.strip():
                        logging.info(f"R OUTPUT: {line}")
                
                # Log all stderr
                for line in result.stderr.split('\n'):
                    if line.strip():
                        logging.error(f"R ERROR: {line}")
                
                # Check exit code
                if result.returncode != 0:
                    logging.error(f"R script exited with code {result.returncode}")
                    return False
                
                return True
                
            except Exception as e:
                logging.error(f"Error running command: {e}")
                return False
        
        # Run the R script
        success = run_r_command_with_conda(script_path)
        
        # Check if output files were created
        output_files = [
            f"{output_prefix}_circular_tips.pdf",
            f"{output_prefix}_rectangular.pdf",
            f"{output_prefix}_heatmap.pdf",
            f"{output_prefix}_circular_tips.png"
        ]
        
        files_exist = any(os.path.exists(f) for f in output_files)
        
        if files_exist:
            logging.info(f"Successfully created ggtree visualizations using minimal functions")
            return True
        else:
            logging.warning(f"No output files were created")
            return False
            
    except Exception as e:
        logging.error(f"Error creating tree visualizations: {e}")
        return False

"""
Improved Phylogenetic Tree Visualization with Proper Spacing
This script fixes the issue where tree lines merge with taxonomic rings by:
1. Creating more space between the inner tree and the outer rings
2. Improving overall layout and appearance
"""

import os
import subprocess
import logging
from MetaMAG.utils import ensure_directory_exists

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_improved_visualization(tree_file, metadata_file, output_prefix):
    """
    Create publication-quality phylogenetic tree visualizations with improved spacing.
    
    Parameters:
    - tree_file: Path to the Newick tree file
    - metadata_file: Path to the metadata CSV file
    - output_prefix: Output directory and prefix for saving files
    
    Returns:
    - Boolean indicating success or failure
    """
    try:
        # Create R script for visualization
        script_path = f"{output_prefix}_improved_vis.R"
        
        # Write R script content with fixes for spacing issues
        r_script = f"""#!/usr/bin/env Rscript

# Load required packages with error handling
tryCatch({{
  # Load common packages first
  library(ape)
  library(RColorBrewer)
  cat("Base packages loaded successfully\\n")
  
  # Try to load ggtree and dependencies
  if (!require("ggplot2")) {{
    install.packages("ggplot2", repos="https://cloud.r-project.org")
    library(ggplot2)
  }}
  
  use_ggtree <- TRUE
  tryCatch({{
    if (!require("ggtree")) {{
      if (!require("BiocManager")) {{
        install.packages("BiocManager", repos="https://cloud.r-project.org")
      }}
      BiocManager::install("ggtree")
      library(ggtree)
    }}
    cat("ggtree packages loaded successfully\\n")
  }}, error = function(e) {{
    cat("Could not load ggtree packages, falling back to base R:", e$message, "\\n")
    use_ggtree <<- FALSE
  }})
}}, error = function(e) {{
  cat("Error loading required packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"

# Ensure reproducible colors
set.seed(42)

# Read tree with error handling
cat("Reading tree from", tree_file, "\\n")
tryCatch({{
  tree <- read.tree(tree_file)
  if (is.null(tree)) {{
    stop("Failed to read tree file")
  }}
}}, error = function(e) {{
  cat("Error reading tree, attempting to fix format...\\n")
  tree_text <- readLines(tree_file)
  # Try to fix common issues
  tree_text <- gsub("\\\\)[a-zA-Z0-9]*:", "\\\\):", tree_text)
  temp_file <- tempfile(fileext = ".tre")
  writeLines(tree_text, temp_file)
  tree <<- read.tree(temp_file)
  
  if (is.null(tree)) {{
    stop("Failed to read tree file after attempted fixes")
  }}
}})

cat("Successfully read tree with", length(tree$tip.label), "tips\\n")

# Ensure tree is rooted
if (!is.rooted(tree)) {{
  cat("Tree is not rooted. Applying midpoint rooting...\\n")
  tree <- midpoint(tree)
}}

# Read and process metadata
cat("Processing metadata from", metadata_file, "\\n")
if (endsWith(metadata_file, ".tsv")) {{
  metadata <- read.delim(metadata_file, sep="\\t", check.names=FALSE)
}} else {{
  metadata <- read.csv(metadata_file, check.names=FALSE)
}}

# Identify ID column
id_cols <- c("user_genome", "genome", "bin_id", "name", "Genome")
id_col <- NULL
for (col in id_cols) {{
  if (col %in% colnames(metadata)) {{
    id_col <- col
    break
  }}
}}

if (is.null(id_col)) {{
  id_col <- colnames(metadata)[1]
  cat("Using", id_col, "as ID column\\n")
}}

# Extract taxonomy information
tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]])

# Create taxonomy columns
for (level in tax_levels) {{
  # Check if column exists in metadata
  if (level %in% colnames(metadata)) {{
    taxa_data[[level]] <- metadata[[level]]
  }} else if (tolower(level) %in% colnames(metadata)) {{
    # Try lowercase version
    taxa_data[[level]] <- metadata[[tolower(level)]]
  }} else if ("classification" %in% colnames(metadata)) {{
    # Initialize with Unknown
    taxa_data[[level]] <- "Unknown"
  }} else {{
    # Default to Unknown
    taxa_data[[level]] <- "Unknown"
  }}
  
  # Clean up values and handle NAs
  taxa_data[[level]] <- as.character(taxa_data[[level]])
  taxa_data[[level]] <- gsub("^[a-z]__", "", taxa_data[[level]])
  taxa_data[[level]][is.na(taxa_data[[level]]) | taxa_data[[level]] == ""] <- "Unknown"
}}

# If classification column exists, parse it
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

# Function to match tree tips with taxa data
match_tip_to_taxa <- function(tip, taxa_data) {{
  # Try exact match first
  match_idx <- which(taxa_data$ID == tip)
  
  if (length(match_idx) == 0) {{
    # Try case-insensitive match
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  
  if (length(match_idx) == 0) {{
    # Try partial matching
    for (i in 1:nrow(taxa_data)) {{
      if (grepl(taxa_data$ID[i], tip, fixed=TRUE) || 
          grepl(tip, taxa_data$ID[i], fixed=TRUE)) {{
        match_idx <- i
        break
      }}
    }}
  }}
  
  if (length(match_idx) > 0) {{
    return(taxa_data[match_idx[1], ])
  }} else {{
    # Return Unknown values if no match
    unknown_row <- taxa_data[1, ]
    unknown_row[] <- c(tip, rep("Unknown", length(tax_levels)))
    return(unknown_row)
  }}
}}

# Match tree tips with taxa data
tip_data <- do.call(rbind, lapply(tree$tip.label, function(tip) {{
  return(match_tip_to_taxa(tip, taxa_data))
}}))

# Get unique taxa
phyla <- unique(tip_data$Phylum)
phyla <- phyla[phyla != "Unknown"]
phyla <- phyla[!is.na(phyla)]

genera <- unique(tip_data$Genus)
genera <- genera[genera != "Unknown"]
genera <- genera[!is.na(genera)]

# Create robust color palettes
if (length(phyla) > 0) {{
  # For phyla, use maximum 9 distinct colors from Set1
  if (length(phyla) <= 9) {{
    phylum_colors <- brewer.pal(length(phyla), "Set1")
  }} else {{
    # For more phyla, create a gradient of colors
    phylum_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(phyla))
  }}
  names(phylum_colors) <- phyla
  
  # Add color for Unknown
  phylum_colors <- c(phylum_colors, Unknown = "lightgray")
}} else {{
  phylum_colors <- c(Unknown = "lightgray")
}}

if (length(genera) > 0) {{
  # For genera, use maximum 8 distinct colors from Set2
  if (length(genera) <= 8) {{
    genus_colors <- brewer.pal(length(genera), "Set2")
  }} else {{
    # For more genera, create a gradient of colors
    genus_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(genera))
  }}
  names(genus_colors) <- genera
  
  # Add color for Unknown
  genus_colors <- c(genus_colors, Unknown = "lightgray")
}} else {{
  genus_colors <- c(Unknown = "lightgray")
}}

# Create maps from tip labels to taxonomy
tip_to_phylum <- setNames(tip_data$Phylum, tip_data$ID)
tip_to_genus <- setNames(tip_data$Genus, tip_data$ID)

# Function to get taxon with proper fallback
get_taxon <- function(tip_label, taxon_map) {{
  if (tip_label %in% names(taxon_map)) {{
    return(taxon_map[tip_label])
  }} else {{
    cat("Warning: Tip", tip_label, "not found in taxonomy map\\n")
    return("Unknown")
  }}
}}

# Function to get color for a taxon with proper error handling
get_color <- function(value, palette) {{
  if (is.null(value) || is.na(value) || value == "") {{
    return(palette["Unknown"])
  }} else if (value %in% names(palette)) {{
    return(palette[value])
  }} else {{
    cat("Warning: No color found for", value, "\\n")
    return(palette["Unknown"])
  }}
}}

##############################################
# VISUALIZATION WITH GGTREE (if available)
##############################################
if (use_ggtree) {{
  cat("Creating improved ggtree visualizations...\\n")
  
  # CIRCULAR TREE WITH IMPROVED SPACING
  circular_file <- paste0(output_prefix, "_improved_circular.pdf")
  pdf(circular_file, width = 14, height = 14)  # Increased size for better visibility
  
  tryCatch({{
    # Create tree data with taxonomic information
    tree_data <- data.frame(
      label = tree$tip.label,
      stringsAsFactors = FALSE
    )
    
    # Add taxonomy columns
    tree_data$Phylum <- NA
    tree_data$Genus <- NA
    
    # Fill in taxonomic info where available
    for (i in 1:nrow(tree_data)) {{
      tip <- tree_data$label[i]
      idx <- match(tip, tip_data$ID)
      if (!is.na(idx)) {{
        tree_data$Phylum[i] <- tip_data$Phylum[idx]
        tree_data$Genus[i] <- tip_data$Genus[idx]
      }}
    }}
    
    # Create a basic tree plot with SIGNIFICANTLY DECREASED outer_radius to create more space
    # The key improvement is setting a much smaller outer_radius (0.5) to leave plenty of space for rings
    p <- ggtree(tree, layout="circular", branch.length="none", open.angle=15, 
               outer.radius=0.5) +  # Further reduced to ensure rings are fully visible
         ggtitle("Phylogenetic Tree with Taxonomic Rings") +
         theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    # Add data
    p <- p %<+% tree_data
    
    # Add clear tip points with consistent colors
    p <- p + geom_tippoint(aes(color=Phylum), size=3, na.rm=TRUE) +
           scale_color_manual(values=phylum_colors, 
                             na.value="gray70", 
                             name="Phylum")
    
    # Create manually positioned rings with clear separation
    # First, get the data needed for creating rings
    dat <- p$data
    tip_data <- dat[dat$isTip, ]
    
    # Calculate angles for each tip
    tip_data$angle <- atan2(tip_data$y, tip_data$x)
    # Fix negative angles
    tip_data$angle[tip_data$angle < 0] <- tip_data$angle[tip_data$angle < 0] + 2 * pi
    
    # Calculate max radius of tree
    max_radius <- max(sqrt(tip_data$x^2 + tip_data$y^2))
    
    # Sort tips by angle
    tip_data <- tip_data[order(tip_data$angle), ]
    
    # Create phylum ring with much more space from tree
    ring_gap = 0.2 * max_radius  # Increased from 10% to 20% gap between tree and rings
    phylum_inner <- max_radius + ring_gap
    phylum_outer <- phylum_inner + 0.18 * max_radius  # Increased from 15% to 18% for wider ring
    
    # Create genus ring with larger gap from phylum ring
    genus_inner <- phylum_outer + 0.08 * max_radius  # Increased from 2% to 8% gap between rings
    genus_outer <- genus_inner + 0.18 * max_radius  # Increased from 15% to 18% for wider ring
    
    # Function to create arc data for a specific ring
    create_arc_data <- function(inner_r, outer_r, tip_data) {{
      n_tips <- nrow(tip_data)
      data <- data.frame()
      
      for (i in 1:(n_tips)) {{
        j <- if (i < n_tips) i + 1 else 1
        
        start_angle <- tip_data$angle[i]
        end_angle <- tip_data$angle[j]
        
        # Ensure proper direction for arcs (clockwise)
        if (end_angle < start_angle) end_angle <- end_angle + 2 * pi
        
        # Create more points for smoother arcs
        angles <- seq(start_angle, end_angle, length.out = 50)
        
        # Inner arc points
        xi <- inner_r * cos(angles)
        yi <- inner_r * sin(angles)
        
        # Outer arc points (in reverse to create a polygon)
        xo <- outer_r * cos(rev(angles))
        yo <- outer_r * sin(rev(angles))
        
        # Create a segment for this tip
        segment <- data.frame(
          x = c(xi, xo),
          y = c(yi, yo),
          id = i,
          taxon = tip_data$Phylum[i],  # Use corresponding taxonomy
          stringsAsFactors = FALSE
        )
        
        data <- rbind(data, segment)
      }}
      
      return(data)
    }}
    
    # Create phylum arcs data
    phylum_data <- create_arc_data(phylum_inner, phylum_outer, tip_data)
    phylum_data$taxon <- sapply(phylum_data$id, function(i) {{
      get_taxon(tip_data$label[i], tip_to_phylum)
    }})
    
    # Create genus arcs data
    genus_data <- create_arc_data(genus_inner, genus_outer, tip_data)
    genus_data$taxon <- sapply(genus_data$id, function(i) {{
      get_taxon(tip_data$label[i], tip_to_genus)
    }})
    
    # Add phylum ring with better spacing
    p <- p + geom_polygon(data=phylum_data, 
                          aes(x=x, y=y, group=id, fill=taxon),
                          alpha=0.8, color="darkgray", size=0.2) +
             scale_fill_manual(values=phylum_colors, name="Phylum")
    
    # Add genus ring with better spacing
    p <- p + geom_polygon(data=genus_data, 
                         aes(x=x, y=y, group=id, fill=taxon),
                         alpha=0.8, color="darkgray", size=0.2) +
             new_scale_fill() +
             scale_fill_manual(values=genus_colors, name="Genus")
    
    # Add labels for the rings
    text_distance <- genus_outer + 0.05 * max_radius
    
    # Add proper legend
    p <- p + theme(legend.position = "right",
                   legend.box = "vertical",
                   legend.title = element_text(size=12, face="bold"),
                   plot.margin = margin(t=20, r=120, b=20, l=20, unit="pt"))
    
    # Extend plot limits to ensure ALL rings are completely visible
    # Use a much larger multiplier to guarantee no cropping
    max_dim <- genus_outer * 1.5  # Increased from 1.2 to 1.5 for full visibility
    p <- p + coord_fixed(xlim=c(-max_dim, max_dim), 
                        ylim=c(-max_dim, max_dim), 
                        clip="off") +
         # Add much larger margins to ensure rings aren't cut off
         theme(plot.margin = margin(t=50, r=150, b=50, l=50, unit="pt"))
    
    # Print plot
    print(p)
    
  }}, error = function(e) {{
    cat("Error in ggtree circular visualization:", e$message, "\\n")
  }})
  
  dev.off()
  
  # RECTANGULAR TREE WITH IMPROVED SPACING
  rectangular_file <- paste0(output_prefix, "_improved_rectangular.pdf")
  height <- max(8, length(tree$tip.label) * 0.25)
  pdf(rectangular_file, width = 14, height = height)
  
  tryCatch({{
    # Create tree data
    tree_data <- data.frame(
      label = tree$tip.label,
      stringsAsFactors = FALSE
    )
    
    # Add taxonomy columns
    tree_data$Phylum <- NA
    tree_data$Genus <- NA
    
    # Fill in taxonomic info
    for (i in 1:nrow(tree_data)) {{
      tip <- tree_data$label[i]
      idx <- match(tip, tip_data$ID)
      if (!is.na(idx)) {{
        tree_data$Phylum[i] <- tip_data$Phylum[idx]
        tree_data$Genus[i] <- tip_data$Genus[idx]
      }}
    }}
    
    # Create rectangular tree with improved spacing
    p <- ggtree(tree, branch.length="none") + 
         ggtitle("Phylogenetic Tree (Rectangular)") +
         theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    # Add data
    p <- p %<+% tree_data
    
    # Add colored tip points with more space
    p <- p + 
      # Note the increased offset for more space
      geom_tippoint(aes(color=Phylum), size=3, na.rm=TRUE) +
      scale_color_manual(values=phylum_colors, na.value="gray70", name="Phylum") +
      # Add labels with more space
      geom_tiplab(align=TRUE, size=3, hjust=-0.1)  # Increased hjust for more space
    
    # Add legend
    p <- p + theme(legend.position = "right",
                  legend.title = element_text(size = 12, face = "bold"),
                  legend.text = element_text(size = 10),
                  plot.margin = margin(t=20, r=50, b=20, l=20, unit="pt"))
    
    # Print the rectangular tree
    print(p)
  }}, error = function(e) {{
    cat("Error in ggtree rectangular visualization:", e$message, "\\n")
  }})
  
  dev.off()
  
  cat("Improved ggtree visualizations completed\\n")
}}

##############################################
# VISUALIZATION WITH BASE R (always created as backup)
##############################################
cat("Creating improved base R visualizations...\\n")

# Improved function to draw a colored arc segment with better visibility
draw_arc_segment <- function(center_x, center_y, radius_inner, radius_outer, 
                           start_angle, end_angle, color) {{
  # Ensure angles are in correct order
  if (start_angle > end_angle) {{
    tmp <- start_angle
    start_angle <- end_angle
    end_angle <- tmp
  }}
  
  # Ensure minimum arc size for visibility
  if ((end_angle - start_angle) < 0.01) {{
    end_angle <- start_angle + 0.01
  }}
  
  # Create more points for smoother arcs
  angles <- seq(start_angle, end_angle, length.out = max(100, ceiling((end_angle - start_angle) * 200/pi)))
  
  x_inner <- center_x + radius_inner * cos(angles)
  y_inner <- center_y + radius_inner * sin(angles)
  x_outer <- center_x + radius_outer * cos(rev(angles))
  y_outer <- center_y + radius_outer * sin(rev(angles))
  
  # Add border for better visualization
  polygon(c(x_inner, x_outer), c(y_inner, y_outer), 
          col = color, border = "darkgray", lwd = 0.5)
}}

# Create circular tree with taxonomic rings - IMPROVED VERSION WITH FULL RING VISIBILITY
circular_base_file <- paste0(output_prefix, "_improved_circular_base.pdf")
pdf(circular_base_file, width = 18, height = 18)  # Further increased size for complete ring visibility

# Set up layout with more space for the tree
layout(matrix(c(1,2), nrow=1), widths=c(5,1))

# Plot area 1: The tree with rings - with much larger margins to prevent cropping
par(mar=c(8,8,8,8), xpd=TRUE)  # Significantly increased margins on all sides

# Draw tree in circular layout with reduced size to create space
tryCatch({{
  # KEY IMPROVEMENT: Reduce tree size to 80% to leave space for rings
  plot(tree, type="fan", show.tip.label=FALSE, main="Phylogenetic Tree with Taxonomic Rings", 
       edge.width=1.2, cex.main=1.5, edge.color="black", 
       no.margin=FALSE, show.grid=FALSE, use.edge.length=FALSE,
       # Add scale factor to reduce tree size (key improvement)
       x.lim=c(-1.2, 1.2), y.lim=c(-1.2, 1.2))
}}, error = function(e) {{
  cat("Error plotting circular tree:", e$message, "\\n")
  plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", 
       main="Error plotting tree", cex.main=1.5)
  return()
}})

# Get the coordinates for tips from the last plotted tree
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords_x <- last_plot$xx[1:length(tree$tip.label)]
tip_coords_y <- last_plot$yy[1:length(tree$tip.label)]

# Calculate center of plot and radius
center_x <- mean(range(tip_coords_x))
center_y <- mean(range(tip_coords_y))
max_radius <- max(sqrt((tip_coords_x - center_x)^2 + (tip_coords_y - center_y)^2))

# Calculate angles for each tip
tip_angles <- atan2(tip_coords_y - center_y, tip_coords_x - center_x)
# Fix negative angles
tip_angles[tip_angles < 0] <- tip_angles[tip_angles < 0] + 2 * pi

# Sort tips by angle
sorted_indices <- order(tip_angles)
sorted_tips <- tree$tip.label[sorted_indices]
sorted_angles <- tip_angles[sorted_indices]

# Complete the circle by adding the first tip at the end
sorted_tips <- c(sorted_tips, sorted_tips[1])
sorted_angles <- c(sorted_angles, sorted_angles[1] + 2*pi)

# Draw rings for Phylum and Genus with GREATLY IMPROVED SPACING
# KEY IMPROVEMENT: Add larger gap between tree and first ring
ring_gap <- 0.2 * max_radius  # Increased from 15% to 20% gap between tree and rings

# Set inner and outer radius for each ring with more space
phylum_inner <- max_radius + ring_gap
phylum_outer <- phylum_inner + 0.18 * max_radius  # Slightly wider phylum ring

# Add larger gap between rings
ring_separation <- 0.08 * max_radius  # Increased from 5% to 8% gap between rings
genus_inner <- phylum_outer + ring_separation
genus_outer <- genus_inner + 0.18 * max_radius  # Slightly wider genus ring too

# Background circles for each ring
angles <- seq(0, 2*pi, length.out=200)
# Phylum ring background
polygon(center_x + phylum_inner * cos(angles), 
        center_y + phylum_inner * sin(angles),
        center_x + phylum_outer * cos(rev(angles)), 
        center_y + phylum_outer * sin(rev(angles)),
        col="lightgray", border="darkgray")
# Genus ring background
polygon(center_x + genus_inner * cos(angles), 
        center_y + genus_inner * sin(angles),
        center_x + genus_outer * cos(rev(angles)), 
        center_y + genus_outer * sin(rev(angles)),
        col="lightgray", border="darkgray")

# Draw segments for each tip
for (i in 1:(length(sorted_tips)-1)) {{
  tryCatch({{
    tip <- sorted_tips[i]
    next_tip <- sorted_tips[i+1]
    start_angle <- sorted_angles[i]
    end_angle <- sorted_angles[i+1]
    
    # Get phylum and genus
    phylum <- get_taxon(tip, tip_to_phylum)
    genus <- get_taxon(tip, tip_to_genus)
    
    # Get colors safely
    phylum_color <- get_color(phylum, phylum_colors)
    genus_color <- get_color(genus, genus_colors)
    
    # Draw phylum segment
    draw_arc_segment(center_x, center_y, phylum_inner, phylum_outer, 
                    start_angle, end_angle, phylum_color)
    
    # Draw genus segment
    draw_arc_segment(center_x, center_y, genus_inner, genus_outer, 
                    start_angle, end_angle, genus_color)
  }}, error = function(e) {{
    cat("Error drawing segment", i, ":", e$message, "\\n")
  }})
}}

# Add colored tip points with more visibility
for (i in 1:length(tree$tip.label)) {{
  tryCatch({{
    tip <- tree$tip.label[i]
    phylum <- get_taxon(tip, tip_to_phylum)
    phylum_color <- get_color(phylum, phylum_colors)
    
    points(tip_coords_x[i], tip_coords_y[i], pch=19, col=phylum_color, cex=1.5)
    # Add outline for better visibility
    points(tip_coords_x[i], tip_coords_y[i], pch=1, col="black", cex=1.5, lwd=0.5)
  }}, error = function(e) {{
    cat("Error drawing tip point", i, ":", e$message, "\\n")
  }})
}}

# Draw connecting lines between tree tips and rings for clarity
# KEY IMPROVEMENT: Add radial lines to better separate tree from rings
for (i in 1:length(sorted_tips) - 1) {{
  tip_idx <- which(tree$tip.label == sorted_tips[i])
  if (length(tip_idx) > 0) {{
    angle <- sorted_angles[i]
    x1 <- tip_coords_x[tip_idx]
    y1 <- tip_coords_y[tip_idx]
    x2 <- center_x + phylum_inner * cos(angle)
    y2 <- center_y + phylum_inner * sin(angle)
    
    # Draw a thin, light gray line connecting tip to ring
    lines(c(x1, x2), c(y1, y2), col="#CCCCCC", lwd=0.5, lty=3)
  }}
}}

# Plot area 2: Legends
par(mar=c(2,1,2,1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Draw phylum legend
text(0.1, 0.95, "Phylum", font=2, pos=4, cex=1.2)
legend_y <- 0.9
for (i in 1:length(phyla)) {{
  rect(0.1, legend_y - 0.015, 0.2, legend_y + 0.015, 
       col=phylum_colors[phyla[i]], border="black")
  text(0.22, legend_y, phyla[i], pos=4)
  legend_y <- legend_y - 0.05
}}
# Add unknown entry
rect(0.1, legend_y - 0.015, 0.2, legend_y + 0.015, 
     col="lightgray", border="black")
text(0.22, legend_y, "Unknown", pos=4)

# Draw genus legend
text(0.1, 0.55, "Genus", font=2, pos=4, cex=1.2)
legend_y <- 0.5
max_genus_to_show <- min(length(genera), 15)  # Limit genera in legend
for (i in 1:max_genus_to_show) {{
  rect(0.1, legend_y - 0.015, 0.2, legend_y + 0.015, 
       col=genus_colors[genera[i]], border="black")
  text(0.22, legend_y, genera[i], pos=4)
  legend_y <- legend_y - 0.03
}}
# Add unknown entry
rect(0.1, legend_y - 0.015, 0.2, legend_y + 0.015, 
     col="lightgray", border="black")
text(0.22, legend_y, "Unknown", pos=4)

# Add note if too many genera
if (length(genera) > max_genus_to_show) {{
  text(0.5, legend_y - 0.05, 
      paste("(+", length(genera) - max_genus_to_show, "more genera)"), 
      cex=0.8, col="darkgray")
}}

dev.off()

# Create rectangular tree with improved spacing
rectangular_base_file <- paste0(output_prefix, "_improved_rectangular_base.pdf")
height <- max(8, length(tree$tip.label) * 0.25)
pdf(rectangular_base_file, width = 14, height = height)

# Set up layout with more space for the tree
layout(matrix(c(1,2), nrow=1), widths=c(4,1))

# Plot tree in rectangular layout with improved spacing
par(mar=c(2,1,3,1), xpd=TRUE)
tryCatch({{
  # KEY IMPROVEMENT: Use direction="right" to get more horizontal space
  plot(tree, type="phylogram", direction="right", 
       show.tip.label=TRUE, main="Phylogenetic Tree (Rectangular)", 
       edge.width=1.2, cex.main=1.5, cex=0.7, 
       label.offset=0.05)  # Increased offset for better spacing
}}, error = function(e) {{
  cat("Error plotting rectangular tree:", e$message, "\\n")
  plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", 
       main="Error plotting tree", cex.main=1.5)
  return()
}})

# Add colored points to tips with more space
for (i in 1:length(tree$tip.label)) {{
  tryCatch({{
    tip <- tree$tip.label[i]
    phylum <- get_taxon(tip, tip_to_phylum)
    phylum_color <- get_color(phylum, phylum_colors)
    
    # KEY IMPROVEMENT: Increased offset for better spacing
    tiplabels(pch=19, col=phylum_color, cex=1.5, offset=0.3, frame="none", 
              tip=which(tree$tip.label == tip))
    # Add outline for better visibility
    tiplabels(pch=1, col="black", cex=1.5, offset=0.3, frame="none", 
              tip=which(tree$tip.label == tip))
  }}, error = function(e) {{
    cat("Error drawing tip", i, ":", e$message, "\\n")
  }})
}}

# Plot legends in second panel
par(mar=c(1,1,3,1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Draw phylum legend
text(0.1, 0.95, "Phylum", font=2, pos=4, cex=1.2)
legend_y <- 0.9
for (i in 1:length(phyla)) {{
  points(0.15, legend_y, pch=19, col=phylum_colors[phyla[i]], cex=1.5)
  points(0.15, legend_y, pch=1, col="black", cex=1.5)
  text(0.22, legend_y, phyla[i], pos=4)
  legend_y <- legend_y - 0.05
}}
# Add unknown entry
points(0.15, legend_y, pch=19, col="lightgray", cex=1.5)
points(0.15, legend_y, pch=1, col="black", cex=1.5)
text(0.22, legend_y, "Unknown", pos=4)

# Draw genus legend
text(0.1, 0.55, "Genus", font=2, pos=4, cex=1.2)
legend_y <- 0.5
max_genus_to_show <- min(length(genera), 15)  # Limit genera in legend
for (i in 1:max_genus_to_show) {{
  points(0.15, legend_y, pch=19, col=genus_colors[genera[i]], cex=1.5)
  points(0.15, legend_y, pch=1, col="black", cex=1.5)
  text(0.22, legend_y, genera[i], pos=4)
  legend_y <- legend_y - 0.03
}}
# Add unknown entry
points(0.15, legend_y, pch=19, col="lightgray", cex=1.5)
points(0.15, legend_y, pch=1, col="black", cex=1.5)
text(0.22, legend_y, "Unknown", pos=4)

# Add note if too many genera
if (length(genera) > max_genus_to_show) {{
  text(0.5, legend_y - 0.05, 
      paste("(+", length(genera) - max_genus_to_show, "more genera)"), 
      cex=0.8, col="darkgray")
}}

dev.off()

cat("Base R visualizations completed successfully!\\n")
cat("Output files:\\n")
if (use_ggtree) {{
  cat("  ", circular_file, "\\n")
  cat("  ", rectangular_file, "\\n")
}}
cat("  ", circular_base_file, "\\n")
cat("  ", rectangular_base_file, "\\n")
"""
        
        # Write the R script to file
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        # Make the script executable
        import os
        os.chmod(script_path, 0o755)
        
        # Run the R script
        import subprocess
        import logging
        logging.info(f"Running improved R script: {script_path}")
        
        # Define a function to run the command with detailed logging        
        success = run_command(f"Rscript {script_path}")
        
        # Check if output files were created
        output_files = [
            f"{output_prefix}_improved_circular_base.pdf",
            f"{output_prefix}_improved_rectangular_base.pdf",
            f"{output_prefix}_improved_circular.pdf",
            f"{output_prefix}_improved_rectangular.pdf"
        ]
        
        files_exist = any(os.path.exists(f) for f in output_files)
        
        if files_exist:
            logging.info(f"Successfully created improved tree visualizations")
            return True
        else:
            logging.warning(f"No output files were created")
            return False
            
    except Exception as e:
        logging.error(f"Error creating tree visualizations: {e}")
        return False

def run_improved_visualization(tree_file, metadata_file, output_dir):
    """
    Main function to run the improved tree visualization.
    
    Parameters:
    - tree_file: Path to the tree file (.tree or .nwk format)
    - metadata_file: Path to the metadata file (CSV or TSV)
    - output_dir: Directory to save output files
    
    Returns:
    - Boolean indicating success or failure
    """
    try:
        # Ensure output directory exists
        ensure_directory_exists(output_dir)
        
        # Build output prefix
        tree_basename = os.path.basename(tree_file).split('.')[0]
        output_prefix = os.path.join(output_dir, f"{tree_basename}_improved")
        
        # Run the visualization
        success = create_improved_visualization(tree_file, metadata_file, output_prefix)
        
        if success:
            logger.info(f"Successfully created improved tree visualizations at {output_prefix}")
            return True
        else:
            logger.error("Failed to create improved tree visualizations")
            return False
            
    except Exception as e:
        logger.error(f"Error in improved visualization: {e}")
        return False

# The function can be used as a drop-in replacement for the original visualization
def create_publication_quality_visualization(tree_file, metadata_file, output_prefix):
    """
    Drop-in replacement for the original visualization function.
    This function adds proper spacing between the tree and rings.
    """
    return create_improved_visualization(tree_file, metadata_file, output_prefix)

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 6:
        print("Usage: python tree_visualization.py <cpus> <memory> <time_limit> <project_input> <project_output> [taxonomy_source]")
        sys.exit(1)
    
    # Parse command line arguments
    cpus = int(sys.argv[1])
    memory = sys.argv[2]
    time_limit = sys.argv[3]
    project_input = sys.argv[4] 
    project_output = sys.argv[5]
    taxonomy_source = None
    if len(sys.argv) > 6:
        taxonomy_source = sys.argv[6]
    
    # Run tree visualization
    success = run_tree_visualization(cpus, memory, time_limit, project_input, project_output, taxonomy_source)
    sys.exit(0 if success else 1)

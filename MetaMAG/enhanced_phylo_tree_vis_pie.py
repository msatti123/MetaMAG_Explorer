"""
Enhanced Phylogenetic Tree Visualization - UPDATED VERSION
- Support for multiple functional annotation types (CAZyme, COG, KEGG)
- Fixed MAG name matching between functional data and tree tips
- Multiple pie chart rings for different annotation types
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
from collections import defaultdict
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import setup_r_environment as config_setup_r, get_r_config, detect_r_environment

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_enhanced_tree_visualization(cpus, memory, time_limit, project_input, project_output, 
                                          taxonomy_source=None, novel_mags_file=None, 
                                          functional_annotations=None, annotation_type="auto",
                                          cazyme_annotations=None, cog_annotations=None, 
                                          kegg_annotations=None):
    """FIXED: Enhanced phylogenetic tree visualization with better error handling and data processing."""
    start_time = time.time()
    logger.info("=" * 60)
    logger.info("STARTING ENHANCED PHYLOGENETIC TREE VISUALIZATION - FIXED VERSION")
    logger.info("=" * 60)
    
    # Validate inputs
    if not project_output or not os.path.exists(project_output):
        logger.error(f"Project output directory does not exist: {project_output}")
        return False
    
    # Log all input parameters for debugging
    logger.info("Input parameters:")
    logger.info(f"  - Project output: {project_output}")
    logger.info(f"  - Taxonomy source: {taxonomy_source}")
    logger.info(f"  - Novel MAGs file: {novel_mags_file}")
    logger.info(f"  - Functional annotations: {functional_annotations}")
    logger.info(f"  - CAZyme annotations: {cazyme_annotations}")
    logger.info(f"  - COG annotations: {cog_annotations}")
    logger.info(f"  - KEGG annotations: {kegg_annotations}")
    
    # Set up R environment
    logger.info("Setting up R environment...")
    from MetaMAG.config import setup_r_environment as config_setup_r, get_r_config
    
    success = config_setup_r()
    if not success:
        logger.error("Failed to set up R environment")
        return False
    
    r_config = get_r_config()
    if r_config:
        logger.info(f"R environment configured: {r_config.get('env_name', 'unknown')}")
    
    # Define paths
    phylogeny_dir = os.path.join(project_output, "Phylogeny", "gtdbtk", "taxonomy2")
    trees_dir = os.path.join(phylogeny_dir, "trees")
    vis_output_dir = os.path.join(phylogeny_dir, "enhanced_visualizations")
    
    # Ensure output directory exists
    from MetaMAG.utils import ensure_directory_exists
    try:
        ensure_directory_exists(vis_output_dir)
        logger.info(f"Output directory: {vis_output_dir}")
    except Exception as e:
        logger.error(f"Failed to create output directory: {e}")
        return False
    
    # Import required functions
    from MetaMAG.enhanced_phylo_tree_vis import (
        find_tree_files, load_novel_mags_list_enhanced, get_tree_tip_names,
        find_taxonomy_files, process_metadata_file, create_metadata_from_tree,
        load_functional_annotations_enhanced, load_cazyme_annotations,
        load_cog_annotations, create_simplified_taxonomic_visualization
    )
    
    # Use the fixed functions from the previous artifact
    from MetaMAG.enhanced_phylo_tree_vis_fixed import (
        load_kegg_annotations_fixed
    )
    
    # Find tree files
    logger.info("Searching for phylogenetic tree files...")
    tree_files = find_tree_files(trees_dir, project_output)
    
    if not tree_files:
        logger.error("No tree files found!")
        return False
    
    logger.info(f"Found {len(tree_files)} tree file(s):")
    for tree_file in tree_files:
        logger.info(f"  - {tree_file}")
        # Verify file exists and is readable
        if not os.path.exists(tree_file):
            logger.error(f"Tree file does not exist: {tree_file}")
            return False
        if os.path.getsize(tree_file) == 0:
            logger.error(f"Tree file is empty: {tree_file}")
            return False
    
    # Process novel MAGs list
    novel_mags = []
    if novel_mags_file and os.path.exists(novel_mags_file):
        novel_mags = load_novel_mags_list_enhanced(novel_mags_file)
        logger.info(f"Loaded {len(novel_mags)} novel MAGs for highlighting")
        if len(novel_mags) > 0:
            logger.info(f"First few novel MAGs: {novel_mags[:5]}")
    else:
        logger.info("No novel MAGs file provided or file doesn't exist")
    
    # Get tree tip names for matching
    tree_tips = get_tree_tip_names(tree_files[0])
    logger.info(f"Found {len(tree_tips)} tree tip names")
    if len(tree_tips) > 0:
        logger.info(f"Sample tree tips: {tree_tips[:min(5, len(tree_tips))]}")
    
    # Process functional annotations - FIXED to handle all three types properly
    functional_data = {}
    
    # Load CAZyme annotations
    if cazyme_annotations and os.path.exists(cazyme_annotations):
        logger.info(f"Loading CAZyme annotations from: {cazyme_annotations}")
        try:
            cazyme_data = load_cazyme_annotations(cazyme_annotations, tree_tips)
            if cazyme_data:
                functional_data['cazyme'] = cazyme_data
                logger.info(f"? Loaded CAZyme data for {len(cazyme_data)} MAGs")
            else:
                logger.warning("No CAZyme data loaded")
        except Exception as e:
            logger.error(f"Error loading CAZyme annotations: {e}")
    
    # Load COG annotations
    if cog_annotations and os.path.exists(cog_annotations):
        logger.info(f"Loading COG annotations from: {cog_annotations}")
        try:
            cog_data = load_cog_annotations(cog_annotations, tree_tips)
            if cog_data:
                functional_data['cog'] = cog_data
                logger.info(f"? Loaded COG data for {len(cog_data)} MAGs")
            else:
                logger.warning("No COG data loaded")
        except Exception as e:
            logger.error(f"Error loading COG annotations: {e}")
    
    # Load KEGG annotations - USING FIXED FUNCTION
    if kegg_annotations and os.path.exists(kegg_annotations):
        logger.info(f"Loading KEGG annotations from: {kegg_annotations}")
        try:
            kegg_data = load_kegg_annotations_fixed(kegg_annotations, tree_tips)
            if kegg_data:
                functional_data['kegg'] = kegg_data
                logger.info(f"? Loaded KEGG data for {len(kegg_data)} MAGs")
            else:
                logger.warning("No KEGG data loaded")
        except Exception as e:
            logger.error(f"Error loading KEGG annotations: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    # Find or create taxonomy metadata
    metadata_file = None
    
    if taxonomy_source and os.path.exists(taxonomy_source):
        logger.info(f"Using provided taxonomy source: {taxonomy_source}")
        try:
            metadata_file = process_metadata_file(
                taxonomy_source, 
                os.path.join(vis_output_dir, "Metadata_For_Tree.csv")
            )
            if metadata_file:
                logger.info(f"? Successfully processed taxonomy metadata: {metadata_file}")
        except Exception as e:
            logger.error(f"Error processing taxonomy source: {e}")
    
    if not metadata_file:
        logger.info("Searching for taxonomy files in project output...")
        try:
            metadata_file = find_taxonomy_files(project_output, vis_output_dir)
            if metadata_file:
                logger.info(f"? Found and processed taxonomy file: {metadata_file}")
        except Exception as e:
            logger.warning(f"Error searching for taxonomy files: {e}")
    
    if not metadata_file:
        logger.warning("No metadata file found. Creating basic metadata from tree file...")
        try:
            metadata_file = create_metadata_from_tree(tree_files[0], vis_output_dir)
            if metadata_file:
                logger.info(f"? Created basic metadata: {metadata_file}")
            else:
                logger.error("Failed to create metadata from tree file")
                return False
        except Exception as e:
            logger.error(f"Error creating metadata from tree: {e}")
            return False
    
    # Process each tree file
    all_success = True
    total_files = len(tree_files)
    
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"Processing {total_files} tree file(s) for visualization...")
    logger.info("Will create the following visualization types:")
    logger.info("  1. Simplified taxonomic (phylum & genus rings)")
    logger.info("  2. Rectangular tree with colored strips")
    if functional_data:
        logger.info(f"  3. Functional visualizations ({', '.join(functional_data.keys())})")
    logger.info("=" * 60)
    
    for i, tree_file in enumerate(tree_files, 1):
        tree_basename = os.path.basename(tree_file)
        output_prefix = os.path.join(vis_output_dir, os.path.splitext(tree_basename)[0])
        
        logger.info(f"\n[{i}/{total_files}] Processing tree: {tree_basename}")
        logger.info("-" * 40)
        
        # Create simplified taxonomic visualization
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
        
        # Create rectangular tree visualization - USING FIXED FUNCTION
        try:
            logger.info("Creating rectangular tree visualization...")
            success = create_rectangular_tree_visualization_publication(
                tree_file, metadata_file, output_prefix, novel_mags
            )
            if success:
                logger.info("? Successfully created rectangular tree visualization")
            else:
                logger.warning("? Failed to create rectangular tree visualization")
                all_success = False
        except Exception as e:
            logger.error(f"? Error creating rectangular tree visualization: {e}")
            import traceback
            logger.error(traceback.format_exc())
            all_success = False
        
        # Create functional visualizations if data is available
        if functional_data:
            # Check if pie module is available
            try:
                from MetaMAG import enhanced_phylo_tree_vis_pie as pie_viz
                PIE_MODULE_AVAILABLE = True
            except ImportError:
                PIE_MODULE_AVAILABLE = False
                logger.warning("Pie module not available, skipping functional visualizations")
            
            if PIE_MODULE_AVAILABLE:
                try:
                    if len(functional_data) > 1:
                        # Multi-functional visualization with all annotation types
                        logger.info(f"Creating multi-functional visualization with {len(functional_data)} annotation types...")
                        logger.info(f"Annotation types: {list(functional_data.keys())}")
                        
                        success = pie_viz.create_multi_functional_visualization(
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
                        
                        func_type = list(functional_data.keys())[0]
                        single_functional_data = functional_data[func_type]
                        
                        success = pie_viz.create_publication_quality_visualization(
                            tree_file, metadata_file, output_prefix, 
                            novel_mags, single_functional_data, func_type
                        )
                        
                        if success:
                            logger.info("? Successfully created single functional visualization")
                        else:
                            logger.warning("? Failed to create single functional visualization")
                            all_success = False
                            
                except Exception as e:
                    logger.error(f"? Error creating functional visualization: {e}")
                    import traceback
                    logger.error(traceback.format_exc())
                    all_success = False
        
        # List all created files for this tree
        logger.info(f"\nCreated files for {tree_basename}:")
        potential_files = [
            f"{output_prefix}_simplified_taxonomic.pdf",
            f"{output_prefix}_rectangular_tree.pdf",
            f"{output_prefix}_multi_functional.pdf",
            f"{output_prefix}_single_functional.pdf",
            f"{output_prefix}_enhanced_beautiful.pdf",
            f"{output_prefix}_enhanced_publication.pdf"
        ]
        
        created_count = 0
        for potential_file in potential_files:
            if os.path.exists(potential_file):
                file_size = os.path.getsize(potential_file) / 1024  # KB
                logger.info(f"  ? {os.path.basename(potential_file)} ({file_size:.1f} KB)")
                created_count += 1
        
        if created_count == 0:
            logger.warning("  ? No visualization files were created!")
    
    # Final summary
    end_time = time.time()
    duration = end_time - start_time
    
    logger.info("")
    logger.info("=" * 60)
    if all_success:
        logger.info("? ENHANCED TREE VISUALIZATION COMPLETED SUCCESSFULLY")
    else:
        logger.warning("? ENHANCED TREE VISUALIZATION COMPLETED WITH SOME ERRORS")
    
    logger.info(f"Total processing time: {duration:.2f} seconds ({duration/60:.1f} minutes)")
    logger.info(f"Output directory: {vis_output_dir}")
    
    return all_success

def create_rectangular_tree_visualization_publication(tree_file, metadata_file, output_prefix, novel_mags):
    """FIXED: Create rectangular tree with PUBLICATION-QUALITY fonts and proper error handling."""
    try:
        script_path = f"{output_prefix}_rectangular_vis.R"
        
        logger.info(f"Creating PUBLICATION-QUALITY rectangular tree visualization")
        logger.info(f"  Tree file: {tree_file}")
        logger.info(f"  Metadata file: {metadata_file}")
        logger.info(f"  Output prefix: {output_prefix}")
        logger.info(f"  Novel MAGs: {len(novel_mags)}")
        
        # Prepare data for R script with better escaping
        novel_mags_json = json.dumps(novel_mags).replace("'", "\\'")
        
        r_script = f'''#!/usr/bin/env Rscript

# PUBLICATION-QUALITY Rectangular Tree Visualization with Enhanced Fonts
# All fonts sized for readability in publications and Word documents

# Enhanced error handling
options(error = function() {{
  traceback()
  quit(status = 1)
}})

# Load required packages with error handling
required_packages <- c("ape", "RColorBrewer", "jsonlite")
for (pkg in required_packages) {{
  if (!requireNamespace(pkg, quietly = TRUE)) {{
    cat("ERROR: Package", pkg, "is not installed\\n")
    quit(status = 1)
  }}
  library(pkg, character.only = TRUE)
  cat("Loaded package:", pkg, "\\n")
}}

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- tryCatch({{
  fromJSON('{novel_mags_json}')
}}, error = function(e) {{
  cat("Warning: Could not parse novel MAGs JSON\\n")
  character(0)
}})

set.seed(42)

cat("\\n=== Starting PUBLICATION-QUALITY Rectangular Tree Visualization ===\\n")
cat("Tree file exists:", file.exists(tree_file), "\\n")
cat("Metadata file exists:", file.exists(metadata_file), "\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")

# Read and process tree with better error handling
tree <- NULL
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded successfully with", length(tree$tip.label), "tips\\n")
  
  # Ladderize for better visualization
  tree <- ladderize(tree, right = FALSE)
  cat("Tree ladderized\\n")
  
  # Root if needed
  if (!is.rooted(tree)) {{
    if (length(tree$tip.label) > 2) {{
      outgroup_tip <- tree$tip.label[1]
      tryCatch({{
        tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
        tree <- ladderize(tree, right = FALSE)
        cat("Tree rooted using tip:", outgroup_tip, "\\n")
      }}, error = function(e) {{
        cat("Warning: Could not root tree, using as-is\\n")
      }})
    }}
  }}
}}, error = function(e) {{
  cat("ERROR loading tree:", e$message, "\\n")
  quit(status = 1)
}})

# Verify tree object
if (is.null(tree)) {{
  cat("ERROR: Tree object is NULL\\n")
  quit(status = 1)
}}

# Read metadata with error handling
metadata <- NULL
tryCatch({{
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns\\n")
  cat("Metadata columns:", paste(colnames(metadata), collapse=", "), "\\n")
}}, error = function(e) {{
  cat("ERROR loading metadata:", e$message, "\\n")
  quit(status = 1)
}})

# Get ID column
id_col <- colnames(metadata)[1]
cat("Using ID column:", id_col, "\\n")

# Process taxonomy data
taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]], stringsAsFactors = FALSE)

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
  # Try exact match
  match_idx <- which(taxa_data$ID == tip)
  
  # Try case-insensitive
  if (length(match_idx) == 0) {{
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  
  # Try with underscore/dot conversion
  if (length(match_idx) == 0) {{
    converted_tip <- gsub("\\\\.", "_", tip)
    match_idx <- which(taxa_data$ID == converted_tip)
  }}
  
  # Try partial matching
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

# Create ordered tip data
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

# Create color palettes - Enhanced for publication quality
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

cat("Unique phyla:", length(phyla), "\\n")
cat("Unique genera:", length(genera), "\\n")

# Enhanced color palettes for publication quality
if (length(phyla) > 0) {{
  if (length(phyla) <= 8) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Set1")[1:length(phyla)]
  }} else if (length(phyla) <= 12) {{
    phylum_colors <- c(brewer.pal(8, "Set1"), brewer.pal(max(3, min(4, length(phyla)-8)), "Dark2"))
  }} else {{
    # Use colorRampPalette for many phyla
    base_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
    phylum_colors <- colorRampPalette(base_colors)(length(phyla))
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
    # Use more distinct colors for many genera
    base_colors <- c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))
    genus_colors <- colorRampPalette(base_colors)(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "#FAFAFA")
}} else {{
  genus_colors <- c("Unknown" = "#FAFAFA")
}}

##############################################
# CREATE PUBLICATION-QUALITY RECTANGULAR TREE
##############################################

rectangular_file <- paste0(output_prefix, "_rectangular_tree.pdf")

# Calculate optimal figure height based on number of tips
fig_height <- max(14, min(50, length(tree$tip.label) * 0.35))
fig_width <- 24

cat("Creating PDF with dimensions:", fig_width, "x", fig_height, "\\n")

tryCatch({{
  pdf(rectangular_file, width = fig_width, height = fig_height)
  
  # Enhanced layout for publication quality
  layout(matrix(c(1,2), nrow=1), widths=c(4, 1.2))
  
  # Main tree plot with adjusted margins
  par(mar=c(4, 1, 5, 14), xpd=TRUE, cex=1.0)
  
  # Plot tree with LARGER tip labels
  plot(tree, type="phylogram", show.tip.label=TRUE, 
       main="", 
       edge.width=3, 
       cex.main=2.5, 
       edge.color="#2C3E50",
       x.lim=c(0, max(node.depth.edgelength(tree)) * 1.45),
       label.offset=0.002,
       cex=1.2,  # INCREASED tip label size
       font=1)   # Regular font for tips
  
  # Add title separately for better control
  title(main="Phylogenetic Tree with Taxonomic Classification", 
        cex.main=3.0, font.main=2, line=3)
  
  # Get tree layout information
  last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tips <- length(tree$tip.label)
  
  cat("Tree plotted with", n_tips, "tips\\n")
  
  # Calculate positions for colored strips
  x_tree_max <- max(last_plot$xx[1:n_tips])
  strip_start_x <- x_tree_max + 0.03 * max(last_plot$xx)
  strip_width <- 0.035 * max(last_plot$xx)
  strip_gap <- 0.01 * max(last_plot$xx)
  
  # Add column headers with LARGER FONTS
  text(strip_start_x + strip_width/2, n_tips + 2, "Phylum", 
       font = 2, cex = 2.0, srt = 0)
  text(strip_start_x + strip_width + strip_gap + strip_width/2, n_tips + 2, "Genus", 
       font = 2, cex = 2.0, srt = 0)
  
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
      
      # Add thin border for clarity
      rect(strip_start_x, y_pos - 0.4, 
           strip_start_x + strip_width, y_pos + 0.4,
           col = NA, border = "gray70", lwd = 0.8)
      
      # Genus strip
      genus_color <- if (genus != "Unknown" && genus %in% names(genus_colors)) {{
        genus_colors[genus]
      }} else {{
        genus_colors["Unknown"]
      }}
      
      rect(strip_start_x + strip_width + strip_gap, y_pos - 0.4, 
           strip_start_x + 2*strip_width + strip_gap, y_pos + 0.4,
           col = genus_color, border = NA)
      
      # Add thin border for clarity
      rect(strip_start_x + strip_width + strip_gap, y_pos - 0.4, 
           strip_start_x + 2*strip_width + strip_gap, y_pos + 0.4,
           col = NA, border = "gray70", lwd = 0.8)
      
      # Novel MAG indicator - make more prominent
      if (is_novel) {{
        # Add thick red border around both strips
        rect(strip_start_x - 0.002 * max(last_plot$xx), y_pos - 0.48, 
             strip_start_x + 2*strip_width + strip_gap + 0.002 * max(last_plot$xx), y_pos + 0.48,
             col = NA, border = "#FF0000", lwd = 4)
        
        # Add star symbol
        points(strip_start_x + 2*strip_width + strip_gap + 0.025 * max(last_plot$xx), 
               y_pos, pch = 8, col = "#FF0000", cex = 3, lwd = 3)
      }}
    }}
  }}
  
  # ENHANCED LEGEND - PUBLICATION QUALITY WITH LARGE FONTS
  par(mar=c(5, 2, 5, 2), xpd=TRUE, cex=1.0)
  plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
  
  # Title with MUCH larger font
  text(0.5, 0.98, "Legend", font=2, cex=3.0, adj=0.5)
  
  # Phylum legend - enhanced size
  text(0.05, 0.90, "Phylum", font=2, pos=4, cex=2.2)
  legend_y <- 0.84
  n_phyla_shown <- min(15, length(phyla))  # Show more items
  
  for (i in 1:n_phyla_shown) {{
    # Larger color boxes
    rect(0.05, legend_y - 0.022, 0.20, legend_y + 0.022, 
         col=phylum_colors[phyla[i]], border="black", lwd=1.5)
    
    # Truncate long names
    phylum_name <- phyla[i]
    if (nchar(phylum_name) > 16) {{
      phylum_name <- paste0(substr(phylum_name, 1, 14), "...")
    }}
    text(0.22, legend_y, phylum_name, pos=4, cex=1.6)  # MUCH larger text
    legend_y <- legend_y - 0.045
  }}
  
  if (length(phyla) > n_phyla_shown) {{
    text(0.22, legend_y, paste0("... +", length(phyla) - n_phyla_shown, " more"), 
         pos=4, cex=1.4, font=3)
    legend_y <- legend_y - 0.045
  }}
  
  # Genus legend - enhanced size
  legend_y <- legend_y - 0.05
  text(0.05, legend_y, "Genus", font=2, pos=4, cex=2.2)
  legend_y <- legend_y - 0.045
  n_genera_shown <- min(10, length(genera))
  
  for (i in 1:n_genera_shown) {{
    # Larger color boxes
    rect(0.05, legend_y - 0.022, 0.20, legend_y + 0.022, 
         col=genus_colors[genera[i]], border="black", lwd=1.5)
    
    # Truncate long names
    genus_name <- genera[i]
    if (nchar(genus_name) > 16) {{
      genus_name <- paste0(substr(genus_name, 1, 14), "...")
    }}
    text(0.22, legend_y, genus_name, pos=4, cex=1.6)  # MUCH larger text
    legend_y <- legend_y - 0.045
  }}
  
  if (length(genera) > n_genera_shown) {{
    text(0.22, legend_y, paste0("... +", length(genera) - n_genera_shown, " more"), 
         pos=4, cex=1.4, font=3)
    legend_y <- legend_y - 0.045
  }}
  
  # Novel MAG indicator - enhanced
  legend_y <- legend_y - 0.06
  text(0.05, legend_y, "Novel MAGs", font=2, pos=4, cex=2.2)
  legend_y <- legend_y - 0.045
  rect(0.05, legend_y - 0.022, 0.20, legend_y + 0.022, 
       col=NA, border="#FF0000", lwd=4)
  points(0.125, legend_y, pch=8, col="#FF0000", cex=3, lwd=3)
  text(0.22, legend_y, paste0("Novel (n=", novel_count, ")"), pos=4, cex=1.6, col="#FF0000")
  
  # Statistics - enhanced
  legend_y <- legend_y - 0.08
  text(0.05, legend_y, "Statistics", font=2, pos=4, cex=2.2)
  legend_y <- legend_y - 0.045
  text(0.05, legend_y, paste0("Total MAGs: ", n_tips), pos=4, cex=1.6)
  legend_y <- legend_y - 0.035
  text(0.05, legend_y, paste0("Phyla: ", length(phyla)), pos=4, cex=1.6)
  legend_y <- legend_y - 0.035
  text(0.05, legend_y, paste0("Genera: ", length(genera)), pos=4, cex=1.6)
  
  dev.off()
  cat("PDF created successfully:", rectangular_file, "\\n")
  
}}, error = function(e) {{
  cat("ERROR creating PDF:", e$message, "\\n")
  if (dev.cur() > 1) dev.off()
  quit(status = 1)
}})

# Verify file was created
if (file.exists(rectangular_file)) {{
  file_size <- file.info(rectangular_file)$size / 1024
  cat("\\n=== Rectangular Tree Visualization Complete ===\\n")
  cat("? Created:", rectangular_file, "\\n")
  cat("? File size:", round(file_size, 1), "KB\\n")
  cat(sprintf("? Summary: %d total MAGs, %d novel MAGs, %d phyla, %d genera\\n", 
             n_tips, novel_count, length(phyla), length(genera)))
}} else {{
  cat("ERROR: PDF file was not created\\n")
  quit(status = 1)
}}

quit(status = 0)
'''
        
        # Write R script
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        os.chmod(script_path, 0o755)
        logger.info(f"R script written to: {script_path}")
        
        # Run the R script with better error capture
        success = run_r_command_with_enhanced_error_handling(script_path)
        
        if success:
            output_file = f"{output_prefix}_rectangular_tree.pdf"
            if os.path.exists(output_file):
                logger.info(f"? Successfully created rectangular tree: {output_file}")
                file_size = os.path.getsize(output_file) / 1024
                logger.info(f"? File size: {file_size:.1f} KB")
                return True
            else:
                logger.error(f"R script ran but output file not created: {output_file}")
                return False
        else:
            logger.error("R script execution failed")
            return False
        
    except Exception as e:
        logger.error(f"Error creating rectangular tree visualization: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False


def check_r_packages():
    """Check if required R packages are installed."""
    import subprocess
    
    required_packages = ['ape', 'RColorBrewer', 'jsonlite', 'ggtree', 'ggplot2']
    
    r_check_script = """
    packages <- c('%s')
    for (pkg in packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            cat(sprintf("? %%s installed\\n", pkg))
        } else {
            cat(sprintf("? %%s NOT installed\\n", pkg))
        }
    }
    """ % "', '".join(required_packages)
    
    try:
        result = subprocess.run(
            ['Rscript', '-e', r_check_script],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        logger.info("R Package Status:")
        for line in result.stdout.split('\n'):
            if line.strip():
                logger.info(f"  {line}")
                
        return result.returncode == 0
    except Exception as e:
        logger.error(f"Error checking R packages: {e}")
        return False


# Function to validate annotation data format
def validate_functional_data(functional_data):
    """Validate that functional data has the correct format."""
    
    for data_type, mags_data in functional_data.items():
        logger.info(f"Validating {data_type} data...")
        
        if not isinstance(mags_data, dict):
            logger.error(f"  ? {data_type} data is not a dictionary")
            return False
        
        if len(mags_data) == 0:
            logger.warning(f"  ? {data_type} data is empty")
            continue
        
        # Check a sample MAG
        sample_mag = list(mags_data.keys())[0]
        sample_data = mags_data[sample_mag]
        
        if 'categories' not in sample_data:
            logger.error(f"  ? {data_type} data missing 'categories' field")
            return False
        
        if not isinstance(sample_data['categories'], dict):
            logger.error(f"  ? {data_type} categories is not a dictionary")
            return False
        
        logger.info(f"  ? {data_type} data format valid ({len(mags_data)} MAGs)")
    
    return True

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
        return load_eggnog_annotations(functional_annotations, project_output)
    elif annotation_type == "kegg":
        return load_kegg_annotations(functional_annotations, project_output)
    else:
        logger.warning(f"Unknown annotation type: {annotation_type}")
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

def load_kegg_annotations_fixed(kegg_file, tree_tips):
    """FIXED: Load KEGG annotations with better pathway handling."""
    try:
        logger.info(f"Reading KEGG Excel file: {kegg_file}")
        df = pd.read_excel(kegg_file)
        logger.info(f"KEGG file shape: {df.shape}")
        logger.info(f"KEGG columns: {list(df.columns)}")
        
        # Check file structure
        if df.shape[1] < 2:
            logger.error("KEGG file has too few columns")
            return {}
        
        # Identify MAG columns (all except first)
        mag_columns = df.columns[1:].tolist()
        logger.info(f"Found {len(mag_columns)} MAG columns")
        
        functional_data = {}
        
        for mag_col in mag_columns:
            # Clean MAG name
            mag_name = clean_mag_column_name(mag_col)
            if not mag_name:
                logger.debug(f"Skipping invalid MAG column: {mag_col}")
                continue
            
            # Extract KEGG pathway counts
            kegg_counts = {}
            
            for idx, row in df.iterrows():
                try:
                    # Get pathway name from first column
                    pathway = str(row.iloc[0]).strip() if not pd.isna(row.iloc[0]) else ""
                    
                    # Get count for this MAG
                    count = row[mag_col]
                    
                    if pd.isna(count) or count == 0:
                        continue
                    
                    count = int(float(count))
                    
                    if count > 0 and pathway and pathway != 'nan':
                        # Enhanced pathway name simplification
                        pathway_short = simplify_kegg_pathway_name_enhanced(pathway)
                        if pathway_short:
                            if pathway_short in kegg_counts:
                                kegg_counts[pathway_short] += count
                            else:
                                kegg_counts[pathway_short] = count
                        
                except (ValueError, TypeError, KeyError) as e:
                    logger.debug(f"Error processing KEGG row {idx}: {e}")
                    continue
            
            # Add to functional data if we have pathways
            if kegg_counts:
                # Sort by count and take top categories for visualization
                sorted_kegg = dict(sorted(kegg_counts.items(), key=lambda x: x[1], reverse=True))
                
                # Limit to top 10 pathways for cleaner visualization
                top_kegg = dict(list(sorted_kegg.items())[:10])
                
                functional_data[mag_name] = {
                    'categories': top_kegg,
                    'total_genes': sum(kegg_counts.values()),
                    'total_pathways': len(kegg_counts),
                    'annotation_type': 'kegg'
                }
                
                logger.debug(f"MAG {mag_name}: {len(kegg_counts)} pathways, {sum(kegg_counts.values())} total genes")
        
        logger.info(f"Processed {len(functional_data)} MAGs with KEGG data")
        
        # Match to tree tips
        matched_data = match_functional_data_to_tree_tips(functional_data, tree_tips)
        logger.info(f"Matched {len(matched_data)} MAGs to tree tips")
        
        return matched_data
        
    except Exception as e:
        logger.error(f"Error loading KEGG annotations: {e}")
        import traceback
        logger.error(traceback.format_exc())
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
    """Enhanced MAG column name cleaning."""
    if not mag_col or mag_col in ['nan', 'NaN', '', 'Unnamed: 0']:
        return None
    
    mag_name = str(mag_col)
    original_name = mag_name
    
    # Remove various suffixes in order
    suffixes_to_remove = [
        '_fa_faa', '_fa.counts', '.counts', '_fa', '.fa', 
        '.fasta', '.fna', '_faa', '.faa', '_contigs', '_scaffolds'
    ]
    
    for suffix in suffixes_to_remove:
        if mag_name.endswith(suffix):
            mag_name = mag_name[:-len(suffix)]
    
    # Convert bin_NUMBER to bin.NUMBER
    mag_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
    
    clean_name = mag_name.strip()
    
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

def create_rectangular_tree_visualization_fixed(tree_file, metadata_file, output_prefix, novel_mags):
    """FIXED: Create rectangular tree with colored strips - enhanced error handling and publication quality."""
    try:
        script_path = f"{output_prefix}_rectangular_vis.R"
        
        logger.info(f"Creating rectangular tree visualization")
        logger.info(f"  Tree file: {tree_file}")
        logger.info(f"  Metadata file: {metadata_file}")
        logger.info(f"  Output prefix: {output_prefix}")
        logger.info(f"  Novel MAGs: {len(novel_mags)}")
        
        # Prepare data for R script with better escaping
        novel_mags_json = json.dumps(novel_mags).replace("'", "\\'")
        
        r_script = f'''#!/usr/bin/env Rscript

# Rectangular Tree Visualization with Taxonomic Strips - PUBLICATION QUALITY

# Enhanced error handling
options(error = function() {{
  traceback()
  quit(status = 1)
}})

# Load required packages with error handling
required_packages <- c("ape", "RColorBrewer", "jsonlite")
for (pkg in required_packages) {{
  if (!requireNamespace(pkg, quietly = TRUE)) {{
    cat("ERROR: Package", pkg, "is not installed\\n")
    quit(status = 1)
  }}
  library(pkg, character.only = TRUE)
  cat("Loaded package:", pkg, "\\n")
}}

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- tryCatch({{
  fromJSON('{novel_mags_json}')
}}, error = function(e) {{
  cat("Warning: Could not parse novel MAGs JSON\\n")
  character(0)
}})

set.seed(42)

cat("\\n=== Starting Rectangular Tree Visualization ===\\n")
cat("Tree file exists:", file.exists(tree_file), "\\n")
cat("Metadata file exists:", file.exists(metadata_file), "\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")

# Read and process tree with better error handling
tree <- NULL
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded successfully with", length(tree$tip.label), "tips\\n")
  
  # Ladderize for better visualization
  tree <- ladderize(tree, right = FALSE)
  cat("Tree ladderized\\n")
  
  # Root if needed
  if (!is.rooted(tree)) {{
    if (length(tree$tip.label) > 2) {{
      outgroup_tip <- tree$tip.label[1]
      tryCatch({{
        tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
        tree <- ladderize(tree, right = FALSE)
        cat("Tree rooted using tip:", outgroup_tip, "\\n")
      }}, error = function(e) {{
        cat("Warning: Could not root tree, using as-is\\n")
      }})
    }}
  }}
}}, error = function(e) {{
  cat("ERROR loading tree:", e$message, "\\n")
  quit(status = 1)
}})

# Verify tree object
if (is.null(tree)) {{
  cat("ERROR: Tree object is NULL\\n")
  quit(status = 1)
}}

# Read metadata with error handling
metadata <- NULL
tryCatch({{
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns\\n")
  cat("Metadata columns:", paste(colnames(metadata), collapse=", "), "\\n")
}}, error = function(e) {{
  cat("ERROR loading metadata:", e$message, "\\n")
  quit(status = 1)
}})

# Get ID column
id_col <- colnames(metadata)[1]
cat("Using ID column:", id_col, "\\n")

# Process taxonomy data
taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]], stringsAsFactors = FALSE)

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
  # Try exact match
  match_idx <- which(taxa_data$ID == tip)
  
  # Try case-insensitive
  if (length(match_idx) == 0) {{
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  
  # Try with underscore/dot conversion
  if (length(match_idx) == 0) {{
    converted_tip <- gsub("\\\\.", "_", tip)
    match_idx <- which(taxa_data$ID == converted_tip)
  }}
  
  # Try partial matching
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

# Create ordered tip data
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

# Create color palettes - Enhanced for publication quality
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

cat("Unique phyla:", length(phyla), "\\n")
cat("Unique genera:", length(genera), "\\n")

# Enhanced color palettes for publication quality
if (length(phyla) > 0) {{
  if (length(phyla) <= 8) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Set1")[1:length(phyla)]
  }} else if (length(phyla) <= 12) {{
    phylum_colors <- c(brewer.pal(8, "Set1"), brewer.pal(max(3, min(4, length(phyla)-8)), "Dark2"))
  }} else {{
    # Use colorRampPalette for many phyla
    base_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
    phylum_colors <- colorRampPalette(base_colors)(length(phyla))
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
    # Use more distinct colors for many genera
    base_colors <- c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))
    genus_colors <- colorRampPalette(base_colors)(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "#FAFAFA")
}} else {{
  genus_colors <- c("Unknown" = "#FAFAFA")
}}

##############################################
# CREATE PUBLICATION-QUALITY RECTANGULAR TREE
##############################################

rectangular_file <- paste0(output_prefix, "_rectangular_tree.pdf")

# Calculate optimal figure height based on number of tips
fig_height <- max(14, min(50, length(tree$tip.label) * 0.3))
fig_width <- 20

cat("Creating PDF with dimensions:", fig_width, "x", fig_height, "\\n")

tryCatch({{
  pdf(rectangular_file, width = fig_width, height = fig_height)
  
  # Enhanced layout for publication quality
  layout(matrix(c(1,2), nrow=1), widths=c(4.5, 1))
  
  # Main tree plot with adjusted margins
  par(mar=c(4, 1, 4, 12), xpd=TRUE, cex=1.0)
  
  # Plot tree
  plot(tree, type="phylogram", show.tip.label=TRUE, 
       main="", 
       edge.width=2.5, 
       cex.main=2.0, 
       edge.color="#2C3E50",
       x.lim=c(0, max(node.depth.edgelength(tree)) * 1.4),
       label.offset=0.002,
       cex=0.8,  # Tip label size
       font=1)   # Regular font for tips
  
  # Add title separately for better control
  title(main="Phylogenetic Tree with Taxonomic Classification", 
        cex.main=2.2, font.main=2, line=2)
  
  # Get tree layout information
  last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tips <- length(tree$tip.label)
  
  cat("Tree plotted with", n_tips, "tips\\n")
  
  # Calculate positions for colored strips
  x_tree_max <- max(last_plot$xx[1:n_tips])
  strip_start_x <- x_tree_max + 0.03 * max(last_plot$xx)
  strip_width <- 0.03 * max(last_plot$xx)
  strip_gap <- 0.008 * max(last_plot$xx)
  
  # Add column headers
  text(strip_start_x + strip_width/2, n_tips + 1.8, "Phylum", 
       font = 2, cex = 1.4, srt = 0)
  text(strip_start_x + strip_width + strip_gap + strip_width/2, n_tips + 1.8, "Genus", 
       font = 2, cex = 1.4, srt = 0)
  
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
      
      rect(strip_start_x, y_pos - 0.35, 
           strip_start_x + strip_width, y_pos + 0.35,
           col = phylum_color, border = NA)
      
      # Add thin border for clarity
      rect(strip_start_x, y_pos - 0.35, 
           strip_start_x + strip_width, y_pos + 0.35,
           col = NA, border = "gray80", lwd = 0.5)
      
      # Genus strip
      genus_color <- if (genus != "Unknown" && genus %in% names(genus_colors)) {{
        genus_colors[genus]
      }} else {{
        genus_colors["Unknown"]
      }}
      
      rect(strip_start_x + strip_width + strip_gap, y_pos - 0.35, 
           strip_start_x + 2*strip_width + strip_gap, y_pos + 0.35,
           col = genus_color, border = NA)
      
      # Add thin border for clarity
      rect(strip_start_x + strip_width + strip_gap, y_pos - 0.35, 
           strip_start_x + 2*strip_width + strip_gap, y_pos + 0.35,
           col = NA, border = "gray80", lwd = 0.5)
      
      # Novel MAG indicator - make more prominent
      if (is_novel) {{
        # Add thick red border around both strips
        rect(strip_start_x - 0.002 * max(last_plot$xx), y_pos - 0.42, 
             strip_start_x + 2*strip_width + strip_gap + 0.002 * max(last_plot$xx), y_pos + 0.42,
             col = NA, border = "#FF0000", lwd = 3)
        
        # Add star symbol
        points(strip_start_x + 2*strip_width + strip_gap + 0.02 * max(last_plot$xx), 
               y_pos, pch = 8, col = "#FF0000", cex = 2, lwd = 2.5)
      }}
    }}
  }}
  
  # ENHANCED LEGEND - PUBLICATION QUALITY
  par(mar=c(5, 2, 5, 2), xpd=TRUE, cex=1.0)
  plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
  
  # Title with larger font
  text(0.5, 0.98, "Legend", font=2, cex=2.0, adj=0.5)
  
  # Phylum legend - enhanced size
  text(0.05, 0.91, "Phylum", font=2, pos=4, cex=1.6)
  legend_y <- 0.85
  n_phyla_shown <- min(20, length(phyla))  # Show more items
  
  for (i in 1:n_phyla_shown) {{
    # Larger color boxes
    rect(0.05, legend_y - 0.018, 0.18, legend_y + 0.018, 
         col=phylum_colors[phyla[i]], border="black", lwd=1.2)
    
    # Truncate long names
    phylum_name <- phyla[i]
    if (nchar(phylum_name) > 18) {{
      phylum_name <- paste0(substr(phylum_name, 1, 16), "...")
    }}
    text(0.20, legend_y, phylum_name, pos=4, cex=1.1)  # Larger text
    legend_y <- legend_y - 0.035
  }}
  
  if (length(phyla) > n_phyla_shown) {{
    text(0.20, legend_y, paste0("... +", length(phyla) - n_phyla_shown, " more"), 
         pos=4, cex=1.0, font=3)
    legend_y <- legend_y - 0.035
  }}
  
  # Genus legend - enhanced size
  legend_y <- legend_y - 0.04
  text(0.05, legend_y, "Genus", font=2, pos=4, cex=1.6)
  legend_y <- legend_y - 0.038
  n_genera_shown <- min(12, length(genera))
  
  for (i in 1:n_genera_shown) {{
    # Larger color boxes
    rect(0.05, legend_y - 0.018, 0.18, legend_y + 0.018, 
         col=genus_colors[genera[i]], border="black", lwd=1.2)
    
    # Truncate long names
    genus_name <- genera[i]
    if (nchar(genus_name) > 18) {{
      genus_name <- paste0(substr(genus_name, 1, 16), "...")
    }}
    text(0.20, legend_y, genus_name, pos=4, cex=1.1)  # Larger text
    legend_y <- legend_y - 0.035
  }}
  
  if (length(genera) > n_genera_shown) {{
    text(0.20, legend_y, paste0("... +", length(genera) - n_genera_shown, " more"), 
         pos=4, cex=1.0, font=3)
    legend_y <- legend_y - 0.035
  }}
  
  # Novel MAG indicator - enhanced
  legend_y <- legend_y - 0.05
  text(0.05, legend_y, "Novel MAGs", font=2, pos=4, cex=1.6)
  legend_y <- legend_y - 0.038
  rect(0.05, legend_y - 0.018, 0.18, legend_y + 0.018, 
       col=NA, border="#FF0000", lwd=3)
  points(0.115, legend_y, pch=8, col="#FF0000", cex=2, lwd=2.5)
  text(0.20, legend_y, paste0("Novel (n=", novel_count, ")"), pos=4, cex=1.1)
  
  # Statistics - enhanced
  legend_y <- legend_y - 0.06
  text(0.05, legend_y, "Statistics", font=2, pos=4, cex=1.6)
  legend_y <- legend_y - 0.038
  text(0.05, legend_y, paste0("Total MAGs: ", n_tips), pos=4, cex=1.1)
  legend_y <- legend_y - 0.03
  text(0.05, legend_y, paste0("Phyla: ", length(phyla)), pos=4, cex=1.1)
  legend_y <- legend_y - 0.03
  text(0.05, legend_y, paste0("Genera: ", length(genera)), pos=4, cex=1.1)
  
  dev.off()
  cat("PDF created successfully:", rectangular_file, "\\n")
  
}}, error = function(e) {{
  cat("ERROR creating PDF:", e$message, "\\n")
  if (dev.cur() > 1) dev.off()
  quit(status = 1)
}})

# Verify file was created
if (file.exists(rectangular_file)) {{
  file_size <- file.info(rectangular_file)$size / 1024
  cat("\\n=== Rectangular Tree Visualization Complete ===\\n")
  cat("Created:", rectangular_file, "\\n")
  cat("File size:", round(file_size, 1), "KB\\n")
  cat(sprintf("Summary: %d total MAGs, %d novel MAGs, %d phyla, %d genera\\n", 
             n_tips, novel_count, length(phyla), length(genera)))
}} else {{
  cat("ERROR: PDF file was not created\\n")
  quit(status = 1)
}}

quit(status = 0)
'''
        
        # Write R script
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        os.chmod(script_path, 0o755)
        logger.info(f"R script written to: {script_path}")
        
        # Run the R script with better error capture
        success = run_r_command_with_enhanced_error_handling(script_path)
        
        if success:
            output_file = f"{output_prefix}_rectangular_tree.pdf"
            if os.path.exists(output_file):
                logger.info(f"? Successfully created rectangular tree: {output_file}")
                return True
            else:
                logger.error(f"R script ran but output file not created: {output_file}")
                return False
        else:
            logger.error("R script execution failed")
            return False
        
    except Exception as e:
        logger.error(f"Error creating rectangular tree visualization: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def run_r_command_with_enhanced_error_handling(script_path):
    """Run R script with enhanced error handling and logging."""
    try:
        # Try multiple R execution methods
        commands = [
            f"Rscript {script_path}",
            f"R CMD BATCH --no-save --no-restore {script_path}",
            f"/usr/local/bin/Rscript {script_path}",
            f"/usr/bin/Rscript {script_path}"
        ]
        
        for cmd in commands:
            logger.info(f"Attempting to run R with command: {cmd}")
            
            result = subprocess.run(
                cmd, 
                shell=True, 
                check=False,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=300  # 5 minute timeout
            )
            
            # Log all output
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if line.strip():
                        if "ERROR" in line:
                            logger.error(f"R ERROR: {line}")
                        elif "Warning" in line:
                            logger.warning(f"R WARNING: {line}")
                        else:
                            logger.info(f"R OUTPUT: {line}")
            
            if result.stderr:
                for line in result.stderr.split('\n'):
                    if line.strip():
                        if "Warning" in line:
                            logger.warning(f"R WARNING: {line}")
                        else:
                            logger.error(f"R ERROR: {line}")
            
            # Check if successful
            if result.returncode == 0:
                logger.info("? R script executed successfully")
                return True
            else:
                logger.warning(f"R script failed with return code: {result.returncode}")
                # Check if output file was created anyway
                if "rectangular_tree.pdf" in script_path:
                    output_file = script_path.replace("_vis.R", ".pdf")
                    if os.path.exists(output_file):
                        logger.info("? Output file created despite non-zero return code")
                        return True
                continue  # Try next command
        
        # All commands failed
        logger.error("All R execution methods failed")
        return False
        
    except subprocess.TimeoutExpired:
        logger.error("R script execution timed out after 5 minutes")
        return False
    except Exception as e:
        logger.error(f"Error running R command: {e}")
        return False

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

def simplify_kegg_pathway_name_enhanced(pathway):
    """Enhanced KEGG pathway name simplification for better visualization."""
    if not pathway or pathway == 'nan':
        return None
    
    # Remove numeric codes more comprehensively
    pathway = re.sub(r'\d{5}', '', pathway)  # Remove 5-digit codes
    pathway = re.sub(r'^\d+\s*', '', pathway)  # Remove leading numbers
    pathway = re.sub(r'\s*\[\d+\]', '', pathway)  # Remove bracketed numbers
    pathway = pathway.strip()
    
    # Enhanced mapping for major KEGG categories
    name_map = {
        # Main categories
        'Metabolism': 'Metabolism',
        'Genetic Information Processing': 'Genetic Info',
        'Environmental Information Processing': 'Environ Info',
        'Cellular Processes': 'Cellular',
        'Organismal Systems': 'Organismal',
        'Human Diseases': 'Diseases',
        'Brite Hierarchies': 'Brite',
        'Not Included in Pathway or Brite': 'Other',
        
        # Subcategories - Metabolism
        'Carbohydrate metabolism': 'Carbohydrate',
        'Energy metabolism': 'Energy',
        'Lipid metabolism': 'Lipid',
        'Nucleotide metabolism': 'Nucleotide',
        'Amino acid metabolism': 'Amino Acid',
        'Metabolism of other amino acids': 'Other AA',
        'Glycan biosynthesis and metabolism': 'Glycan',
        'Metabolism of cofactors and vitamins': 'Cofactors',
        'Metabolism of terpenoids and polyketides': 'Terpenoids',
        'Biosynthesis of other secondary metabolites': 'Secondary',
        'Xenobiotics biodegradation and metabolism': 'Xenobiotics',
        
        # Subcategories - Genetic Information
        'Transcription': 'Transcription',
        'Translation': 'Translation',
        'Folding, sorting and degradation': 'Protein Proc',
        'Replication and repair': 'DNA Repair',
        
        # Subcategories - Environmental Information
        'Membrane transport': 'Transport',
        'Signal transduction': 'Signaling',
        'Signaling molecules and interaction': 'Signal Mol',
        
        # Subcategories - Cellular Processes
        'Transport and catabolism': 'Transport',
        'Cell growth and death': 'Cell Cycle',
        'Cellular community': 'Community',
        'Cell motility': 'Motility',
    }
    
    # Try exact match first
    if pathway in name_map:
        return name_map[pathway]
    
    # Try partial matches
    for key, value in name_map.items():
        if key.lower() in pathway.lower():
            return value
    
    # If no match, truncate long names
    if len(pathway) > 15:
        return pathway[:12] + "..."
    
    return pathway if pathway else None

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
def clean_mag_column_name(mag_col):
    """Enhanced MAG column name cleaning."""
    if not mag_col or mag_col in ['nan', 'NaN', '', 'Unnamed: 0']:
        return None
    
    mag_name = str(mag_col)
    original_name = mag_name
    
    # Remove various suffixes in order
    suffixes_to_remove = [
        '_fa_faa', '_fa.counts', '.counts', '_fa', '.fa', 
        '.fasta', '.fna', '_faa', '.faa', '_contigs', '_scaffolds'
    ]
    
    for suffix in suffixes_to_remove:
        if mag_name.endswith(suffix):
            mag_name = mag_name[:-len(suffix)]
    
    # Convert bin_NUMBER to bin.NUMBER
    mag_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
    
    clean_name = mag_name.strip()
    
    if original_name != clean_name:
        logger.debug(f"Name conversion: {original_name} -> {clean_name}")
    
    return clean_name if clean_name else None


def match_functional_data_to_tree_tips(raw_functional_data, tree_tips):
    """Enhanced matching with better debugging."""
    matched_data = {}
    
    logger.info(f"Matching {len(raw_functional_data)} functional entries to {len(tree_tips)} tree tips")
    
    if len(raw_functional_data) > 0:
        logger.debug(f"Sample functional MAGs: {list(raw_functional_data.keys())[:5]}")
    if len(tree_tips) > 0:
        logger.debug(f"Sample tree tips: {tree_tips[:5]}")
    
    # Multiple matching strategies
    exact_matches = 0
    conversion_matches = 0
    fuzzy_matches = 0
    
    for mag_name, data in raw_functional_data.items():
        matched = False
        
        # Strategy 1: Direct exact match
        if mag_name in tree_tips:
            matched_data[mag_name] = data
            exact_matches += 1
            matched = True
            continue
        
        # Strategy 2: Convert underscores to dots
        converted_name = re.sub(r'_bin_(\d+)', r'_bin.\1', mag_name)
        if converted_name != mag_name and converted_name in tree_tips:
            matched_data[converted_name] = data
            conversion_matches += 1
            matched = True
            continue
        
        # Strategy 3: Try with _sub suffix
        if converted_name + '_sub' in tree_tips:
            matched_data[converted_name + '_sub'] = data
            conversion_matches += 1
            matched = True
            continue
        
        # Strategy 4: Fuzzy matching
        if not matched:
            normalized_mag = normalize_mag_name(mag_name)
            for tree_tip in tree_tips:
                if tree_tip not in matched_data:
                    normalized_tip = normalize_mag_name(tree_tip)
                    if normalized_mag == normalized_tip or \
                       (len(normalized_mag) > 6 and len(normalized_tip) > 6 and 
                        (normalized_mag in normalized_tip or normalized_tip in normalized_mag)):
                        matched_data[tree_tip] = data
                        fuzzy_matches += 1
                        matched = True
                        break
    
    logger.info(f"Matching results: {exact_matches} exact, {conversion_matches} conversions, {fuzzy_matches} fuzzy")
    logger.info(f"Total matched: {len(matched_data)} out of {len(raw_functional_data)}")
    
    return matched_data


def normalize_mag_name(name):
    """Normalize MAG name for matching."""
    if not name:
        return ""
    
    norm = str(name).lower()
    
    # Remove common suffixes
    for suffix in ['_fa', '_fasta', '.fa', '.fasta', '.counts', '_contigs', '_scaffolds']:
        norm = norm.replace(suffix, '')
    
    # Standardize separators
    norm = norm.replace('.', '_').replace('-', '_')
    
    # Remove multiple underscores
    while '__' in norm:
        norm = norm.replace('__', '_')
    
    return norm.strip('_')

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

def process_standard_cazy_data_enhanced(df):
    """Enhanced processing of standard CAZy data format."""
    # This would implement the standard format processing
    # For now, return empty since your data appears to be transposed
    logger.info("Standard format processing not implemented yet")
    return {}

def create_rectangular_tree_with_pies(tree_file, metadata_file, output_prefix, 
                                      novel_mags, functional_data):
    """Create rectangular tree with pie charts for functional annotations - FIXED VERSION."""
    logger.info("="*60)
    logger.info("ENTERING create_rectangular_tree_with_pies FUNCTION")
    logger.info(f"Tree file: {tree_file}")
    logger.info(f"Metadata file: {metadata_file}")
    logger.info(f"Output prefix: {output_prefix}")
    logger.info(f"Number of novel MAGs: {len(novel_mags)}")
    logger.info(f"Functional data types: {list(functional_data.keys())}")
    logger.info("="*60)

    try:
        script_path = f"{output_prefix}_rectangular_pies_vis.R"
        
        logger.info(f"Creating rectangular tree with functional pie charts")
        logger.info(f"Functional annotation types: {list(functional_data.keys())}")
        logger.info(f"Novel MAGs to highlight: {len(novel_mags)}")
        
        # Prepare data for R script
        novel_mags_json = json.dumps(novel_mags)
        functional_data_json = json.dumps(functional_data)
        
        r_script = f'''#!/usr/bin/env Rscript

# Rectangular Tree with Functional Pie Charts Visualization - FIXED VERSION
# Combines rectangular tree layout with CAZyme, COG, and KEGG pie charts

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
functional_data <- fromJSON('{functional_data_json}')

set.seed(42)

cat("=== Rectangular Tree with Functional Pie Charts ===\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")
for (func_type in names(functional_data)) {{
  cat(paste("Functional data for", toupper(func_type), ":", length(functional_data[[func_type]]), "MAGs\\n"))
}}

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

# Add novel MAG and functional data
tip_data$is_novel <- FALSE
tip_data$cazyme_data <- rep(list(NULL), nrow(tip_data))
tip_data$cog_data <- rep(list(NULL), nrow(tip_data))
tip_data$kegg_data <- rep(list(NULL), nrow(tip_data))

# Match novel MAGs
cat("=== MATCHING NOVEL MAGs ===\\n")
novel_matches <- 0
for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  if (tip %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    next
  }}
  
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    next
  }}
  
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
  }}
}}

cat("Novel MAGs matched:", novel_matches, "out of", length(novel_mags), "\\n")

# FIXED: Define novel_count here
novel_count <- sum(tip_data$is_novel)
cat("Total novel MAGs count:", novel_count, "\\n")

# Match functional data for each type
for (func_type in names(functional_data)) {{
  cat(paste("=== MATCHING", toupper(func_type), "DATA ===\\n"))
  func_matches <- 0
  
  for (i in 1:nrow(tip_data)) {{
    tip <- tip_data$ID[i]
    col_name <- paste0(func_type, "_data")
    
    # Direct match
    if (tip %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip]]$categories
      func_matches <- func_matches + 1
      next
    }}
    
    # Conversion match
    tip_converted <- gsub("\\\\.", "_", tip)
    if (tip_converted %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip_converted]]$categories
      func_matches <- func_matches + 1
      next
    }}
    
    # Base match
    tip_base <- gsub("_sub$", "", tip)
    if (tip_base %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip_base]]$categories
      func_matches <- func_matches + 1
    }}
  }}
  
  cat(paste(toupper(func_type), "data matched:", func_matches, "out of", length(functional_data[[func_type]]), "\\n"))
}}

# Calculate functional data counts for legend
cazyme_count <- if("cazyme" %in% names(functional_data)) sum(sapply(tip_data[["cazyme_data"]], function(x) !is.null(x))) else 0
cog_count <- if("cog" %in% names(functional_data)) sum(sapply(tip_data[["cog_data"]], function(x) !is.null(x))) else 0
kegg_count <- if("kegg" %in% names(functional_data)) sum(sapply(tip_data[["kegg_data"]], function(x) !is.null(x))) else 0

cat(sprintf("Functional data summary: CAZyme=%d, COG=%d, KEGG=%d\\n", cazyme_count, cog_count, kegg_count))

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

# Define functional annotation color schemes
cazyme_colors <- c("GT" = "#E31A1C", "CE" = "#1F78B4", "GH" = "#33A02C", 
                  "CBM" = "#FF7F00", "PL" = "#6A3D9A", "AA" = "#FFD700")

cog_colors <- c("D" = "#E31A1C", "M" = "#1F78B4", "N" = "#33A02C", "O" = "#FF7F00",
               "T" = "#6A3D9A", "U" = "#FFD700", "V" = "#FB9A99", "W" = "#A6CEE3", 
               "Y" = "#B2DF8A", "Z" = "#FDBF6F")

kegg_colors <- c("Metabolism" = "#E31A1C", "Genetic Info" = "#1F78B4", "Environ Info" = "#33A02C",
                "Cellular" = "#FF7F00", "Organismal" = "#6A3D9A", "Diseases" = "#FFD700",
                "Brite" = "#FB9A99", "Other" = "#A6CEE3")

# Helper function to draw pie charts
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
          col = color, border = "white", lwd = 0.3)
}}

# Draw pie charts for functional annotations
draw_functional_pie <- function(center_x, center_y, radius, func_data, color_scheme) {{
  if (is.null(func_data) || length(func_data) == 0) return()
  
  values <- sapply(names(color_scheme), function(cat) {{
    val <- func_data[[cat]]
    if (is.null(val) || is.na(val)) 0 else as.numeric(val)
  }})
  
  if (sum(values) > 0) {{
    angles <- cumsum(values) / sum(values) * 2 * pi
    start_angles <- c(0, angles[-length(angles)])
    
    for (i in 1:length(names(color_scheme))) {{
      if (values[i] > 0) {{
        # PIE CHART (not donut) - inner radius is 0
        draw_arc_segment(center_x, center_y, 0, radius, 
                        start_angles[i], angles[i], color_scheme[names(color_scheme)[i]])
      }}
    }}
    
    # Add border around the pie
    symbols(center_x, center_y, circles = radius, add = TRUE, 
            inches = FALSE, fg = "black", lwd = 0.8)
  }}
}}

##############################################
# CREATE RECTANGULAR TREE WITH PIE CHARTS
##############################################

rectangular_file <- paste0(output_prefix, "_rectangular_pies.pdf")

# Calculate page height based on number of tips
page_height <- max(14, length(tree$tip.label) * 0.3)
pdf(rectangular_file, width = 24, height = page_height)

# Layout for tree and legend
layout(matrix(c(1,2), nrow=1), widths=c(4,1))

# Main tree plot with more space for pie charts
par(mar=c(3,1,3,16), xpd=TRUE)

# Plot rectangular tree without tip labels
plot(tree, type="phylogram", show.tip.label=FALSE, 
     main="Rectangular Phylogenetic Tree with Functional Annotation Pie Charts", 
     edge.width=2, cex.main=1.6, edge.color="#2C3E50",
     x.lim=c(0, max(node.depth.edgelength(tree)) * 1.6))

# Get tree layout information
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
n_tips <- length(tree$tip.label)

# Calculate positions for elements
x_tree_max <- max(last_plot$xx[1:n_tips])
strip_start_x <- x_tree_max + 0.02 * max(last_plot$xx)
strip_width <- 0.025 * max(last_plot$xx)
strip_gap <- 0.005 * max(last_plot$xx)

# Calculate pie chart positions (after taxonomic strips)
pie_radius <- min(0.4, 10 / n_tips)  # Scale pie size based on number of tips
pie_gap <- 0.015 * max(last_plot$xx)

# Position pie charts
cazyme_x <- strip_start_x + 2*strip_width + strip_gap + 3*pie_gap + pie_radius
cog_x <- cazyme_x + 2*pie_radius + pie_gap
kegg_x <- cog_x + 2*pie_radius + pie_gap

# Add tip labels first (before strips and pies)
for (i in 1:n_tips) {{
  tip <- tree$tip.label[i]
  y_pos <- last_plot$yy[i]
  
  # Draw tip label
  text(x_tree_max + 0.005 * max(last_plot$xx), y_pos, tip, 
       pos = 4, cex = 0.7, font = 1)
}}

# Draw elements for each tip
for (i in 1:n_tips) {{
  tip <- tree$tip.label[i]
  y_pos <- last_plot$yy[i]
  
  # Get data for this tip
  tip_idx <- which(tip_data$ID == tip)
  if (length(tip_idx) > 0) {{
    phylum <- tip_data$Phylum[tip_idx[1]]
    genus <- tip_data$Genus[tip_idx[1]]
    is_novel <- tip_data$is_novel[tip_idx[1]]
    
    # Draw taxonomic strips
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
    
    # Draw functional annotation pie charts
    # CAZyme pie
    if ("cazyme" %in% names(functional_data)) {{
      draw_functional_pie(cazyme_x, y_pos, pie_radius, 
                         tip_data$cazyme_data[[tip_idx[1]]], cazyme_colors)
    }}
    
    # COG pie
    if ("cog" %in% names(functional_data)) {{
      draw_functional_pie(cog_x, y_pos, pie_radius, 
                         tip_data$cog_data[[tip_idx[1]]], cog_colors)
    }}
    
    # KEGG pie
    if ("kegg" %in% names(functional_data)) {{
      draw_functional_pie(kegg_x, y_pos, pie_radius, 
                         tip_data$kegg_data[[tip_idx[1]]], kegg_colors)
    }}
    
    # Novel MAG indicator
    if (is_novel) {{
      # Add red border around taxonomic strips
      rect(strip_start_x - 0.001 * max(last_plot$xx), y_pos - 0.45, 
           strip_start_x + 2*strip_width + strip_gap + 0.001 * max(last_plot$xx), y_pos + 0.45,
           col = NA, border = "#FF0000", lwd = 2)
      
      # Add star symbol before pie charts
      points(strip_start_x + 2*strip_width + strip_gap + 0.015 * max(last_plot$xx), 
             y_pos, pch = 8, col = "#FF0000", cex = 1.5, lwd = 2)
    }}
  }}
}}

# Add column headers
header_y <- n_tips + 1.5
text(strip_start_x + strip_width/2, header_y, "Phylum", 
     font = 2, cex = 1.2, srt = 0)
text(strip_start_x + strip_width + strip_gap + strip_width/2, header_y, "Genus", 
     font = 2, cex = 1.2, srt = 0)

if ("cazyme" %in% names(functional_data)) {{
  text(cazyme_x, header_y, "CAZyme", font = 2, cex = 1.2, srt = 0)
}}
if ("cog" %in% names(functional_data)) {{
  text(cog_x, header_y, "COG", font = 2, cex = 1.2, srt = 0)
}}
if ("kegg" %in% names(functional_data)) {{
  text(kegg_x, header_y, "KEGG", font = 2, cex = 1.2, srt = 0)
}}

# COMPREHENSIVE LEGEND - FIXED AND COMPLETE
par(mar=c(5,2,5,2), xpd=TRUE)
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Title
text(0.5, 0.98, "Legend", font=2, cex=1.6, adj=0.5)

# SECTION 1: Element Descriptions
text(0.05, 0.92, "ELEMENT DESCRIPTIONS", font=2, pos=4, cex=1.2, col="darkblue")
legend_y <- 0.88

text(0.05, legend_y, "Column 1:", font=2, pos=4, cex=0.95)
text(0.25, legend_y, "Phylum", pos=4, cex=0.9)
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Column 2:", font=2, pos=4, cex=0.95)
text(0.25, legend_y, "Genus", pos=4, cex=0.9)
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Pie Charts:", font=2, pos=4, cex=0.95)
text(0.25, legend_y, "Functional Categories", pos=4, cex=0.9)
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Red Border:", font=2, pos=4, cex=0.95)
text(0.25, legend_y, "Novel MAGs", pos=4, cex=0.9, col="#FF0000")

# SECTION 2: Taxonomy Colors
legend_y <- legend_y - 0.05
text(0.05, legend_y, "TAXONOMY", font=2, pos=4, cex=1.2, col="darkblue")
legend_y <- legend_y - 0.04

# Phylum colors
text(0.05, legend_y, "Phylum:", font=2, pos=4, cex=1.0)
legend_y <- legend_y - 0.03
n_phyla_shown <- min(8, length(phyla))
for (i in 1:n_phyla_shown) {{
  rect(0.05, legend_y - 0.012, 0.10, legend_y + 0.012, 
       col=phylum_colors[phyla[i]], border="black", lwd=0.8)
  phylum_name <- phyla[i]
  if (nchar(phylum_name) > 20) {{
    phylum_name <- paste0(substr(phylum_name, 1, 18), "...")
  }}
  text(0.12, legend_y, phylum_name, pos=4, cex=0.8)
  legend_y <- legend_y - 0.028
}}
if (length(phyla) > n_phyla_shown) {{
  text(0.05, legend_y, paste0("... and ", length(phyla) - n_phyla_shown, " more"), 
       pos=4, cex=0.7, font=3)
  legend_y <- legend_y - 0.025
}}

# Genus colors (fewer shown)
legend_y <- legend_y - 0.03
text(0.05, legend_y, "Genus:", font=2, pos=4, cex=1.0)
legend_y <- legend_y - 0.03
n_genera_shown <- min(5, length(genera))
for (i in 1:n_genera_shown) {{
  rect(0.05, legend_y - 0.012, 0.10, legend_y + 0.012, 
       col=genus_colors[genera[i]], border="black", lwd=0.8)
  genus_name <- genera[i]
  if (nchar(genus_name) > 20) {{
    genus_name <- paste0(substr(genus_name, 1, 18), "...")
  }}
  text(0.12, legend_y, genus_name, pos=4, cex=0.8)
  legend_y <- legend_y - 0.028
}}
if (length(genera) > n_genera_shown) {{
  text(0.05, legend_y, paste0("... and ", length(genera) - n_genera_shown, " more"), 
       pos=4, cex=0.7, font=3)
  legend_y <- legend_y - 0.025
}}

# Novel MAG indicator
legend_y <- legend_y - 0.04
text(0.05, legend_y, "Novel MAGs:", font=2, pos=4, cex=1.0)
legend_y <- legend_y - 0.03
rect(0.05, legend_y - 0.012, 0.10, legend_y + 0.012, 
     col=NA, border="#FF0000", lwd=2)
points(0.075, legend_y, pch=8, col="#FF0000", cex=1.2, lwd=2)
text(0.12, legend_y, paste0("Novel (n=", novel_count, ")"), pos=4, cex=0.8)

# SECTION 3: Functional Annotations
legend_y <- legend_y - 0.05
text(0.05, legend_y, "FUNCTIONAL PIE CHARTS", font=2, pos=4, cex=1.2, col="darkblue")
legend_y <- legend_y - 0.04

# CAZyme categories
if ("cazyme" %in% names(functional_data)) {{
  text(0.05, legend_y, "CAZyme Categories:", font=2, pos=4, cex=1.0)
  legend_y <- legend_y - 0.03
  
  cazyme_names <- c("GT" = "Glycosyltransferases", 
                   "CE" = "Carbohydr. Esterases", 
                   "GH" = "Glycoside Hydrolases", 
                   "CBM" = "Carbohydr.-Binding",
                   "PL" = "Polysaccharide Lyases", 
                   "AA" = "Auxiliary Activities")
  
  for (i in 1:length(cazyme_colors)) {{
    cat_code <- names(cazyme_colors)[i]
    rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
         col=cazyme_colors[cat_code], border="black", lwd=0.6)
    text(0.09, legend_y, paste(cat_code, "-", cazyme_names[cat_code]), pos=4, cex=0.7)
    legend_y <- legend_y - 0.023
  }}
  legend_y <- legend_y - 0.02
}}

# COG categories
if ("cog" %in% names(functional_data)) {{
  text(0.05, legend_y, "COG Categories:", font=2, pos=4, cex=1.0)
  legend_y <- legend_y - 0.03
  
  cog_names <- c("D" = "Cell cycle control", 
                "M" = "Cell wall/membrane", 
                "N" = "Cell motility", 
                "O" = "Post-translational modif.",
                "T" = "Signal transduction", 
                "U" = "Intracellular trafficking", 
                "V" = "Defense mechanisms", 
                "W" = "Extracellular structures",
                "Y" = "Nuclear structure",
                "Z" = "Cytoskeleton")
  
  shown_cogs <- 0
  for (cat_code in names(cog_colors)) {{
    if (cat_code %in% names(cog_names) && shown_cogs < 6) {{
      rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
           col=cog_colors[cat_code], border="black", lwd=0.6)
      text(0.09, legend_y, paste(cat_code, "-", cog_names[cat_code]), pos=4, cex=0.7)
      legend_y <- legend_y - 0.023
      shown_cogs <- shown_cogs + 1
    }}
  }}
  legend_y <- legend_y - 0.02
}}

# KEGG pathways
if ("kegg" %in% names(functional_data)) {{
  text(0.05, legend_y, "KEGG Pathways:", font=2, pos=4, cex=1.0)
  legend_y <- legend_y - 0.03
  
  for (pathway in names(kegg_colors)) {{
    rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
         col=kegg_colors[pathway], border="black", lwd=0.6)
    text(0.09, legend_y, pathway, pos=4, cex=0.7)
    legend_y <- legend_y - 0.023
    if (legend_y < 0.15) break  # Stop if running out of space
  }}
}}

# SECTION 4: Statistics at bottom
legend_y <- 0.12
text(0.05, legend_y, "STATISTICS", font=2, pos=4, cex=1.2, col="darkblue")
legend_y <- legend_y - 0.03

text(0.05, legend_y, paste0("Total MAGs: ", n_tips), pos=4, cex=0.85)
legend_y <- legend_y - 0.022
text(0.05, legend_y, paste0("Novel MAGs: ", novel_count), pos=4, cex=0.85, col="#FF0000")
legend_y <- legend_y - 0.022
text(0.05, legend_y, paste0("Phyla: ", length(phyla)), pos=4, cex=0.85)
legend_y <- legend_y - 0.022
text(0.05, legend_y, paste0("Genera: ", length(genera)), pos=4, cex=0.85)

if ("cazyme" %in% names(functional_data)) {{
  legend_y <- legend_y - 0.022
  text(0.05, legend_y, paste0("With CAZyme data: ", cazyme_count), pos=4, cex=0.85)
}}
if ("cog" %in% names(functional_data)) {{
  legend_y <- legend_y - 0.022
  text(0.05, legend_y, paste0("With COG data: ", cog_count), pos=4, cex=0.85)
}}
if ("kegg" %in% names(functional_data)) {{
  legend_y <- legend_y - 0.022
  text(0.05, legend_y, paste0("With KEGG data: ", kegg_count), pos=4, cex=0.85)
}}

dev.off()

cat("\\n=== Rectangular Tree with Pie Charts Complete ===\\n")
cat("Created:", rectangular_file, "\\n")
cat(sprintf("Summary: %d total MAGs, %d novel MAGs\\n", n_tips, novel_count))
for (func_type in names(functional_data)) {{
  func_count <- sum(sapply(tip_data[[paste0(func_type, "_data")]], function(x) !is.null(x)))
  cat(sprintf("%s data: %d MAGs\\n", toupper(func_type), func_count))
}}

quit(status=0)
'''
        
        # Write and execute R script
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        os.chmod(script_path, 0o755)
        
        # Run the R script
        success = run_r_command_with_conda(script_path)
        
        if success:
            logger.info("Successfully created rectangular tree with pie charts")
        else:
            logger.error("Failed to create rectangular tree visualization")
        
        return success
        
    except Exception as e:
        logger.error(f"Error creating rectangular tree with pie charts: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_multi_functional_visualization(tree_file, metadata_file, output_prefix, 
                                        novel_mags, functional_data):
    """Create multi-functional visualization with CAZyme, COG, and KEGG pie charts."""
    try:
        script_path = f"{output_prefix}_multi_functional_vis.R"
        
        # Show what data we have
        logger.info(f"Creating multi-functional visualization with:")
        for func_type, data in functional_data.items():
            logger.info(f"  - {len(data)} {func_type.upper()} entries")
        logger.info(f"  - {len(novel_mags)} novel MAG entries")
        
        # Prepare data for R script
        novel_mags_json = json.dumps(novel_mags)
        functional_data_json = json.dumps(functional_data)
        
        r_script = f'''#!/usr/bin/env Rscript

# Multi-Functional Enhanced Phylogenetic Tree Visualization
# Support for CAZyme, COG, and KEGG annotations
# FIXED: Changed donuts to PIE CHARTS by setting inner radius to 0

# Load required packages
tryCatch({{
  library(ape)
  library(RColorBrewer)
  library(jsonlite)
  cat("Base packages loaded successfully\\n")
  
  use_ggtree <- TRUE
  tryCatch({{
    library(ggplot2)
    library(ggtree)
    cat("ggtree loaded successfully\\n")
  }}, error = function(e) {{
    cat("ggtree not available, using base R visualization\\n")
    use_ggtree <<- FALSE
  }})
}}, error = function(e) {{
  cat("Error loading packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- fromJSON('{novel_mags_json}')
functional_data <- fromJSON('{functional_data_json}')

set.seed(42)

cat("=== Multi-Functional Enhanced Tree Visualization (PIE CHARTS) ===\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")
for (func_type in names(functional_data)) {{
  cat(paste("Functional data for", toupper(func_type), ":", length(functional_data[[func_type]]), "MAGs\\n"))
}}

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

# Tip matching function
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

# Match tree tips with taxa data
tip_data <- do.call(rbind, lapply(tree$tip.label, function(tip) {{
  return(match_tip_to_taxa(tip, taxa_data))
}}))

# Add novel MAG and functional data
tip_data$is_novel <- FALSE
tip_data$cazyme_data <- rep(list(NULL), nrow(tip_data))
tip_data$cog_data <- rep(list(NULL), nrow(tip_data))
tip_data$kegg_data <- rep(list(NULL), nrow(tip_data))

# Match novel MAGs
cat("=== MATCHING NOVEL MAGs ===\\n")
novel_matches <- 0
for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  if (tip %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    next
  }}
  
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    next
  }}
  
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
  }}
}}

cat("Novel MAGs matched:", novel_matches, "out of", length(novel_mags), "\\n")

# FIXED: Define novel_count here
novel_count <- sum(tip_data$is_novel)

# Match functional data for each type
for (func_type in names(functional_data)) {{
  cat(paste("=== MATCHING", toupper(func_type), "DATA ===\\n"))
  func_matches <- 0
  
  for (i in 1:nrow(tip_data)) {{
    tip <- tip_data$ID[i]
    col_name <- paste0(func_type, "_data")
    
    # Direct match
    if (tip %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip]]$categories
      func_matches <- func_matches + 1
      next
    }}
    
    # Conversion match
    tip_converted <- gsub("\\\\.", "_", tip)
    if (tip_converted %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip_converted]]$categories
      func_matches <- func_matches + 1
      next
    }}
    
    # Base match
    tip_base <- gsub("_sub$", "", tip)
    if (tip_base %in% names(functional_data[[func_type]])) {{
      tip_data[[col_name]][[i]] <- functional_data[[func_type]][[tip_base]]$categories
      func_matches <- func_matches + 1
    }}
  }}
  
  cat(paste(toupper(func_type), "data matched:", func_matches, "out of", length(functional_data[[func_type]]), "\\n"))
}}

# Create color palettes
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

if (length(phyla) > 0) {{
  if (length(phyla) <= 9) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Set1")[1:length(phyla)]
  }} else {{
    phylum_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(phyla))
  }}
  names(phylum_colors) <- phyla
  phylum_colors <- c(phylum_colors, "Unknown" = "lightgray")
}} else {{
  phylum_colors <- c("Unknown" = "lightgray")
}}

if (length(genera) > 0) {{
  if (length(genera) <= 8) {{
    genus_colors <- brewer.pal(max(3, length(genera)), "Set2")[1:length(genera)]
  }} else {{
    genus_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "lightgray")
}} else {{
  genus_colors <- c("Unknown" = "lightgray")
}}

# Define functional annotation color schemes
cazyme_colors <- c("GT" = "#E31A1C", "CE" = "#1F78B4", "GH" = "#33A02C", 
                  "CBM" = "#FF7F00", "PL" = "#6A3D9A", "AA" = "#FFD700")

cog_colors <- c("D" = "#E31A1C", "M" = "#1F78B4", "N" = "#33A02C", "O" = "#FF7F00",
               "T" = "#6A3D9A", "U" = "#FFD700", "V" = "#FB9A99", "W" = "#A6CEE3", "Z" = "#B2DF8A")

kegg_colors <- c("Metabolism" = "#E31A1C", "Genetic Info" = "#1F78B4", "Environ Info" = "#33A02C",
                "Cellular" = "#FF7F00", "Organismal" = "#6A3D9A", "Diseases" = "#FFD700",
                "Brite" = "#FB9A99", "Other" = "#A6CEE3")

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
          col = color, border = "white", lwd = 0.3)
}}

# FIXED: Draw pie charts (NOT donuts) for functional annotations
draw_functional_pie <- function(center_x, center_y, radius, func_data, color_scheme) {{
  if (is.null(func_data) || length(func_data) == 0) return()
  
  values <- sapply(names(color_scheme), function(cat) {{
    val <- func_data[[cat]]
    if (is.null(val) || is.na(val)) 0 else as.numeric(val)
  }})
  
  if (sum(values) > 0) {{
    angles <- cumsum(values) / sum(values) * 2 * pi
    start_angles <- c(0, angles[-length(angles)])
    
    for (i in 1:length(names(color_scheme))) {{
      if (values[i] > 0) {{
        # FIXED: Changed radius * 0.6 to 0 to make PIE CHART instead of donut
        draw_arc_segment(center_x, center_y, 0, radius, 
                        start_angles[i], angles[i], color_scheme[names(color_scheme)[i]])
      }}
    }}
    
    # Add border around the pie
    symbols(center_x, center_y, circles = radius, add = TRUE, 
            inches = FALSE, fg = "black", lwd = 0.8)
  }}
}}

##############################################
# CREATE MULTI-FUNCTIONAL VISUALIZATION WITH PIE CHARTS
##############################################

enhanced_file <- paste0(output_prefix, "_multi_functional.pdf")
# Larger size to accommodate multiple functional rings
pdf(enhanced_file, width = 24, height = 16)

# Layout with more space for comprehensive legend
layout(matrix(c(1,2), nrow=1), widths=c(2.8,1))
par(mar=c(0.5,0.5,2,0.5), oma=c(0,0,0,0))

# Draw tree with wider limits for multiple pie chart rings
plot(tree, type="fan", show.tip.label=FALSE, 
     main="Multi-Functional Phylogenetic Tree (CAZyme, COG, KEGG) - PIE CHARTS", 
     edge.width=2.2, cex.main=1.8, edge.color="black", 
     no.margin=FALSE, use.edge.length=FALSE,
     x.lim=c(-2.2, 2.2), y.lim=c(-2.2, 2.2))

# Get coordinates
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords_x <- last_plot$xx[1:length(tree$tip.label)]
tip_coords_y <- last_plot$yy[1:length(tree$tip.label)]

center_x <- mean(range(tip_coords_x))
center_y <- mean(range(tip_coords_y))
max_radius <- max(sqrt((tip_coords_x - center_x)^2 + (tip_coords_y - center_y)^2))

# Calculate angles for each tip
tip_angles <- atan2(tip_coords_y - center_y, tip_coords_x - center_x)
tip_angles[tip_angles < 0] <- tip_angles[tip_angles < 0] + 2 * pi

# Sort tips by angle
sorted_indices <- order(tip_angles)
sorted_tips <- tree$tip.label[sorted_indices]
sorted_angles <- tip_angles[sorted_indices]
sorted_tips <- c(sorted_tips, sorted_tips[1])
sorted_angles <- c(sorted_angles, sorted_angles[1] + 2*pi)

# Define ring positions with space for 6 rings total
ring_gap <- 0.12 * max_radius
ring_width <- 0.08 * max_radius
ring_separation <- 0.04 * max_radius

# Ring positions
phylum_inner <- max_radius + ring_gap
phylum_outer <- phylum_inner + ring_width

genus_inner <- phylum_outer + ring_separation
genus_outer <- genus_inner + ring_width

novel_inner <- genus_outer + ring_separation
novel_outer <- novel_inner + ring_width * 0.6

# Functional annotation pie chart positions
cazyme_radius <- 0.05 * max_radius
cog_radius <- 0.05 * max_radius
kegg_radius <- 0.05 * max_radius

cazyme_distance <- novel_outer + ring_separation + cazyme_radius
cog_distance <- cazyme_distance + ring_separation + cazyme_radius + cog_radius
kegg_distance <- cog_distance + ring_separation + cog_radius + kegg_radius

# Draw basic rings (phylum, genus, novel)
for (i in 1:(length(sorted_tips)-1)) {{
  tryCatch({{
    tip <- sorted_tips[i]
    start_angle <- sorted_angles[i]
    end_angle <- sorted_angles[i+1]
    
    # Ring 1: Phylum color
    phylum <- get_taxon(tip, tip_to_phylum)
    phylum_color <- get_color(phylum, phylum_colors)
    draw_arc_segment(center_x, center_y, phylum_inner, phylum_outer, 
                    start_angle, end_angle, phylum_color)
    
    # Ring 2: Genus color
    genus <- get_taxon(tip, tip_to_genus)
    genus_color <- get_color(genus, genus_colors)
    draw_arc_segment(center_x, center_y, genus_inner, genus_outer, 
                    start_angle, end_angle, genus_color)
    
    # Ring 3: Novel MAG indicator
    is_novel <- get_taxon(tip, tip_to_novel)
    novel_color <- if (is_novel) "#FF0000" else "#F0F0F0"
    draw_arc_segment(center_x, center_y, novel_inner, novel_outer, 
                    start_angle, end_angle, novel_color)
    
  }}, error = function(e) {{
    cat("Error drawing basic rings for segment", i, ":", e$message, "\\n")
  }})
}}

# Draw functional annotation pie charts (NOT donuts)
for (i in 1:(length(sorted_tips)-1)) {{
  tip <- sorted_tips[i]
  angle <- (sorted_angles[i] + sorted_angles[i+1]) / 2
  
  tip_idx <- which(tip_data$ID == tip)
  if (length(tip_idx) > 0) {{
    tip_row <- tip_data[tip_idx[1], ]
    
    # CAZyme pie chart
    if ("cazyme" %in% names(functional_data)) {{
      cazyme_x <- center_x + cazyme_distance * cos(angle)
      cazyme_y <- center_y + cazyme_distance * sin(angle)
      draw_functional_pie(cazyme_x, cazyme_y, cazyme_radius, tip_row$cazyme_data[[1]], cazyme_colors)
    }}
    
    # COG pie chart
    if ("cog" %in% names(functional_data)) {{
      cog_x <- center_x + cog_distance * cos(angle)
      cog_y <- center_y + cog_distance * sin(angle)
      draw_functional_pie(cog_x, cog_y, cog_radius, tip_row$cog_data[[1]], cog_colors)
    }}
    
    # KEGG pie chart
    if ("kegg" %in% names(functional_data)) {{
      kegg_x <- center_x + kegg_distance * cos(angle)
      kegg_y <- center_y + kegg_distance * sin(angle)
      draw_functional_pie(kegg_x, kegg_y, kegg_radius, tip_row$kegg_data[[1]], kegg_colors)
    }}
  }}
}}

# Add tip points with enhanced styling
for (i in 1:length(tree$tip.label)) {{
  tip <- tree$tip.label[i]
  phylum <- get_taxon(tip, tip_to_phylum)
  phylum_color <- get_color(phylum, phylum_colors)
  
  points(tip_coords_x[i], tip_coords_y[i], pch=19, col=phylum_color, cex=2.0)
  points(tip_coords_x[i], tip_coords_y[i], pch=1, col="black", cex=2.0, lwd=1.5)
  
  if (get_taxon(tip, tip_to_novel)) {{
    points(tip_coords_x[i], tip_coords_y[i], pch=1, col="#FF0000", cex=5, lwd=3)
    points(tip_coords_x[i], tip_coords_y[i], pch=8, col="#FF0000", cex=4, lwd=2)
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
    
    lines(c(x1, x2), c(y1, y2), col="#D0D0D0", lwd=1.2, lty=1)
  }}
}}

# COMPREHENSIVE LEGEND - WITH FONT HIERARCHY: 1.8, 1.6, 1.5, 1.4
par(mar=c(2,1,3,1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Title - 1.8
text(0.5, 0.98, "Multi-Functional Analysis Legend", font=2, cex=1.8, adj=0.5)

# SECTION 1: Ring Descriptions - HEADING 1.6
text(0.05, 0.93, "RING DESCRIPTIONS", font=2, pos=4, cex=1.6, col="darkblue")
legend_y <- 0.89

# Ring descriptions - TEXT 1.4
text(0.05, legend_y, "Ring 1 (Inner):", font=2, pos=4, cex=1.4)
text(0.45, legend_y, "Phylum", pos=4, cex=1.4, col="darkgreen")
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Ring 2 (Middle):", font=2, pos=4, cex=1.4)
text(0.45, legend_y, "Genus", pos=4, cex=1.4, col="darkgreen")
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Ring 3 (Outer):", font=2, pos=4, cex=1.4)
text(0.45, legend_y, "Novel MAG Status", pos=4, cex=1.4, col="darkgreen")
legend_y <- legend_y - 0.025

text(0.05, legend_y, "Pie Charts:", font=2, pos=4, cex=1.4)
text(0.45, legend_y, "Functional Categories", pos=4, cex=1.4, col="darkgreen")

# SECTION 2: Taxonomy Colors - HEADING 1.6
legend_y <- legend_y - 0.04
text(0.05, legend_y, "TAXONOMY", font=2, pos=4, cex=1.6, col="darkblue")
legend_y <- legend_y - 0.025

# Genus - SUBHEADING 1.5
text(0.05, legend_y, "Genus:", font=2, pos=4, cex=1.5)
legend_y <- legend_y - 0.02
n_genera_to_show <- min(8, length(genera))  # Show only 4 to save space

genera_col1 <- ceiling(n_genera_to_show / 2)
genera_col2 <- n_genera_to_show - genera_col1

# Column 1 - TEXT 1.4
for (i in 1:genera_col1) {{
  rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
       col=genus_colors[genera[i]], border="black", lwd=0.8)
  genus_display <- if(nchar(genera[i]) > 10) paste0(substr(genera[i], 1, 8), "..") else genera[i]
  text(0.09, legend_y, genus_display, pos=4, cex=1.4)
  legend_y <- legend_y - 0.023
}}

# Column 2 - TEXT 1.4
if (genera_col2 > 0) {{
  legend_y_col2 <- legend_y + (genera_col1 * 0.023)
  for (i in (genera_col1+1):n_genera_to_show) {{
    rect(0.40, legend_y_col2 - 0.01, 0.43, legend_y_col2 + 0.01, 
         col=genus_colors[genera[i]], border="black", lwd=0.8)
    genus_display <- if(nchar(genera[i]) > 10) paste0(substr(genera[i], 1, 8), "..") else genera[i]
    text(0.44, legend_y_col2, genus_display, pos=4, cex=1.4)
    legend_y_col2 <- legend_y_col2 - 0.023
  }}
}}

if (length(genera) > n_genera_to_show) {{
  text(0.05, legend_y, paste0("+", length(genera) - n_genera_to_show, " more"), 
       pos=4, cex=1.4, font=3)
  legend_y <- legend_y - 0.02
}}

# Phylum - SUBHEADING 1.5
legend_y <- legend_y - 0.02
text(0.05, legend_y, "Phylum:", font=2, pos=4, cex=1.5)
legend_y <- legend_y - 0.02
n_phyla_to_show <- min(6, length(phyla))  # Show only 4 to save space

phyla_col1 <- ceiling(n_phyla_to_show / 2)
phyla_col2 <- n_phyla_to_show - phyla_col1

# Column 1 - TEXT 1.4
for (i in 1:phyla_col1) {{
  rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
       col=phylum_colors[phyla[i]], border="black", lwd=0.8)
  phylum_display <- if(nchar(phyla[i]) > 10) paste0(substr(phyla[i], 1, 8), "..") else phyla[i]
  text(0.09, legend_y, phylum_display, pos=4, cex=1.4)
  legend_y <- legend_y - 0.023
}}

# Column 2 - TEXT 1.4
if (phyla_col2 > 0) {{
  legend_y_col2 <- legend_y + (phyla_col1 * 0.023)
  for (i in (phyla_col1+1):n_phyla_to_show) {{
    rect(0.40, legend_y_col2 - 0.01, 0.43, legend_y_col2 + 0.01, 
         col=phylum_colors[phyla[i]], border="black", lwd=0.8)
    phylum_display <- if(nchar(phyla[i]) > 10) paste0(substr(phyla[i], 1, 8), "..") else phyla[i]
    text(0.44, legend_y_col2, phylum_display, pos=4, cex=1.4)
    legend_y_col2 <- legend_y_col2 - 0.023
  }}
}}

if (length(phyla) > n_phyla_to_show) {{
  text(0.05, legend_y, paste0("+", length(phyla) - n_phyla_to_show, " more"), 
       pos=4, cex=1.4, font=3)
  legend_y <- legend_y - 0.02
}}

# Novel MAG indicator - SUBHEADING 1.5
legend_y <- legend_y - 0.02
text(0.05, legend_y, "Novel MAGs:", font=2, pos=4, cex=1.5)
legend_y <- legend_y - 0.02
rect(0.05, legend_y - 0.01, 0.08, legend_y + 0.01, 
     col="#FF0000", border="black", lwd=1.0)
points(0.065, legend_y, pch=8, col="#FF0000", cex=1.5)
text(0.09, legend_y, paste0("Novel (n=", novel_count, ")"), pos=4, cex=1.4, col="#FF0000")

# SECTION 3: Functional Annotations - HEADING 1.6
legend_y <- legend_y - 0.05
text(0.05, legend_y, "FUNCTIONAL PIE CHARTS", font=2, pos=4, cex=1.6, col="darkblue")
legend_y <- legend_y - 0.025

# CAZyme categories - SUBHEADING 1.5
if ("cazyme" %in% names(functional_data)) {{
  text(0.05, legend_y, "CAZyme:", font=2, pos=4, cex=1.5)
  legend_y <- legend_y - 0.02
  
  cazyme_abbrev <- c("GT"="Glycosyltransf.", 
                     "CE"="Carb. Esterases", 
                     "GH"="Glyc. Hydrolases", 
                     "CBM"="Carb. Binding",
                     "PL"="Polysac. Lyases", 
                     "AA"="Aux. Activities")
  
  # Two columns - TEXT 1.4
  legend_y_start <- legend_y
  for (i in 1:3) {{
    cat_code <- names(cazyme_colors)[i]
    rect(0.05, legend_y - 0.008, 0.07, legend_y + 0.008, 
         col=cazyme_colors[cat_code], border="black", lwd=0.8)
    text(0.075, legend_y, paste0(cat_code, ":", cazyme_abbrev[cat_code]), pos=4, cex=1.4)
    legend_y <- legend_y - 0.022
  }}
  
  legend_y <- legend_y_start
  for (i in 4:6) {{
    if (i <= length(cazyme_colors)) {{
      cat_code <- names(cazyme_colors)[i]
      rect(0.55, legend_y - 0.008, 0.57, legend_y + 0.008, 
           col=cazyme_colors[cat_code], border="black", lwd=0.8)
      text(0.575, legend_y, paste0(cat_code, ":", cazyme_abbrev[cat_code]), pos=4, cex=1.4)
      legend_y <- legend_y - 0.022
    }}
  }}
  legend_y <- legend_y - 0.025
}}

# COG categories - SUBHEADING 1.5
if ("cog" %in% names(functional_data)) {{
  text(0.05, legend_y, "COG:", font=2, pos=4, cex=1.5)
  legend_y <- legend_y - 0.02
  
  cog_abbrev <- c("D"="Cell cycle", "M"="Cell wall", "N"="Motility", 
                 "O"="Protein mod.", "T"="Signal trans.", "U"="Secretion",
                 "V"="Defense", "W"="Extracell.", "Z"="Cytoskel.")
  
  n_cogs_to_show <- min(9, length(cog_colors))
  cogs_to_display <- names(cog_colors)[names(cog_colors) %in% names(cog_abbrev)][1:n_cogs_to_show]
  
  # Two columns - TEXT 1.4
  legend_y_start <- legend_y
  col1_items <- ceiling(length(cogs_to_display) / 2)
  
  for (i in 1:col1_items) {{
    if (i <= length(cogs_to_display)) {{
      cat_code <- cogs_to_display[i]
      rect(0.05, legend_y - 0.008, 0.07, legend_y + 0.008, 
           col=cog_colors[cat_code], border="black", lwd=0.8)
      text(0.075, legend_y, paste0(cat_code, ":", cog_abbrev[cat_code]), pos=4, cex=1.4)
      legend_y <- legend_y - 0.022
    }}
  }}
  
  legend_y <- legend_y_start
  for (i in (col1_items+1):length(cogs_to_display)) {{
    if (i <= length(cogs_to_display)) {{
      cat_code <- cogs_to_display[i]
      rect(0.40, legend_y - 0.008, 0.42, legend_y + 0.008, 
           col=cog_colors[cat_code], border="black", lwd=0.8)
      text(0.425, legend_y, paste0(cat_code, ":", cog_abbrev[cat_code]), pos=4, cex=1.4)
      legend_y <- legend_y - 0.022
    }}
  }}
  
  legend_y <- min(legend_y, legend_y_start - (col1_items * 0.022)) - 0.025
}}

# KEGG pathways - SUBHEADING 1.5
if ("kegg" %in% names(functional_data)) {{
  text(0.05, legend_y, "KEGG:", font=2, pos=4, cex=1.5)
  legend_y <- legend_y - 0.02
  
  kegg_full_names <- c("Metabolism"="Metabolic", 
                       "Genetic Info"="Genetic", 
                       "Environ Info"="Environ.",
                       "Cellular"="Cellular", 
                       "Organismal"="Organism.", 
                       "Diseases"="Diseases")
  
  n_kegg_to_show <- min(6, length(kegg_colors))
  kegg_pathways <- names(kegg_colors)[1:n_kegg_to_show]
  
  # Two columns - TEXT 1.4
  legend_y_start <- legend_y
  col1_items <- ceiling(n_kegg_to_show / 2)
  
  for (i in 1:col1_items) {{
    if (i <= length(kegg_pathways)) {{
      pathway <- kegg_pathways[i]
      rect(0.05, legend_y - 0.008, 0.07, legend_y + 0.008, 
           col=kegg_colors[pathway], border="black", lwd=0.8)
      display_name <- ifelse(pathway %in% names(kegg_full_names), 
                            kegg_full_names[pathway], pathway)
      text(0.075, legend_y, display_name, pos=4, cex=1.4)
      legend_y <- legend_y - 0.022
    }}
  }}
  
  legend_y <- legend_y_start
  for (i in (col1_items+1):n_kegg_to_show) {{
    if (i <= length(kegg_pathways)) {{
      pathway <- kegg_pathways[i]
      rect(0.40, legend_y - 0.008, 0.42, legend_y + 0.008, 
           col=kegg_colors[pathway], border="black", lwd=0.8)
      display_name <- ifelse(pathway %in% names(kegg_full_names), 
                            kegg_full_names[pathway], pathway)
      text(0.425, legend_y, display_name, pos=4, cex=1.4)
      legend_y <- legend_y - 0.022
    }}
  }}
  
  legend_y <- min(legend_y, legend_y_start - (col1_items * 0.022)) - 0.025
}}

# Move debug OUTSIDE the if statement
cat("DEBUG: 'kegg' in functional_data?", "kegg" %in% names(functional_data), "\\n")
cat("DEBUG: functional_data names:", paste(names(functional_data), collapse=", "), "\\n")

dev.off()

cat("\\n=== Multi-Functional Visualization Complete (PIE CHARTS) ===\\n")
cat("Created:", enhanced_file, "\\n")

# Final statistics
cat(sprintf("SUMMARY: %d total MAGs, %d novel MAGs highlighted\\n", nrow(tip_data), novel_count))
for (func_type in names(functional_data)) {{
  func_count <- sum(sapply(tip_data[[paste0(func_type, "_data")]], function(x) !is.null(x)))
  cat(sprintf("%s data: %d MAGs\\n", toupper(func_type), func_count))
}}

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
        logger.error(f"Error creating multi-functional visualization: {e}")
        return False

def create_publication_quality_visualization(tree_file, metadata_file, output_prefix, 
                                           novel_mags, functional_data, annotation_type):
    """Create publication-quality visualization with FIXED name matching and PIE CHARTS."""
    try:
        script_path = f"{output_prefix}_enhanced_vis.R"
        
        # Show matching debug information
        logger.info(f"Creating visualization with:")
        logger.info(f"  - {len(functional_data)} functional data entries")
        logger.info(f"  - {len(novel_mags)} novel MAG entries")
        
        # Show sample data for debugging
        if functional_data:
            sample_func = list(functional_data.keys())[:3]
            logger.info(f"  - Sample functional data keys: {sample_func}")
        
        if novel_mags:
            sample_novel = novel_mags[:3] 
            logger.info(f"  - Sample novel MAG keys: {sample_novel}")
        
        # Prepare data for R script
        novel_mags_json = json.dumps(novel_mags)
        functional_data_json = json.dumps(functional_data)
        
        r_script = f'''#!/usr/bin/env Rscript

# FIXED Publication-quality enhanced phylogenetic tree visualization
# Improvements: Fixed name matching, better layout, CAZy category PIE CHARTS (not donuts)

# Load required packages
tryCatch({{
  library(ape)
  library(RColorBrewer)
  library(jsonlite)
  cat("Base packages loaded successfully\\n")
  
  # Try to load ggtree
  use_ggtree <- TRUE
  tryCatch({{
    library(ggplot2)
    library(ggtree)
    cat("ggtree loaded successfully\\n")
  }}, error = function(e) {{
    cat("ggtree not available, using base R visualization\\n")
    use_ggtree <<- FALSE
  }})
}}, error = function(e) {{
  cat("Error loading packages:", e$message, "\\n")
  quit(status = 1)
}})

# Set file paths and data
tree_file <- "{tree_file}"
metadata_file <- "{metadata_file}"
output_prefix <- "{output_prefix}"
novel_mags <- fromJSON('{novel_mags_json}')
functional_data <- fromJSON('{functional_data_json}')
annotation_type <- "{annotation_type}"

set.seed(42)

cat("=== FIXED Enhanced Tree Visualization (PIE CHARTS) ===\\n")
cat("Novel MAGs to highlight:", length(novel_mags), "\\n")
cat("Functional data available for:", length(functional_data), "MAGs\\n")
cat("Annotation type:", annotation_type, "\\n")

# DEBUG: Show what data we have
if (length(novel_mags) > 0) {{
  cat("Sample novel MAGs:", paste(head(novel_mags, 3), collapse=", "), "\\n")
}}
if (length(functional_data) > 0) {{
  cat("Sample functional data keys:", paste(head(names(functional_data), 3), collapse=", "), "\\n")
}}

# Read and process tree
tryCatch({{
  tree <- read.tree(tree_file)
  cat("Tree loaded successfully with", length(tree$tip.label), "tips\\n")
  
  # DEBUG: Show tree tip names
  cat("Tree tip names:", paste(head(tree$tip.label, 5), collapse=", "), "\\n")
  
  # Apply rooting if needed
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
  }} else {{
    cat("Tree is already rooted\\n")
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
  cat("Using ID column:", id_col, "\\n")
}}, error = function(e) {{
  cat("Error loading metadata:", e$message, "\\n")
  quit(status=1)
}})

# Process taxonomy data
taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_data <- data.frame(ID = metadata[[id_col]])

# Create taxonomy columns
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

# FIXED: Better tip matching function
match_tip_to_taxa <- function(tip, taxa_data) {{
  # Try exact match first
  match_idx <- which(taxa_data$ID == tip)
  if (length(match_idx) == 0) {{
    # Try case-insensitive match
    match_idx <- which(tolower(taxa_data$ID) == tolower(tip))
  }}
  if (length(match_idx) == 0) {{
    # Try with underscore/dot conversion
    converted_tip <- gsub("\\\\.", "_", tip)  # Convert dots to underscores
    match_idx <- which(taxa_data$ID == converted_tip)
  }}
  if (length(match_idx) == 0) {{
    # Try partial matching
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
    # Return Unknown values if no match
    unknown_row <- taxa_data[1, ]
    unknown_row[] <- c(tip, rep("Unknown", length(taxonomic_levels)))
    return(unknown_row)
  }}
}}

# Match tree tips with taxa data
tip_data <- do.call(rbind, lapply(tree$tip.label, function(tip) {{
  return(match_tip_to_taxa(tip, taxa_data))
}}))

# CRITICAL FIX: Add novel MAG and CAZy category data with proper matching
tip_data$is_novel <- FALSE
tip_data$cazyme_categories <- rep(list(NULL), nrow(tip_data))

# Match novel MAGs - FIXED to handle name variations
cat("=== MATCHING NOVEL MAGs ===\\n")
novel_matches <- 0
for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  # Direct match
  if (tip %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    cat("Novel MAG direct match:", tip, "\\n")
    next
  }}
  
  # Try with underscore conversion
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    cat("Novel MAG converted match:", tip, "<->", tip_converted, "\\n")
    next
  }}
  
  # Try without _sub suffix
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% novel_mags) {{
    tip_data$is_novel[i] <- TRUE
    novel_matches <- novel_matches + 1
    cat("Novel MAG base match:", tip, "<->", tip_base, "\\n")
    next
  }}
}}

cat("Novel MAGs matched:", novel_matches, "out of", length(novel_mags), "\\n")

# CRITICAL FIX: Extract CAZy category data instead of individual families
cat("=== MATCHING FUNCTIONAL DATA ===\\n")
functional_matches <- 0

for (i in 1:nrow(tip_data)) {{
  tip <- tip_data$ID[i]
  
  # Direct match
  if (tip %in% names(functional_data)) {{
    func_data <- functional_data[[tip]]
    if (!is.null(func_data) && !is.null(func_data$cazy_families)) {{
      # Extract category totals (GT, CE, GH, CBM, PL, AA)
      categories <- c("GT", "CE", "GH", "CBM", "PL", "AA")
      category_counts <- list()
      
      for (cat in categories) {{
        total_key <- paste0(cat, "_total")
        if (total_key %in% names(func_data$cazy_families)) {{
          category_counts[[cat]] <- func_data$cazy_families[[total_key]]
        }} else {{
          category_counts[[cat]] <- 0
        }}
      }}
      
      tip_data$cazyme_categories[[i]] <- category_counts
      functional_matches <- functional_matches + 1
      cat("Functional direct match:", tip, "-> CAZy categories\\n")
    }}
    next
  }}
  
  # Try with underscore conversion
  tip_converted <- gsub("\\\\.", "_", tip)
  if (tip_converted %in% names(functional_data)) {{
    func_data <- functional_data[[tip_converted]]
    if (!is.null(func_data) && !is.null(func_data$cazy_families)) {{
      categories <- c("GT", "CE", "GH", "CBM", "PL", "AA")
      category_counts <- list()
      
      for (cat in categories) {{
        total_key <- paste0(cat, "_total")
        if (total_key %in% names(func_data$cazy_families)) {{
          category_counts[[cat]] <- func_data$cazy_families[[total_key]]
        }} else {{
          category_counts[[cat]] <- 0
        }}
      }}
      
      tip_data$cazyme_categories[[i]] <- category_counts
      functional_matches <- functional_matches + 1
      cat("Functional converted match:", tip, "<->", tip_converted, "-> CAZy categories\\n")
    }}
    next
  }}
  
  # Try without _sub suffix
  tip_base <- gsub("_sub$", "", tip)
  if (tip_base %in% names(functional_data)) {{
    func_data <- functional_data[[tip_base]]
    if (!is.null(func_data) && !is.null(func_data$cazy_families)) {{
      categories <- c("GT", "CE", "GH", "CBM", "PL", "AA")
      category_counts <- list()
      
      for (cat in categories) {{
        total_key <- paste0(cat, "_total")
        if (total_key %in% names(func_data$cazy_families)) {{
          category_counts[[cat]] <- func_data$cazy_families[[total_key]]
        }} else {{
          category_counts[[cat]] <- 0
        }}
      }}
      
      tip_data$cazyme_categories[[i]] <- category_counts
      functional_matches <- functional_matches + 1
      cat("Functional base match:", tip, "<->", tip_base, "-> CAZy categories\\n")
    }}
    next
  }}
}}

cat("Functional data matched:", functional_matches, "out of", length(functional_data), "\\n")

# Count and report results
novel_count <- sum(tip_data$is_novel)
functional_count <- sum(sapply(tip_data$cazyme_categories, function(x) !is.null(x)))
cat("FINAL RESULTS:\\n")
cat("Novel MAGs found:", novel_count, "\\n")
cat("MAGs with functional data:", functional_count, "\\n")

# Get unique taxa for color palettes
phyla <- unique(tip_data$Phylum[tip_data$Phylum != "Unknown"])
phyla <- phyla[!is.na(phyla)]
genera <- unique(tip_data$Genus[tip_data$Genus != "Unknown"])
genera <- genera[!is.na(genera)]

# Create robust color palettes
if (length(phyla) > 0) {{
  if (length(phyla) <= 9) {{
    phylum_colors <- brewer.pal(max(3, length(phyla)), "Set1")[1:length(phyla)]
  }} else {{
    phylum_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(phyla))
  }}
  names(phylum_colors) <- phyla
  phylum_colors <- c(phylum_colors, "Unknown" = "lightgray")
}} else {{
  phylum_colors <- c("Unknown" = "lightgray")
}}

if (length(genera) > 0) {{
  if (length(genera) <= 8) {{
    genus_colors <- brewer.pal(max(3, length(genera)), "Set2")[1:length(genera)]
  }} else {{
    genus_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(genera))
  }}
  names(genus_colors) <- genera
  genus_colors <- c(genus_colors, "Unknown" = "lightgray")
}} else {{
  genus_colors <- c("Unknown" = "lightgray")
}}

# Create maps from tip labels to taxonomy and functional data
tip_to_phylum <- setNames(tip_data$Phylum, tip_data$ID)
tip_to_genus <- setNames(tip_data$Genus, tip_data$ID)
tip_to_novel <- setNames(tip_data$is_novel, tip_data$ID)

cat("Created color palettes for", length(phyla), "phyla and", length(genera), "genera\\n")

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
          col = color, border = "white", lwd = 0.3)
}}

# FIXED: Draw PIE CHARTS (not donuts) for CAZyme categories
draw_cazyme_pie <- function(center_x, center_y, radius, cazyme_data) {{
  if (is.null(cazyme_data) || length(cazyme_data) == 0) return()
  
  categories <- c("GT", "CE", "GH", "CBM", "PL", "AA")
  category_colors <- c("GT" = "#E31A1C", "CE" = "#1F78B4", "GH" = "#33A02C", 
                      "CBM" = "#FF7F00", "PL" = "#6A3D9A", "AA" = "#FFD700")
  
  values <- sapply(categories, function(cat) {{
    val <- cazyme_data[[cat]]
    if (is.null(val) || is.na(val)) 0 else as.numeric(val)
  }})
  
  if (sum(values) > 0) {{
    angles <- cumsum(values) / sum(values) * 2 * pi
    start_angles <- c(0, angles[-length(angles)])
    
    for (i in 1:length(categories)) {{
      if (values[i] > 0) {{
        # FIXED: Changed radius * 0.6 to 0 to make PIE CHART instead of donut
        draw_arc_segment(center_x, center_y, 0, radius, 
                        start_angles[i], angles[i], category_colors[categories[i]])
      }}
    }}
    
    # Add border around the pie
    symbols(center_x, center_y, circles = radius, add = TRUE, 
            inches = FALSE, fg = "black", lwd = 1)
  }}
}}

##############################################
# CREATE IMPROVED PUBLICATION-QUALITY VISUALIZATION WITH PIE CHARTS
##############################################

enhanced_file <- paste0(output_prefix, "_enhanced_beautiful.pdf")
# FIXED: Better page fitting - larger size to accommodate all elements
pdf(enhanced_file, width = 20, height = 14)  # Wider format to fit everything

# Improved layout - more space for legend
layout(matrix(c(1,2), nrow=1), widths=c(2.5,1))  # Adjust proportions
par(mar=c(0.5,0.5,2,0.5), oma=c(0,0,0,0))

# Draw tree with adjusted limits to ensure everything fits
plot(tree, type="fan", show.tip.label=FALSE, 
     main="Enhanced Phylogenetic Tree with Novel MAGs and CAZyme Categories - PIE CHARTS", 
     edge.width=2.2, cex.main=1.6, edge.color="black", 
     no.margin=FALSE, use.edge.length=FALSE,
     x.lim=c(-1.8, 1.8), y.lim=c(-1.8, 1.8))  # Wider limits to fit pie charts

# Get coordinates
last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords_x <- last_plot$xx[1:length(tree$tip.label)]
tip_coords_y <- last_plot$yy[1:length(tree$tip.label)]

center_x <- mean(range(tip_coords_x))
center_y <- mean(range(tip_coords_y))
max_radius <- max(sqrt((tip_coords_x - center_x)^2 + (tip_coords_y - center_y)^2))

# Calculate angles for each tip
tip_angles <- atan2(tip_coords_y - center_y, tip_coords_x - center_x)
tip_angles[tip_angles < 0] <- tip_angles[tip_angles < 0] + 2 * pi

# Sort tips by angle
sorted_indices <- order(tip_angles)
sorted_tips <- tree$tip.label[sorted_indices]
sorted_angles <- tip_angles[sorted_indices]
sorted_tips <- c(sorted_tips, sorted_tips[1])
sorted_angles <- c(sorted_angles, sorted_angles[1] + 2*pi)

# Define improved ring positions with better spacing
ring_gap <- 0.15 * max_radius        # Gap from tree
ring_width <- 0.10 * max_radius      # Ring width
ring_separation <- 0.05 * max_radius # Separation between rings

# Ring positions
phylum_inner <- max_radius + ring_gap
phylum_outer <- phylum_inner + ring_width

genus_inner <- phylum_outer + ring_separation
genus_outer <- genus_inner + ring_width

novel_inner <- genus_outer + ring_separation
novel_outer <- novel_inner + ring_width * 0.8

# Draw segments for each tip (only rings, not pie charts here)
for (i in 1:(length(sorted_tips)-1)) {{
  tryCatch({{
    tip <- sorted_tips[i]
    start_angle <- sorted_angles[i]
    end_angle <- sorted_angles[i+1]
    
    # Ring 1: Phylum color
    phylum <- get_taxon(tip, tip_to_phylum)
    phylum_color <- get_color(phylum, phylum_colors)
    draw_arc_segment(center_x, center_y, phylum_inner, phylum_outer, 
                    start_angle, end_angle, phylum_color)
    
    # Ring 2: Genus color
    genus <- get_taxon(tip, tip_to_genus)
    genus_color <- get_color(genus, genus_colors)
    draw_arc_segment(center_x, center_y, genus_inner, genus_outer, 
                    start_angle, end_angle, genus_color)
    
    # Ring 3: Novel MAG indicator
    is_novel <- get_taxon(tip, tip_to_novel)
    novel_color <- if (is_novel) "#FF0000" else "#F0F0F0"
    draw_arc_segment(center_x, center_y, novel_inner, novel_outer, 
                    start_angle, end_angle, novel_color)
    
  }}, error = function(e) {{
    cat("Error drawing segment", i, ":", e$message, "\\n")
  }})
}}

# FIXED: Draw PIE CHARTS (not donuts) for each tip instead of functional ring
pie_radius <- 0.06 * max_radius
for (i in 1:(length(sorted_tips)-1)) {{
  tip <- sorted_tips[i]
  angle <- (sorted_angles[i] + sorted_angles[i+1]) / 2
  
  pie_x <- center_x + (novel_outer + ring_separation + pie_radius) * cos(angle)
  pie_y <- center_y + (novel_outer + ring_separation + pie_radius) * sin(angle)
  
  tip_idx <- which(tip_data$ID == tip)
  if (length(tip_idx) > 0) {{
    draw_cazyme_pie(pie_x, pie_y, pie_radius, tip_data$cazyme_categories[[tip_idx[1]]])
  }}
}}

# Add tip points with enhanced styling
for (i in 1:length(tree$tip.label)) {{
  tip <- tree$tip.label[i]
  phylum <- get_taxon(tip, tip_to_phylum)
  phylum_color <- get_color(phylum, phylum_colors)
  
  # Regular tip point
  points(tip_coords_x[i], tip_coords_y[i], pch=19, col=phylum_color, cex=2.0)
  points(tip_coords_x[i], tip_coords_y[i], pch=1, col="black", cex=2.0, lwd=1.5)
  
  # Novel MAG highlighting
  if (get_taxon(tip, tip_to_novel)) {{
    # Multiple layers for emphasis
    points(tip_coords_x[i], tip_coords_y[i], pch=1, col="#FF0000", cex=5, lwd=3)
    points(tip_coords_x[i], tip_coords_y[i], pch=8, col="#FF0000", cex=4, lwd=2)
  }}
}}

# Add connecting lines with better visibility
for (i in 1:length(sorted_tips) - 1) {{
  tip_idx <- which(tree$tip.label == sorted_tips[i])
  if (length(tip_idx) > 0) {{
    angle <- sorted_angles[i]
    x1 <- tip_coords_x[tip_idx]
    y1 <- tip_coords_y[tip_idx]
    x2 <- center_x + phylum_inner * cos(angle)
    y2 <- center_y + phylum_inner * sin(angle)
    
    lines(c(x1, x2), c(y1, y2), col="#D0D0D0", lwd=1.2, lty=1)
  }}
}}

# Enhanced legends with better formatting
par(mar=c(2,1,3,1))
plot(0, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

# Phylum legend
text(0.05, 0.95, "Phylum", font=2, pos=4, cex=1.6)
legend_y <- 0.88
for (i in 1:min(10, length(phyla))) {{
  rect(0.05, legend_y - 0.02, 0.15, legend_y + 0.02, 
       col=phylum_colors[phyla[i]], border="black", lwd=1)
  text(0.18, legend_y, phyla[i], pos=4, cex=1.1)
  legend_y <- legend_y - 0.04
}}

# FIXED: Add missing genus legend
legend_y <- legend_y - 0.06
text(0.05, legend_y, "Genus", font=2, pos=4, cex=1.6)
legend_y <- legend_y - 0.05
for (i in 1:min(6, length(genera))) {{
  rect(0.05, legend_y - 0.02, 0.15, legend_y + 0.02, 
       col=genus_colors[genera[i]], border="black", lwd=1)
  text(0.18, legend_y, genera[i], pos=4, cex=1.1)
  legend_y <- legend_y - 0.04
}}

# Novel MAG legend
legend_y <- legend_y - 0.06
text(0.05, legend_y, "Novel MAGs", font=2, pos=4, cex=1.6)
legend_y <- legend_y - 0.05
rect(0.05, legend_y - 0.02, 0.15, legend_y + 0.02, 
     col="#FF0000", border="black", lwd=1)
text(0.18, legend_y, paste("Novel (", novel_count, ")", sep=""), pos=4, cex=1.1)
legend_y <- legend_y - 0.04
rect(0.05, legend_y - 0.02, 0.15, legend_y + 0.02, 
     col="#F0F0F0", border="black", lwd=1)
text(0.18, legend_y, paste("Regular (", nrow(tip_data) - novel_count, ")", sep=""), pos=4, cex=1.1)

# FIXED: CAZyme categories legend
legend_y <- legend_y - 0.06
text(0.05, legend_y, "CAZyme Categories (PIE CHARTS)", font=2, pos=4, cex=1.6)
legend_y <- legend_y - 0.05

categories <- c("GT", "CE", "GH", "CBM", "PL", "AA")
category_colors <- c("GT" = "#E31A1C", "CE" = "#1F78B4", "GH" = "#33A02C", 
                    "CBM" = "#FF7F00", "PL" = "#6A3D9A", "AA" = "#FFD700")
category_names <- c("GT" = "Glycosyltransferases", "CE" = "Carbohydrate Esterases", 
                   "GH" = "Glycoside Hydrolases", "CBM" = "Carbohydrate-Binding Modules",
                   "PL" = "Polysaccharide Lyases", "AA" = "Auxiliary Activities")

for (i in 1:length(categories)) {{
  rect(0.05, legend_y - 0.02, 0.15, legend_y + 0.02, 
       col=category_colors[categories[i]], border="black", lwd=1)
  text(0.18, legend_y, paste(categories[i], "-", category_names[categories[i]]), pos=4, cex=0.9)
  legend_y <- legend_y - 0.035
}}

dev.off()

cat("\\n=== FIXED Enhanced Visualization Complete (PIE CHARTS) ===\\n")
cat("Created:", enhanced_file, "\\n")

# Final statistics
cat(sprintf("FINAL SUMMARY: %d total MAGs, %d novel MAGs highlighted, %d with functional data\\n",
           nrow(tip_data), novel_count, functional_count))

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
        logger.error(f"Error creating publication quality visualization: {e}")
        return False

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

# Import helper functions from original module
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
                import glob
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

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 6:
        print("Usage: python enhanced_phylo_tree_vis.py <cpus> <memory> <time_limit> <project_input> <project_output> [taxonomy_source] [novel_mags_file] [functional_annotations] [annotation_type] [cazyme_annotations] [cog_annotations] [kegg_annotations]")
        sys.exit(1)
    
    cpus = int(sys.argv[1])
    memory = sys.argv[2]
    time_limit = sys.argv[3]
    project_input = sys.argv[4]
    project_output = sys.argv[5]
    taxonomy_source = sys.argv[6] if len(sys.argv) > 6 else None
    novel_mags_file = sys.argv[7] if len(sys.argv) > 7 else None
    functional_annotations = sys.argv[8] if len(sys.argv) > 8 else None
    annotation_type = sys.argv[9] if len(sys.argv) > 9 else "auto"
    cazyme_annotations = sys.argv[10] if len(sys.argv) > 10 else None
    cog_annotations = sys.argv[11] if len(sys.argv) > 11 else None
    kegg_annotations = sys.argv[12] if len(sys.argv) > 12 else None
    
    success = run_enhanced_tree_visualization(
        cpus, memory, time_limit, project_input, project_output, 
        taxonomy_source, novel_mags_file, functional_annotations, annotation_type,
        cazyme_annotations, cog_annotations, kegg_annotations
    )
    sys.exit(0 if success else 1)
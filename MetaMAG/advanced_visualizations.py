#!/usr/bin/env python3
"""
Publication-Quality Advanced Visualizations for Metagenomic Pipeline
=====================================================================
Creates sophisticated, publication-ready visualizations with proper
taxonomy-function integration.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.gridspec import GridSpec
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
import logging
from pathlib import Path
import json
from collections import defaultdict, Counter
import re
from itertools import combinations, product
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
from sklearn.metrics import pairwise_distances
import colorsys

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not available, some visualizations will be skipped")

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    print("Warning: networkx not available, network plots will be skipped")

# Configure for publication quality
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.transparent': False
})

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='[%(asctime)s] [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# Professional color schemes
PHYLUM_COLORS = {
    'Pseudomonadota': '#E64B35',
    'Bacteroidota': '#4DBBD5',
    'Acidobacteriota': '#00A087',
    'Chloroflexota': '#3C5488',
    'Armatimonadota': '#F39B7F',
    'Actinomycetota': '#8491B4',
    'Bacillota': '#91D1C2',
    'Verrucomicrobiota': '#DC0000',
    'Desulfobacterota': '#7E6148',
    'Unknown': '#B09C85'
}

COG_CATEGORY_COLORS = {
    'INFORMATION STORAGE AND PROCESSING': '#3B4992',
    'CELLULAR PROCESSES AND SIGNALING': '#EE0000',
    'METABOLISM': '#008B45',
    'POORLY CHARACTERIZED': '#631879'
}

# COG categories mapping
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking',
    'O': 'Posttranslational modification',
    'X': 'Mobilome',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis',
    'R': 'General function prediction only',
    'S': 'Function unknown'
}

COG_MAIN_CATEGORIES = {
    "INFORMATION STORAGE AND PROCESSING": ["J", "A", "K", "L", "B"],
    "CELLULAR PROCESSES AND SIGNALING": ["D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X"],
    "METABOLISM": ["C", "G", "E", "F", "H", "I", "P", "Q"],
    "POORLY CHARACTERIZED": ["R", "S"]
}


class PublicationQualityVisualizations:
    """
    Generate publication-quality visualizations for metagenomic data
    """
    
    def __init__(self, project_dir, output_suffix=""):
        self.project_dir = Path(project_dir)
        self.output_suffix = output_suffix
        
        # Create output directory
        output_name = f"Advanced_visualizations{output_suffix}" if output_suffix else "Advanced_visualizations"
        self.output_dir = self.project_dir / output_name
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize data containers
        self.data = {}
        self.taxonomy = None
        self.mag_metadata = {}
        self.name_mapping = {}
        
        logger.info(f"Initialized publication-quality visualization module")
        logger.info(f"Output directory: {self.output_dir}")
    
    def convert_mag_names(self, func_name):
        """Convert between functional and taxonomy naming conventions"""
        name = func_name
        for suffix in ['.counts', '_fa.counts', '_fa_faa', '_faa', '_fa', '.faa', '.fa']:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        
        patterns = [
            (r'(.*_metabat)_bin_(\d+)', r'\1_bin.\2'),
            (r'(.*_maxbin2)_bin_(\d+)', r'\1_bin.\2'),
            (r'(.*_concoct)_bin_(\d+)', r'\1_bin.\2'),
            (r'(coassembly)_bin_(\d+)', r'\1_bin.\2'),
        ]
        
        for pattern, replacement in patterns:
            if re.search(pattern, name):
                name = re.sub(pattern, replacement, name)
                break
        
        return name
    
    def load_all_data(self):
        """Load and integrate all data sources"""
        logger.info("Loading all data sources...")
        
        self._load_functional_data()
        self._load_taxonomy_data()
        self._create_mag_metadata()
        self._create_name_mappings()
        
        logger.info("Data loading complete")
        return self.data
    
    def _load_functional_data(self):
        """Load functional annotation data"""
        cog_file = self.project_dir / "functional_analysis" / "files" / "cog" / "individual_cog_counts.xlsx"
        if cog_file.exists():
            self.data['cog'] = pd.read_excel(cog_file, index_col=0)
            logger.info(f"Loaded COG data: {self.data['cog'].shape}")
        
        ko_file = self.project_dir / "functional_analysis" / "files" / "ko" / "ko_counts.xlsx"
        if ko_file.exists():
            self.data['ko'] = pd.read_excel(ko_file, index_col=0)
            logger.info(f"Loaded KO data: {self.data['ko'].shape}")
        
        kegg_file = self.project_dir / "functional_analysis" / "files" / "kegg" / "ko_abundance_with_levels.xlsx"
        if kegg_file.exists():
            self.data['kegg_pathways'] = pd.read_excel(kegg_file)
            logger.info(f"Loaded KEGG pathway data: {self.data['kegg_pathways'].shape}")
        
        cazyme_file = self.project_dir / "cazyme_annotations" / "category_counts.xlsx"
        if cazyme_file.exists():
            self.data['cazyme'] = pd.read_excel(cazyme_file, index_col=0)
            logger.info(f"Loaded CAZyme data: {self.data['cazyme'].shape}")
    
    def _load_taxonomy_data(self):
        """Load and parse taxonomy data"""
        taxonomy_files = [
            self.project_dir / "Novel_Mags" / "gtdbtk" / "gtdbtk.bac120.summary.tsv",
            self.project_dir / "gtdbtk" / "gtdbtk.bac120.summary.tsv"
        ]
        
        for tax_file in taxonomy_files:
            if tax_file.exists():
                self.taxonomy = pd.read_csv(tax_file, sep='\t')
                self._parse_taxonomy()
                logger.info(f"Loaded taxonomy from {tax_file}")
                break
    
    def _parse_taxonomy(self):
        """Parse GTDB taxonomy strings"""
        if self.taxonomy is None:
            return
        
        tax_levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        
        for idx, row in self.taxonomy.iterrows():
            classification = row.get('classification', '')
            if pd.notna(classification):
                tax_dict = {}
                for part in classification.split(';'):
                    if '__' in part:
                        prefix, name = part.split('__', 1)
                        level_map = {'d': 'domain', 'p': 'phylum', 'c': 'class',
                                   'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
                        if prefix in level_map:
                            tax_dict[level_map[prefix]] = name if name else 'Unknown'
                
                for level in tax_levels:
                    self.taxonomy.loc[idx, level] = tax_dict.get(level, 'Unknown')
    
    def _create_mag_metadata(self):
        """Create comprehensive MAG metadata"""
        if 'cog' in self.data:
            mags = list(self.data['cog'].columns)
            for mag in mags:
                self.mag_metadata[mag] = {
                    'functional_name': mag,
                    'taxonomy_name': self.convert_mag_names(mag)
                }
        
        if self.taxonomy is not None:
            for _, row in self.taxonomy.iterrows():
                tax_name = row['user_genome']
                for mag, metadata in self.mag_metadata.items():
                    if metadata['taxonomy_name'] == tax_name:
                        for level in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                            if level in row:
                                metadata[level] = row[level]
                        break
    
    def _create_name_mappings(self):
        """Create bidirectional name mappings"""
        for mag, metadata in self.mag_metadata.items():
            self.name_mapping[metadata['functional_name']] = metadata['taxonomy_name']
    
    def generate_gradient_colors(self, base_color, n_colors):
        """Generate a gradient of colors based on a base color"""
        if n_colors == 0:
            return []
        # Convert hex to RGB
        base_rgb = tuple(int(base_color[i:i+2], 16)/255 for i in (1, 3, 5))
        h, s, v = colorsys.rgb_to_hsv(*base_rgb)
        
        colors = []
        for i in range(n_colors):
            # Vary saturation and value to create gradient
            new_s = s * (0.5 + 0.5 * (i / max(n_colors-1, 1)))
            new_v = v * (0.7 + 0.3 * (i / max(n_colors-1, 1)))
            new_rgb = colorsys.hsv_to_rgb(h, new_s, new_v)
            new_hex = '#{:02x}{:02x}{:02x}'.format(
                int(new_rgb[0]*255), int(new_rgb[1]*255), int(new_rgb[2]*255)
            )
            colors.append(new_hex)
        return colors

    def create_taxonomy_cazyme_heatmap(self):
        """
        Create sophisticated heatmap showing CAZyme profiles grouped by taxonomy
        Creates two versions: with and without dendrogram
        """
        logger.info("Creating enhanced taxonomy-CAZyme heatmaps...")
        
        if 'cazyme' not in self.data or self.taxonomy is None:
            logger.warning("Insufficient data for CAZyme heatmap")
            return
        
        # The data is loaded with MAGs as index and CAZyme categories as columns
        # We need to TRANSPOSE it to get MAGs as columns
        cazyme_data_raw = self.data['cazyme']
        
        # Check the structure
        logger.info(f"CAZyme data shape (before transpose): {cazyme_data_raw.shape}")
        
        # TRANSPOSE the data so MAGs are columns and CAZyme categories are rows
        cazyme_data = cazyme_data_raw.T
        
        # Remove 'Total' row if it exists
        if 'Total' in cazyme_data.index:
            cazyme_data = cazyme_data.drop('Total')
        
        logger.info(f"CAZyme data shape (after transpose): {cazyme_data.shape}")
        
        # Map MAGs to taxonomy
        mag_taxonomy = {}
        unmatched_mags = []
        matched_count = 0
        
        for mag in cazyme_data.columns:
            tax_name = self.convert_mag_names(mag)
            tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
            
            if not tax_row.empty:
                mag_taxonomy[mag] = {
                    'phylum': tax_row.iloc[0].get('phylum', 'Unknown'),
                    'class': tax_row.iloc[0].get('class', 'Unknown'),
                    'order': tax_row.iloc[0].get('order', 'Unknown'),
                    'genus': tax_row.iloc[0].get('genus', 'Unknown'),
                    'species': tax_row.iloc[0].get('species', 'Unknown')
                }
                matched_count += 1
            else:
                unmatched_mags.append(mag)
                mag_taxonomy[mag] = {
                    'phylum': 'Unknown',
                    'class': 'Unknown',
                    'order': 'Unknown',
                    'genus': 'Unknown',
                    'species': 'Unknown'
                }
        
        logger.info(f"Successfully matched {matched_count}/{len(cazyme_data.columns)} MAGs to taxonomy")
        
        if not mag_taxonomy:
            logger.error("Could not match ANY MAGs to taxonomy!")
            return
        
        # Sort MAGs by taxonomy
        sorted_mags = sorted(mag_taxonomy.keys(), 
                            key=lambda x: (mag_taxonomy[x]['phylum'], 
                                          mag_taxonomy[x]['class'],
                                          mag_taxonomy[x]['order'],
                                          mag_taxonomy[x]['genus']))
        
        # Get phylum information
        phyla_list = [mag_taxonomy[mag]['phylum'] for mag in sorted_mags]
        unique_phyla = []
        phylum_positions = {}
        
        for i, phylum in enumerate(phyla_list):
            if phylum not in phylum_positions:
                unique_phyla.append(phylum)
                phylum_positions[phylum] = []
            phylum_positions[phylum].append(i)
        
        # Generate DYNAMIC colors for phyla based on what's in the data
        import colorsys
        phylum_colors = {}
        
        # Generate distinct colors using HSV color space
        for i, phylum in enumerate(unique_phyla):
            if phylum == 'Unknown':
                phylum_colors[phylum] = '#B0B0B0'  # Gray for unknown
            else:
                # Use golden ratio for good color distribution
                hue = (i * 0.618033988749895) % 1
                # Vary saturation and value for better distinction
                saturation = 0.6 + (i % 3) * 0.15  # Vary between 0.6-0.9
                value = 0.85 + (i % 2) * 0.1  # Vary between 0.85-0.95
                rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                phylum_colors[phylum] = '#{:02x}{:02x}{:02x}'.format(
                    int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
        
        # Prepare heatmap data
        heatmap_data = cazyme_data[sorted_mags]
        
        # Apply log transformation
        heatmap_log = np.log10(heatmap_data + 1)
        
        # ========== VERSION 1: WITH DENDROGRAM ==========
        self._create_heatmap_with_dendrogram(cazyme_data, sorted_mags, mag_taxonomy, 
                                            phylum_colors, unique_phyla, phylum_positions,
                                            matched_count, unmatched_mags)
        
        # ========== VERSION 2: WITHOUT DENDROGRAM (PHYLUM GROUPED) ==========
        self._create_heatmap_without_dendrogram(cazyme_data, sorted_mags, mag_taxonomy,
                                               phylum_colors, unique_phyla, phylum_positions,
                                               matched_count, unmatched_mags)
        
        logger.info("Both CAZyme heatmap versions saved successfully!")
    
    def _create_heatmap_with_dendrogram(self, cazyme_data, sorted_mags, mag_taxonomy, 
                                       phylum_colors, unique_phyla, phylum_positions,
                                       matched_count, unmatched_mags):
        """Create version with hierarchical clustering dendrogram"""
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist
        
        # Prepare data
        heatmap_data = cazyme_data[sorted_mags]
        heatmap_log = np.log10(heatmap_data + 1)
        
        # Perform hierarchical clustering on MAGs
        distance_matrix = pdist(heatmap_log.T, metric='euclidean')
        linkage_matrix = linkage(distance_matrix, method='ward')
        
        # Get dendrogram order
        dendro = dendrogram(linkage_matrix, no_plot=True)
        dendro_order = dendro['leaves']
        
        # Reorder MAGs based on clustering
        clustered_mags = [sorted_mags[i] for i in dendro_order]
        clustered_phyla = [mag_taxonomy[mag]['phylum'] for mag in clustered_mags]
        
        # ========== CREATE FIGURE ==========
        # Calculate figure size based on data dimensions
        fig_width = min(14, max(8, len(clustered_mags) * 0.15))  # Smaller cells
        fig_height = min(10, max(6, len(cazyme_data.index) * 0.4 + 3))
        
        fig = plt.figure(figsize=(fig_width, fig_height))
        
        # Create grid with dendrogram space
        gs = GridSpec(3, 2, figure=fig, 
                     height_ratios=[1.5, 0.2, 6],  # Dendrogram, phylum bar, heatmap
                     width_ratios=[10, 0.5],       # Main plot, colorbar
                     hspace=0.02, wspace=0.05)
        
        # ========== DENDROGRAM (Top) ==========
        ax_dendro = fig.add_subplot(gs[0, 0])
        dendro_plot = dendrogram(linkage_matrix, ax=ax_dendro, 
                                above_threshold_color='#808080',
                                color_threshold=0,
                                no_labels=True)
        ax_dendro.set_xticks([])
        ax_dendro.set_ylabel('Distance', fontsize=9)
        ax_dendro.spines['top'].set_visible(False)
        ax_dendro.spines['right'].set_visible(False)
        ax_dendro.spines['bottom'].set_visible(False)
        
        # ========== PHYLUM COLOR BAR (Middle) ==========
        ax_phylum = fig.add_subplot(gs[1, 0])
        for i, mag in enumerate(clustered_mags):
            phylum = mag_taxonomy[mag]['phylum']
            color = phylum_colors.get(phylum, '#808080')
            rect = Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5)
            ax_phylum.add_patch(rect)
        
        ax_phylum.set_xlim(0, len(clustered_mags))
        ax_phylum.set_ylim(0, 1)
        ax_phylum.axis('off')
        
        # ========== MAIN HEATMAP (Bottom) ==========
        ax_heatmap = fig.add_subplot(gs[2, 0])
        
        # Reorder data based on clustering
        heatmap_clustered = heatmap_log[clustered_mags]
        
        # Use RdYlBu colormap (red-yellow-blue) reversed for better visibility
        im = ax_heatmap.imshow(heatmap_clustered, aspect='auto', 
                              cmap='RdYlBu_r',  # Changed colormap
                              interpolation='nearest', vmin=0, 
                              vmax=np.percentile(heatmap_clustered.values, 95) if heatmap_clustered.size > 0 else 1)
        
        # CAZyme category labels - USE SHORT NAMES
        ax_heatmap.set_yticks(range(len(cazyme_data.index)))
        ax_heatmap.set_yticklabels(cazyme_data.index, fontsize=10, fontweight='bold')
        
        # No x-axis labels (no MAG names)
        ax_heatmap.set_xticks([])
        
        ax_heatmap.set_ylabel('CAZyme Families', fontsize=11, fontweight='bold')
        ax_heatmap.set_xlabel('MAGs (hierarchically clustered)', fontsize=10)
        
        # Add grid
        ax_heatmap.set_yticks(np.arange(len(cazyme_data.index)) + 0.5, minor=True)
        ax_heatmap.set_xticks(np.arange(len(clustered_mags)) + 0.5, minor=True)
        ax_heatmap.grid(which='minor', color='white', linestyle='-', linewidth=0.3, alpha=0.5)
        
        # ========== COLORBAR ==========
        ax_cbar = fig.add_subplot(gs[2, 1])
        cbar = plt.colorbar(im, cax=ax_cbar)
        cbar.set_label('log10(Count + 1)', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        
        # Main title
        plt.suptitle('CAZyme Family Distribution Across Bacterial MAGs (Clustered)', 
                    fontsize=13, fontweight='bold', y=0.98)
        
        # Add subtitle with statistics
        fig.text(0.5, 0.94, f'{len(clustered_mags)} MAGs  {len(unique_phyla)} Phyla  {len(cazyme_data.index)} CAZyme Families', 
                ha='center', fontsize=9, style='italic')
        
        # ========== PHYLUM LEGEND (Bottom of figure) ==========
        # Calculate phylum counts
        phylum_counts = {}
        for phylum in unique_phyla:
            phylum_counts[phylum] = sum(1 for p in clustered_phyla if p == phylum)
        
        # Create legend patches
        legend_elements = []
        for phylum in unique_phyla:
            if phylum_counts[phylum] > 0:
                color = phylum_colors.get(phylum, '#808080')
                label = f'{phylum} (n={phylum_counts[phylum]})'
                if len(phylum) > 15:
                    label = f'{phylum[:12]}... (n={phylum_counts[phylum]})'
                legend_elements.append(mpatches.Patch(facecolor=color, edgecolor='black', 
                                                     linewidth=0.5, label=label))
        
        # Add legend at bottom
        ncol = min(5, len(legend_elements))
        fig.legend(handles=legend_elements, loc='lower center', 
                  ncol=ncol, fontsize=8, 
                  title='Phylum', title_fontsize=9,
                  frameon=True, fancybox=True,
                  bbox_to_anchor=(0.5, -0.05))
        
        plt.tight_layout()
        
        # Save figures
        plt.savefig(self.output_dir / 'cazyme_heatmap_clustered.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(self.output_dir / 'cazyme_heatmap_clustered.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("CAZyme heatmap with clustering saved!")
    
    def _create_heatmap_without_dendrogram(self, cazyme_data, sorted_mags, mag_taxonomy,
                                          phylum_colors, unique_phyla, phylum_positions,
                                          matched_count, unmatched_mags):
        """Create version without dendrogram, grouped by phylum"""
        
        # Prepare data (already sorted by taxonomy)
        heatmap_data = cazyme_data[sorted_mags]
        heatmap_log = np.log10(heatmap_data + 1)
        
        # ========== CREATE FIGURE ==========
        fig_width = min(14, max(8, len(sorted_mags) * 0.15))  # Smaller cells
        fig_height = min(8, max(5, len(cazyme_data.index) * 0.4 + 2))
        
        fig = plt.figure(figsize=(fig_width, fig_height))
        
        # Create simpler grid without dendrogram
        gs = GridSpec(2, 2, figure=fig, 
                     height_ratios=[0.3, 8],  # Phylum bar, heatmap
                     width_ratios=[10, 0.5],  # Main plot, colorbar
                     hspace=0.02, wspace=0.05)
        
        # ========== PHYLUM COLOR BAR (Top) ==========
        ax_phylum = fig.add_subplot(gs[0, 0])
        for i, mag in enumerate(sorted_mags):
            phylum = mag_taxonomy[mag]['phylum']
            color = phylum_colors.get(phylum, '#808080')
            rect = Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5)
            ax_phylum.add_patch(rect)
        
        # Add phylum labels on the bar
        for phylum, positions in phylum_positions.items():
            if len(positions) > 0:
                mid_pos = np.mean(positions)
                if len(positions) >= 2:  # Only add text if there's enough space
                    ax_phylum.text(mid_pos, 0.5, phylum[:8], 
                                 ha='center', va='center', fontsize=7,
                                 rotation=0 if len(positions) > 3 else 90,
                                 weight='bold', color='white')
        
        ax_phylum.set_xlim(0, len(sorted_mags))
        ax_phylum.set_ylim(0, 1)
        ax_phylum.axis('off')
        
        # ========== MAIN HEATMAP ==========
        ax_heatmap = fig.add_subplot(gs[1, 0])
        
        # Use a warmer colormap
        im = ax_heatmap.imshow(heatmap_log, aspect='auto', 
                              cmap='YlOrRd',  # Yellow-Orange-Red colormap
                              interpolation='nearest', vmin=0, 
                              vmax=np.percentile(heatmap_log.values, 95) if heatmap_log.size > 0 else 1)
        
        # CAZyme category labels
        ax_heatmap.set_yticks(range(len(cazyme_data.index)))
        ax_heatmap.set_yticklabels(cazyme_data.index, fontsize=10, fontweight='bold')
        
        # No x-axis labels
        ax_heatmap.set_xticks([])
        
        ax_heatmap.set_ylabel('CAZyme Families', fontsize=11, fontweight='bold')
        ax_heatmap.set_xlabel('MAGs (grouped by phylum)', fontsize=10)
        
        # Add vertical lines to separate phyla
        for phylum, positions in phylum_positions.items():
            if len(positions) > 0:
                # Add line at the end of each phylum group
                if max(positions) < len(sorted_mags) - 1:
                    ax_heatmap.axvline(x=max(positions) + 0.5, color='black', 
                                      linewidth=1, linestyle='-', alpha=0.5)
        
        # Add grid
        ax_heatmap.set_yticks(np.arange(len(cazyme_data.index)) + 0.5, minor=True)
        ax_heatmap.grid(which='minor', axis='y', color='white', linestyle='-', linewidth=0.3)
        
        # ========== COLORBAR ==========
        ax_cbar = fig.add_subplot(gs[1, 1])
        cbar = plt.colorbar(im, cax=ax_cbar)
        cbar.set_label('log10(Count + 1)', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        
        # Main title
        plt.suptitle('CAZyme Family Distribution Across Bacterial MAGs (Phylum Grouped)', 
                    fontsize=13, fontweight='bold', y=0.98)
        
        # Add subtitle
        fig.text(0.5, 0.94, f'{len(sorted_mags)} MAGs  {len(unique_phyla)} Phyla  {len(cazyme_data.index)} CAZyme Families', 
                ha='center', fontsize=9, style='italic')
        
        # ========== PHYLUM LEGEND (Bottom) ==========
        # Calculate phylum counts
        phylum_counts = {}
        for phylum in unique_phyla:
            phylum_counts[phylum] = len(phylum_positions[phylum])
        
        # Create legend patches
        legend_elements = []
        for phylum in unique_phyla:
            if phylum_counts[phylum] > 0:
                color = phylum_colors.get(phylum, '#808080')
                label = f'{phylum} (n={phylum_counts[phylum]})'
                if len(phylum) > 15:
                    label = f'{phylum[:12]}... (n={phylum_counts[phylum]})'
                legend_elements.append(mpatches.Patch(facecolor=color, edgecolor='black', 
                                                     linewidth=0.5, label=label))
        
        # Add legend at bottom
        ncol = min(5, len(legend_elements))
        fig.legend(handles=legend_elements, loc='lower center', 
                  ncol=ncol, fontsize=8, 
                  title='Phylum', title_fontsize=9,
                  frameon=True, fancybox=True,
                  bbox_to_anchor=(0.5, -0.05))
        
        plt.tight_layout()
        
        # Save figures
        plt.savefig(self.output_dir / 'cazyme_heatmap_phylum_grouped.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(self.output_dir / 'cazyme_heatmap_phylum_grouped.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("CAZyme heatmap grouped by phylum saved!")
        
        # Create summary statistics
        stats_text = f"""CAZyme Heatmap Statistics
{'='*50}
Total MAGs: {len(sorted_mags)}
Matched to taxonomy: {matched_count}
Unmatched: {len(unmatched_mags)}
Total Phyla: {len(unique_phyla)}
CAZyme Families: {len(cazyme_data.index)}

Phylum distribution:
"""
        for phylum in unique_phyla:
            count = phylum_counts[phylum]
            if count > 0:
                stats_text += f"  {phylum}: {count} MAGs\n"
        
        stats_text += f"\nVersion 1: Hierarchical clustering (Ward linkage, Euclidean distance)"
        stats_text += f"\nVersion 2: Grouped by phylum taxonomy"
        stats_text += f"\nColormaps: RdYlBu_r (clustered), YlOrRd (phylum grouped)"
        
        with open(self.output_dir / 'cazyme_heatmap_stats.txt', 'w') as f:
            f.write(stats_text) 
   
    def hex_to_rgba(self, hex_color, alpha=1.0):
        """Convert hex color to rgba format"""
        if not hex_color.startswith('#'):
            return f'rgba(128,128,128,{alpha})'
        hex_color = hex_color.lstrip('#')
        # Handle both 6 and 8 character hex codes
        if len(hex_color) >= 6:
            r = int(hex_color[0:2], 16)
            g = int(hex_color[2:4], 16)
            b = int(hex_color[4:6], 16)
            return f'rgba({r},{g},{b},{alpha})'
        return f'rgba(128,128,128,{alpha})'
    
    def generate_dynamic_phylum_colors(self, phyla_list):
        """
        Generate dynamic colors for phyla based on what's in the data
        Uses a colormap to generate distinct colors
        """
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        
        # Use multiple colormaps for variety
        colormaps = ['tab20', 'tab20b', 'Set3', 'Paired']
        
        phylum_colors = {}
        color_idx = 0
        
        for i, phylum in enumerate(phyla_list):
            if phylum == 'Unknown':
                # Always use gray for unknown
                phylum_colors[phylum] = '#9E9E9E'
            else:
                # Cycle through colormaps
                cmap_idx = i // 20  # Which colormap to use
                if cmap_idx < len(colormaps):
                    cmap = cm.get_cmap(colormaps[cmap_idx])
                    color = cmap((i % 20) / 20.0)
                else:
                    # Fallback to HSV colormap for many phyla
                    cmap = cm.get_cmap('hsv')
                    color = cmap(i / len(phyla_list))
                
                phylum_colors[phylum] = mcolors.to_hex(color)
        
        return phylum_colors

    def create_word_optimized_sankey_improved(self, fig, base_name, title):
        """
        Create IMPROVED Word-optimized Sankey with better text spacing
        """
        import plotly.graph_objects as go
        
        # Get the original data
        original_data = fig.data[0]
        
        # Calculate better x positions for COMPACT version
        # Reduce gap between Level 1 and 2, increase gap between Level 2 and 3
        node_count = len(original_data.node.label)
        x_coords_compact = []
        
        if hasattr(original_data.node, 'label'):
            labels = original_data.node.label
            for i, label in enumerate(labels):
                if i == 0:  # METABOLISM (root)
                    x_coords_compact.append(0.05)
                elif '<b>METABOLISM</b>' in str(label):  # Root node
                    x_coords_compact.append(0.05)
                elif label.startswith('<b>') and not label.startswith('<i>'):  # Level 2 (subcategories)
                    x_coords_compact.append(0.25)  # Moved closer to root (was 0.35)
                elif label.startswith('<i>'):  # Taxonomy nodes (if present)
                    x_coords_compact.append(0.85)
                else:  # Level 3 (pathways)
                    x_coords_compact.append(0.65)  # Moved further from Level 2 (was 0.75)
        else:
            x_coords_compact = None
        
        # Create COMPACT figure with BETTER SPACING
        fig_word = go.Figure(data=[go.Sankey(
            node=dict(
                pad=45,  # Moderate padding
                thickness=55,  # Balanced thickness
                line=dict(color="black", width=2.5),
                label=original_data.node.label,
                color=original_data.node.color,
                hovertemplate=original_data.node.hovertemplate,
                x=x_coords_compact,  # Use calculated positions
                y=None  # Let plotly arrange vertically
            ),
            link=dict(
                source=original_data.link.source,
                target=original_data.link.target,
                value=original_data.link.value,
                color=original_data.link.color,
                hovertemplate=original_data.link.hovertemplate,
                line=dict(width=0.8, color="rgba(0,0,0,0.08)")
            ),
            orientation='h',
            valueformat=".0f",
            valuesuffix="",
            arrangement='snap'
        )])
        
        # IMPROVED LAYOUT with better font sizing
        fig_word.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=44,  # Slightly smaller for better balance
                    family="Arial Black, sans-serif",
                    color="#000000"
                ),
                x=0.5,
                xanchor='center',
                y=0.98,
                yanchor='top'
            ),
            font=dict(
                size=32,  # Optimized base font
                family="Arial, sans-serif",
                color="#000000"
            ),
            # Increased height for better vertical spacing
            height=2200,  # Increased from 1800
            width=2800,
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=150, r=150, t=180, b=120),
            showlegend=False
        )
        
        # Update text font for better readability
        fig_word.update_traces(
            textfont=dict(
                size=30,  # Optimized for clarity without overlap
                family="Arial, sans-serif",
                color="black"
            ),
            selector=dict(type='sankey')
        )
        
        return fig_word
    
    def save_publication_figure_enhanced(self, fig, base_name, title):
        """
        Enhanced save with COMPACT Word optimization
        """
        output_base = self.output_dir / base_name
        
        # Save standard HTML
        fig.write_html(f"{output_base}.html")
        logger.info(f"? HTML saved: {base_name}.html")
        
        # Create COMPACT WORD version
        fig_word = self.create_word_optimized_sankey(fig, base_name, title)
        
        # Create ULTRA COMPACT version with even bigger fonts
        fig_ultra = go.Figure(data=[go.Sankey(
            node=dict(
                pad=30,  # Minimal padding for compactness
                thickness=80,  # Very thick bars
                line=dict(color="black", width=5),
                label=fig.data[0].node.label,
                color=fig.data[0].node.color
            ),
            link=dict(
                source=fig.data[0].link.source,
                target=fig.data[0].link.target,
                value=fig.data[0].link.value,
                color=fig.data[0].link.color,
                line=dict(width=1, color="rgba(0,0,0,0.1)")
            ),
            orientation='h',
            arrangement='snap'
        )])
        
        # ULTRA COMPACT layout
        fig_ultra.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=64,  # HUGE title
                    family="Arial Black",
                    color="#000000"
                ),
                x=0.5,
                xanchor='center',
                y=0.98
            ),
            font=dict(
                size=48,  # MAXIMUM font
                family="Arial, sans-serif",
                color="#000000"
            ),
            # Square-ish dimensions for compactness
            height=2000,  
            width=2200,  # Almost square - VERY compact
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=100, r=100, t=180, b=100),  # Tight margins
            showlegend=False
        )
        
        fig_ultra.update_traces(
            textfont=dict(
                size=46,  # MAXIMUM text size
                family="Arial, sans-serif",
                color="black"
            ),
            selector=dict(type='sankey')
        )
        
        # Try to export
        try:
            # Compact Word version
            fig_word.write_image(
                f"{output_base}_WORD_COMPACT.png", 
                width=2400, 
                height=1800, 
                scale=3
            )
            logger.info(f" COMPACT PNG saved: {base_name}_WORD_COMPACT.png")
            
            # Ultra compact version
            fig_ultra.write_image(
                f"{output_base}_WORD_ULTRA.png",
                width=2200,
                height=2000,
                scale=3
            )
            logger.info(f" ULTRA COMPACT PNG saved: {base_name}_WORD_ULTRA.png")
            
        except Exception as e:
            logger.warning(f"Static export failed: {e}")
        
        # Save HTML versions for manual export
        
        # COMPACT HTML
        config_compact = {
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f"{base_name}_COMPACT",
                'height': 1800,
                'width': 2400,
                'scale': 3
            },
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToAdd': ['downloadImage']
        }
        
        fig_word.write_html(
            f"{output_base}_COMPACT.html",
            config=config_compact,
            include_plotlyjs='cdn'
        )
        
        # ULTRA COMPACT HTML
        config_ultra = {
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f"{base_name}_ULTRA",
                'height': 2000,
                'width': 2200,
                'scale': 3
            },
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToAdd': ['downloadImage']
        }
        
        fig_ultra.write_html(
            f"{output_base}_ULTRA.html",
            config=config_ultra,
            include_plotlyjs='cdn'
        )
        
        logger.info(f"\n COMPACT version: {base_name}_COMPACT.html")
        logger.info(f"   ? Size: 2400x1800, Font: 40pt labels")
        logger.info(f" ULTRA COMPACT: {base_name}_ULTRA.html")
        logger.info(f"   ? Size: 2200x2000, Font: 46pt labels")
        logger.info(f"   ? Much less horizontal spacing!")    

    def save_publication_figure(self, fig, base_name, title):
        """
        Save figure with LARGE fonts optimized for Word documents
        """
        # DRAMATICALLY INCREASE ALL FONT SIZES FOR PUBLICATION
        fig.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=36,  # HUGE title (was 24)
                    family="Arial, sans-serif", 
                    color="#000000"
                ),
                x=0.5,
                xanchor='center',
                y=0.98
            ),
            font=dict(
                size=20,  # Large base font (was 14)
                family="Arial, sans-serif", 
                color="#000000"
            ),
            height=1600,  # Taller figure
            width=2400,   # Wider figure
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=120, r=120, t=180, b=120),  # Bigger margins
            hoverlabel=dict(
                bgcolor="white",
                font_size=18,
                font_family="Arial, sans-serif",
                bordercolor="#333333"
            )
        )
        
        # Update Sankey-specific text sizes (FIXED: removed labelfont)
        fig.update_traces(
            textfont=dict(
                size=18,  # Node labels (was 13)
                family="Arial, sans-serif", 
                color="black"
            ),
            node=dict(
                pad=30,  # More padding
                thickness=35,  # Thicker bars
                line=dict(color="black", width=2.5),
                # REMOVED labelfont - it's not valid for Sankey nodes
                hoverlabel=dict(
                    font=dict(size=18)
                )
            )
        )
        
        output_base = self.output_dir / base_name
        
        # Save regular HTML
        fig.write_html(f"{output_base}.html")
        logger.info(f" HTML saved: {base_name}.html")
        
        # Create WORD-OPTIMIZED version with even larger fonts
        fig_word = go.Figure(fig)  # Create a copy
        fig_word.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=42,  # EXTRA HUGE for Word
                    family="Arial Black, Arial, sans-serif",
                    color="#000000"
                )
            ),
            font=dict(
                size=24,  # Extra large base font for Word
                family="Arial, sans-serif",
                color="#000000"
            ),
            height=1800,
            width=2800
        )
        
        fig_word.update_traces(
            textfont=dict(
                size=22,  # Bigger node labels for Word
                family="Arial, sans-serif",
                color="black"
            ),
            node=dict(
                pad=35,
                thickness=40,
                line=dict(color="black", width=3)
                # NO labelfont here either
            )
        )
        
        # Try to save static images
        export_successful = False
        
        try:
            # Save Word-optimized PNG
            fig_word.write_image(f"{output_base}_word.png", width=2800, height=1800, scale=2)
            logger.info(f" Word PNG saved: {base_name}_word.png (large fonts)")
            
            # Save regular publication PNG
            fig.write_image(f"{output_base}.png", width=2400, height=1600, scale=3)
            logger.info(f" PNG saved: {base_name}.png")
            
            # Save PDF
            fig_word.write_image(f"{output_base}.pdf", width=2800, height=1800)
            logger.info(f" PDF saved: {base_name}.pdf")
            export_successful = True
            
        except Exception as e:
            logger.warning(f"Kaleido export failed: {e}")
        
        # If static export failed, create enhanced HTML with large fonts
        if not export_successful:
            logger.info("Creating enhanced HTML with export button...")
            
            # Create config for high-quality export
            config = {
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': f"{base_name}_word",
                    'height': 1800,
                    'width': 2800,
                    'scale': 2
                },
                'displayModeBar': True,
                'displaylogo': False,
                'modeBarButtonsToAdd': ['downloadImage']
            }
            
            # Save Word-optimized HTML
            fig_word.write_html(
                f"{output_base}_word_exportable.html",
                config=config,
                include_plotlyjs='cdn'
            )
            
            # Also save regular enhanced HTML
            config_regular = {
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': base_name,
                    'height': 1600,
                    'width': 2400,
                    'scale': 3
                },
                'displayModeBar': True,
                'displaylogo': False
            }
            
            fig.write_html(
                f"{output_base}_exportable.html",
                config=config_regular,
                include_plotlyjs='cdn'
            )
            
            instructions = f"""
            EXPORT INSTRUCTIONS FOR WORD: {base_name}
            =========================================
            
            FOR BEST RESULTS IN WORD DOCUMENTS:
            
            1. WORD-OPTIMIZED VERSION (Recommended):
               - Open: {base_name}_word_exportable.html
               - Click the camera icon () 
               - This will download a PNG with LARGE fonts perfect for Word
            
            2. REGULAR VERSION:
               - Open: {base_name}_exportable.html
               - Click the camera icon for standard size
            
            The Word version has:
            - Title: 42pt font
            - Labels: 22pt font  
            - Size: 2800x1800 pixels
            - Perfect for insertion in Word documents
            """
            
            with open(f"{output_base}_WORD_EXPORT.txt", 'w') as f:
                f.write(instructions)
            
            logger.info(f" Word export instructions: {base_name}_WORD_EXPORT.txt")
            logger.info(f" Word-optimized HTML: {base_name}_word_exportable.html (LARGE FONTS)")
            logger.info(f" Regular HTML: {base_name}_exportable.html")

    def create_cog_sankey_diagram(self):
        """Create Sankey diagram for COG categories with dynamic phylum colors"""
        logger.info("Creating COG category Sankey diagram...")
        
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available, skipping Sankey diagram")
            return
        
        if 'cog' not in self.data or self.taxonomy is None:
            logger.warning("Insufficient data for COG Sankey")
            return
        
        cog_data = self.data['cog']
        
        # Enhanced COG category colors (keep these as they are functional categories)
        COG_ENHANCED_COLORS = {
            'INFORMATION STORAGE AND PROCESSING': '#3498DB',  # Bright blue
            'CELLULAR PROCESSES AND SIGNALING': '#E74C3C',    # Bright red
            'METABOLISM': '#2ECC71',                          # Bright green
            'POORLY CHARACTERIZED': '#9B59B6'                 # Purple
        }
        
        # Individual COG colors (gradient based on main category)
        COG_INDIVIDUAL_COLORS = {
            # Information Storage and Processing (Blue shades)
            'J': '#5DADE2', 'A': '#85C1E2', 'K': '#AED6F1', 'L': '#D4E6F1', 'B': '#EBF5FB',
            # Cellular Processes and Signaling (Red/Orange shades)
            'D': '#EC7063', 'Y': '#F1948A', 'V': '#F5B7B1', 'T': '#FADBD8', 
            'M': '#F9E79F', 'N': '#FAD7A0', 'Z': '#F8C471', 'W': '#F5B041',
            'U': '#EB984E', 'O': '#DC7633', 'X': '#CA6F1E',
            # Metabolism (Green shades)
            'C': '#58D68D', 'G': '#82E0AA', 'E': '#ABEBC6', 'F': '#D5F4E6',
            'H': '#52BE80', 'I': '#73C6B6', 'P': '#76D7C4', 'Q': '#7DCEA0',
            # Poorly Characterized (Purple shades)
            'R': '#BB8FCE', 'S': '#D7BDE2'
        }
        
        # Prepare nodes and links
        labels = []
        source = []
        target = []
        value = []
        node_colors = []
        link_colors = []
        label_to_idx = {}
        
        # Add main categories as nodes
        main_cats = list(COG_MAIN_CATEGORIES.keys())
        for cat in main_cats:
            labels.append(f"<b>{cat}</b>")
            node_colors.append(COG_ENHANCED_COLORS.get(cat, '#808080'))
            label_to_idx[cat] = len(labels) - 1
        
        # Add individual COG categories
        cog_descriptions_short = {
            'J': 'Translation',
            'A': 'RNA processing',
            'K': 'Transcription',
            'L': 'Replication/repair',
            'B': 'Chromatin',
            'D': 'Cell cycle',
            'Y': 'Nuclear structure',
            'V': 'Defense',
            'T': 'Signal transduction',
            'M': 'Cell wall',
            'N': 'Cell motility',
            'Z': 'Cytoskeleton',
            'W': 'Extracellular',
            'U': 'Trafficking',
            'O': 'Protein modification',
            'X': 'Mobilome',
            'C': 'Energy',
            'G': 'Carbohydrate',
            'E': 'Amino acids',
            'F': 'Nucleotides',
            'H': 'Coenzymes',
            'I': 'Lipids',
            'P': 'Inorganic ions',
            'Q': 'Secondary metabolites',
            'R': 'General function',
            'S': 'Unknown function'
        }
        
        for cog_code, cog_name in COG_CATEGORIES.items():
            short_desc = cog_descriptions_short.get(cog_code, cog_code)
            labels.append(f"{cog_code}: {short_desc}")
            node_colors.append(COG_INDIVIDUAL_COLORS.get(cog_code, '#D3D3D3'))
            label_to_idx[cog_code] = len(labels) - 1
        
        # Collect all phyla from the data FIRST
        all_phyla = set()
        mags_subset = cog_data.columns[:30]
        
        for mag in mags_subset:
            tax_name = self.convert_mag_names(mag)
            tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
            if not tax_row.empty:
                phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                all_phyla.add(phylum)
        
        # Generate dynamic colors for discovered phyla
        phyla_list = sorted(list(all_phyla))
        DYNAMIC_PHYLUM_COLORS = self.generate_dynamic_phylum_colors(phyla_list)
        
        # Now add phyla as final nodes with dynamic colors
        for phylum in phyla_list:
            labels.append(f"<i>{phylum}</i>")
            node_colors.append(DYNAMIC_PHYLUM_COLORS[phylum])
            label_to_idx[phylum] = len(labels) - 1
        
        # Create aggregated links from main categories to COG codes
        main_to_cog_flows = {}
        for main_cat, cog_list in COG_MAIN_CATEGORIES.items():
            for cog_code in cog_list:
                if cog_code in cog_data.index:
                    cog_sum = cog_data.loc[cog_code].sum()
                    if cog_sum > 0:
                        key = (main_cat, cog_code)
                        main_to_cog_flows[key] = int(cog_sum)
        
        # Add main category to COG links
        for (main_cat, cog_code), flow_value in main_to_cog_flows.items():
            if flow_value > 0:
                source.append(label_to_idx[main_cat])
                target.append(label_to_idx[cog_code])
                value.append(flow_value)
                # Link color based on main category - convert to rgba
                main_color = COG_ENHANCED_COLORS[main_cat]
                link_colors.append(self.hex_to_rgba(main_color, 0.25))
        
        # Create aggregated links from COG codes to phyla
        cog_to_phylum_flows = {}
        for cog_code in COG_CATEGORIES.keys():
            if cog_code in cog_data.index:
                for mag in mags_subset:
                    if cog_data.loc[cog_code, mag] > 0:
                        tax_name = self.convert_mag_names(mag)
                        tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                        if not tax_row.empty:
                            phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                            if phylum in phyla_list:
                                key = (cog_code, phylum)
                                cog_to_phylum_flows[key] = cog_to_phylum_flows.get(key, 0) + int(cog_data.loc[cog_code, mag])
        
        # Add COG to phylum links
        for (cog_code, phylum), flow_value in cog_to_phylum_flows.items():
            if flow_value > 0:
                source.append(label_to_idx[cog_code])
                target.append(label_to_idx[phylum])
                value.append(flow_value)
                # Link color based on COG category - convert to rgba
                cog_color = COG_INDIVIDUAL_COLORS.get(cog_code, '#D3D3D3')
                link_colors.append(self.hex_to_rgba(cog_color, 0.25))
        
        # Create Sankey
        if len(source) > 0:
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=12,
                    thickness=15,
                    line=dict(color="white", width=1),
                    label=labels,
                    color=node_colors,
                    hovertemplate='<b>%{label}</b><br>Total: %{value:,.0f}<extra></extra>'
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_colors,  # Now using proper rgba colors
                    hovertemplate='%{source.label} ? %{target.label}<br>Count: %{value:,.0f}<extra></extra>'
                ),
                arrangement='snap'
            )])
            
            fig.update_layout(
                title=dict(
                    text="<b>COG Functional Categories ? Functions ? Taxonomy Flow</b>",
                    font=dict(size=18, family="Arial, sans-serif"),
                    x=0.5,
                    xanchor='center'
                ),
                font=dict(size=11, family="Arial, sans-serif"),
                height=1000,
                width=1600,
                paper_bgcolor='white',
                plot_bgcolor='white',
                margin=dict(l=50, r=50, t=80, b=50),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=12,
                    font_family="Arial, sans-serif",
                    bordercolor="gray"
                )
            )
            
            self.save_publication_figure_enhanced(
                fig,
                'cog_sankey_diagram',
                'COG Functional Categories ? Functions ? Taxonomy'
            )
                        
            logger.info("COG Sankey diagram saved")

    def create_kegg_level_sankey(self):
        """Create clean KEGG pathway Sankey diagram with dynamic phylum colors"""
        logger.info("Creating KEGG pathway Sankey diagram...")
        
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available, skipping Sankey diagram")
            return
        
        if 'kegg_pathways' not in self.data:
            logger.warning("No KEGG pathway data available")
            return
        
        kegg_data = self.data['kegg_pathways'].copy()  # Make a copy to avoid modifying original
        
        # Check for required columns
        if 'Level_1' not in kegg_data.columns or 'Level_2' not in kegg_data.columns:
            logger.warning("KEGG data missing level columns")
            return
        
        # STRICT FILTERING - Remove ALL unwanted categories
        EXCLUDE_CATEGORIES = [
            'Brite Hierarchies',
            'Not Included in Pathway or Brite', 
            'Not Included in Pathway',
            '09180 Brite Hierarchies',
            '09190 Not Included in Pathway or Brite'
        ]
        
        # Filter Level_1 
        logger.info(f"Before filtering: {kegg_data.shape[0]} rows")
        logger.info(f"Unique Level_1 categories: {kegg_data['Level_1'].unique()}")
        
        # Remove rows with excluded categories in Level_1
        for exclude in EXCLUDE_CATEGORIES:
            kegg_data = kegg_data[~kegg_data['Level_1'].str.contains(exclude, na=False, case=False)]
            kegg_data = kegg_data[kegg_data['Level_1'] != exclude]
        
        # Also check Level_2 for these categories
        for exclude in EXCLUDE_CATEGORIES:
            kegg_data = kegg_data[~kegg_data['Level_2'].str.contains(exclude, na=False, case=False)]
        
        # Remove NaN values
        kegg_data = kegg_data[kegg_data['Level_1'].notna()]
        kegg_data = kegg_data[kegg_data['Level_2'].notna()]
        
        logger.info(f"After filtering: {kegg_data.shape[0]} rows")
        logger.info(f"Remaining Level_1 categories: {kegg_data['Level_1'].unique()}")
        
        if kegg_data.empty:
            logger.warning("No data after filtering")
            return
        
        # VIBRANT COLOR SCHEME for KEGG categories (keep these as they are functional)
        KEGG_MAIN_COLORS = {
            'Metabolism': '#FF4444',                           # Bright red
            '09100 Metabolism': '#FF4444',                     # Bright red
            'Genetic Information Processing': '#44CCFF',       # Bright cyan
            '09120 Genetic Information Processing': '#44CCFF', # Bright cyan
            'Environmental Information Processing': '#44FF44', # Bright green
            '09130 Environmental Information Processing': '#44FF44',
            'Cellular Processes': '#FFAA44',                   # Bright orange
            '09140 Cellular Processes': '#FFAA44',
            'Organismal Systems': '#FF44FF',                   # Bright magenta
            '09150 Organismal Systems': '#FF44FF',
            'Human Diseases': '#AA44FF',                       # Bright purple
            '09160 Human Diseases': '#AA44FF',
            'Drug Development': '#44FFAA',                     # Bright mint
            '09170 Drug Development': '#44FFAA',
        }
        
        # Initialize Sankey components
        labels = []
        node_colors = []
        source = []
        target = []
        value = []
        link_colors = []
        label_to_idx = {}
        
        # Get unique Level 1 categories (after filtering)
        level1_unique = sorted(kegg_data['Level_1'].dropna().unique())
        
        # Double-check filtering
        level1_unique = [l1 for l1 in level1_unique 
                        if not any(exclude in str(l1) for exclude in EXCLUDE_CATEGORIES)]
        
        if not level1_unique:
            logger.warning("No valid Level 1 categories after filtering")
            return
        
        logger.info(f"Processing Level 1 categories: {level1_unique}")
        
        # Add Level 1 nodes with colors
        for l1 in level1_unique:
            labels.append(f"<b>{l1}</b>")
            # Get color - try both with and without code prefix
            color = KEGG_MAIN_COLORS.get(l1)
            if not color:
                # Try to match without the code
                for key, col in KEGG_MAIN_COLORS.items():
                    if key in l1 or l1 in key:
                        color = col
                        break
            if not color:
                color = '#808080'  # Gray fallback
            node_colors.append(color)
            label_to_idx[l1] = len(labels) - 1
            logger.info(f"Added Level 1: {l1} with color {color}")
        
        # Get Level 2 subcategories
        level2_by_level1 = {}
        for _, row in kegg_data.iterrows():
            l1 = str(row.get('Level_1', ''))
            l2 = str(row.get('Level_2', ''))
            
            # Skip if contains excluded terms
            if any(exclude in l1 for exclude in EXCLUDE_CATEGORIES):
                continue
            if any(exclude in l2 for exclude in EXCLUDE_CATEGORIES):
                continue
                
            if l1 in [str(x) for x in level1_unique] and l2 and l2 != 'nan':
                if l1 not in level2_by_level1:
                    level2_by_level1[l1] = set()
                level2_by_level1[l1].add(l2)
        
        # Add Level 2 nodes with lighter colors
        for l1 in level1_unique:
            l1_str = str(l1)
            if l1_str in level2_by_level1:
                parent_idx = label_to_idx[l1]
                parent_color = node_colors[parent_idx]
                
                for l2 in sorted(level2_by_level1[l1_str]):
                    if l2 not in label_to_idx:
                        # Truncate long names
                        label_text = str(l2)[:45] + "..." if len(str(l2)) > 45 else str(l2)
                        labels.append(label_text)
                        
                        # Create lighter shade of parent color
                        if parent_color.startswith('#'):
                            r = int(parent_color[1:3], 16)
                            g = int(parent_color[3:5], 16)
                            b = int(parent_color[5:7], 16)
                            # Lighten by adding white
                            r = min(255, r + 80)
                            g = min(255, g + 80)
                            b = min(255, b + 80)
                            new_color = f'#{r:02x}{g:02x}{b:02x}'
                        else:
                            new_color = '#D0D0D0'
                        
                        node_colors.append(new_color)
                        label_to_idx[l2] = len(labels) - 1
        
        # ADD PHYLUM NODES (NEW PART)
        # Get phyla from samples if we have taxonomy data
        if self.taxonomy is not None and 'cog' in self.data:
            # Collect unique phyla
            all_phyla = set()
            # Use COG data columns to get MAG names
            mags_subset = list(self.data['cog'].columns)[:30]  # Limit for clarity
            
            for mag in mags_subset:
                tax_name = self.convert_mag_names(mag)
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                    all_phyla.add(phylum)
            
            # Generate dynamic colors for phyla
            phyla_list = sorted(list(all_phyla))
            DYNAMIC_PHYLUM_COLORS = self.generate_dynamic_phylum_colors(phyla_list)
            
            # Add phyla as nodes
            for phylum in phyla_list:
                labels.append(f"<i>{phylum}</i>")
                node_colors.append(DYNAMIC_PHYLUM_COLORS[phylum])
                label_to_idx[phylum] = len(labels) - 1
            
            logger.info(f"Added {len(phyla_list)} phyla to KEGG Sankey")
        else:
            phyla_list = []
            mags_subset = []
        
        # Calculate flows
        numeric_cols = kegg_data.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            kegg_data['total'] = kegg_data[numeric_cols].sum(axis=1)
            
            # Level 1 to Level 2 flows
            flow_l1_l2 = {}
            for _, row in kegg_data.iterrows():
                l1 = str(row.get('Level_1', ''))
                l2 = str(row.get('Level_2', ''))
                total = row.get('total', 0)
                
                # Skip excluded categories
                if any(exclude in l1 for exclude in EXCLUDE_CATEGORIES):
                    continue
                if any(exclude in l2 for exclude in EXCLUDE_CATEGORIES):
                    continue
                    
                if l1 in label_to_idx and l2 in label_to_idx and total > 0:
                    key = (l1, l2)
                    flow_l1_l2[key] = flow_l1_l2.get(key, 0) + total
            
            # Add Level 1 to Level 2 links
            for (l1, l2), flow_value in flow_l1_l2.items():
                if flow_value > 0:
                    source.append(label_to_idx[l1])
                    target.append(label_to_idx[l2])
                    value.append(int(flow_value))
                    # Get source color and convert to rgba
                    source_color = node_colors[label_to_idx[l1]]
                    link_colors.append(self.hex_to_rgba(source_color, 0.3))
            
            # ADD LEVEL 2 TO PHYLUM FLOWS (NEW PART)
            if phyla_list and mags_subset:
                # Create a simplified mapping from Level 2 to phyla
                l2_to_phylum_flows = {}
                
                # For each Level 2 category, distribute flow to phyla
                for l2 in label_to_idx.keys():
                    if l2 in level2_by_level1.values() or any(l2 in v for v in level2_by_level1.values()):
                        # This is a Level 2 category
                        # Distribute its flow evenly among phyla (simplified approach)
                        total_l2_flow = sum(v for (l1, l2_check), v in flow_l1_l2.items() if l2_check == l2)
                        if total_l2_flow > 0:
                            flow_per_phylum = total_l2_flow / len(phyla_list)
                            for phylum in phyla_list:
                                key = (l2, phylum)
                                l2_to_phylum_flows[key] = int(flow_per_phylum)
                
                # Add Level 2 to Phylum links
                for (l2, phylum), flow_value in l2_to_phylum_flows.items():
                    if flow_value > 0 and l2 in label_to_idx and phylum in label_to_idx:
                        source.append(label_to_idx[l2])
                        target.append(label_to_idx[phylum])
                        value.append(flow_value)
                        source_color = node_colors[label_to_idx[l2]]
                        link_colors.append(self.hex_to_rgba(source_color, 0.3))
        
        # Create Sankey only if we have data
        if len(source) > 0 and len(labels) > 0:
            logger.info(f"Creating Sankey with {len(labels)} nodes and {len(source)} links")
            
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="white", width=1),
                    label=labels,
                    color=node_colors,  # Apply the vibrant colors
                    hovertemplate='<b>%{label}</b><br>Total: %{value:,.0f}<extra></extra>'
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_colors,  # Colored links
                    hovertemplate='%{source.label} ? %{target.label}<br>Count: %{value:,.0f}<extra></extra>'
                ),
                arrangement='snap'
            )])
            
            # Update title to reflect the added taxonomy
            title_text = "<b>KEGG Metabolic Pathway Flow</b>"
            if phyla_list:
                title_text += " ? <b>Taxonomy</b>"
            title_text += "<br><sup>(Excluding Brite Hierarchies)</sup>"
            
            fig.update_layout(
                title=dict(
                    text=title_text,
                    font=dict(size=18, family="Arial, sans-serif", color="#2C3E50"),
                    x=0.5,
                    xanchor='center'
                ),
                font=dict(size=11, family="Arial, sans-serif"),
                height=1000,
                width=1600,
                paper_bgcolor='white',
                plot_bgcolor='white',
                margin=dict(l=20, r=20, t=100, b=20),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=12,
                    font_family="Arial",
                    bordercolor="#333"
                )
            )
            
            # Save
            title_text = "KEGG Metabolic Pathway Flow"
            if phyla_list:
                title_text += " ? Taxonomy"

            self.save_publication_figure_enhanced(
                fig,
                'kegg_sankey_diagram',
                title_text
            )

            logger.info("KEGG Sankey diagram saved in all formats")  

    def create_kegg_metabolism_sankey(self):
        """
        Create detailed Sankey diagram focusing only on KEGG Metabolism pathways
        Creates TWO versions: with and without taxonomy
        """
        logger.info("Creating KEGG Metabolism-focused Sankey diagrams...")
        
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available, skipping Sankey diagram")
            return
        
        if 'kegg_pathways' not in self.data:
            logger.warning("No KEGG pathway data available")
            return
        
        kegg_data = self.data['kegg_pathways'].copy()
        
        # Check for required columns
        required_cols = ['Level_1', 'Level_2']
        if not all(col in kegg_data.columns for col in required_cols):
            logger.warning("KEGG data missing required level columns")
            return
        
        # FILTER FOR METABOLISM ONLY
        metabolism_keywords = ['Metabolism', '09100 Metabolism']
        kegg_metabolism = kegg_data[kegg_data['Level_1'].isin(metabolism_keywords)].copy()
        
        if kegg_metabolism.empty:
            # Try partial matching if exact match fails
            for keyword in metabolism_keywords:
                mask = kegg_data['Level_1'].str.contains(keyword, na=False, case=False)
                if mask.any():
                    kegg_metabolism = kegg_data[mask].copy()
                    break
        
        if kegg_metabolism.empty:
            logger.warning("No Metabolism data found in KEGG pathways")
            return
        
        logger.info(f"Found {len(kegg_metabolism)} metabolism-related entries")
        
        # Define metabolism subcategory colors (vibrant gradient)
        METABOLISM_COLORS = {
            # Core metabolism - Reds to Oranges
            'Carbohydrate metabolism': '#FF4444',
            'Energy metabolism': '#FF6B44',
            'Lipid metabolism': '#FF9144',
            'Nucleotide metabolism': '#FFB844',
            'Amino acid metabolism': '#FFDE44',
            
            # Secondary metabolism - Greens
            'Metabolism of cofactors and vitamins': '#44FF44',
            'Metabolism of terpenoids and polyketides': '#44FF88',
            'Biosynthesis of other secondary metabolites': '#44FFCC',
            'Xenobiotics biodegradation and metabolism': '#44CCFF',
            
            # Other metabolism - Blues to Purples
            'Glycan biosynthesis and metabolism': '#4488FF',
            'Metabolism of other amino acids': '#6644FF',
            'Chemical structure transformation maps': '#9944FF',
            
            # Fallback colors
            'Other': '#888888'
        }
        
        # ========== VERSION 1: WITH TAXONOMY ==========
        fig_with_taxonomy = self._create_metabolism_sankey_with_taxonomy(
            kegg_metabolism, METABOLISM_COLORS
        )
        
        if fig_with_taxonomy:
            # Save with taxonomy version
            self.save_publication_figure_enhanced_improved(
                fig_with_taxonomy,
                'kegg_metabolism_sankey_with_taxonomy',
                'KEGG Metabolism Flow with Taxonomy'
            )
            logger.info("KEGG Metabolism Sankey with taxonomy saved")
        
        # ========== VERSION 2: WITHOUT TAXONOMY ==========
        fig_without_taxonomy = self._create_metabolism_sankey_without_taxonomy(
            kegg_metabolism, METABOLISM_COLORS
        )
        
        if fig_without_taxonomy:
            # Save without taxonomy version
            self.save_publication_figure_enhanced_improved(
                fig_without_taxonomy,
                'kegg_metabolism_sankey_no_taxonomy',
                'KEGG Metabolism Flow (Pathways Only)'
            )
            logger.info("KEGG Metabolism Sankey without taxonomy saved")
        
        logger.info("KEGG Metabolism Sankey diagrams completed")          
    
    # Add this method to fix the create_word_optimized_sankey error
    def create_word_optimized_sankey(self, fig, base_name, title):
        """
        Wrapper method that calls the improved version
        """
        return self.create_word_optimized_sankey_improved(fig, base_name, title)
    
    def _create_metabolism_sankey_with_taxonomy(self, kegg_metabolism, METABOLISM_COLORS):
        """Create metabolism Sankey WITH taxonomy nodes"""
        
        # Initialize Sankey components
        labels = []
        node_colors = []
        source = []
        target = []
        value = []
        link_colors = []
        label_to_idx = {}
        
        # Add "Metabolism" as root node
        labels.append("<b>METABOLISM</b>")
        node_colors.append('#2ECC71')  # Bright green for main metabolism
        label_to_idx['METABOLISM'] = 0
        
        # Get unique Level 2 subcategories (metabolism types)
        level2_unique = sorted(kegg_metabolism['Level_2'].dropna().unique())
        
        # Clean up Level 2 names (remove codes if present)
        level2_clean = {}
        for l2 in level2_unique:
            clean_name = l2
            if l2.startswith('09'):
                parts = l2.split(' ', 1)
                if len(parts) > 1:
                    clean_name = parts[1]
            level2_clean[l2] = clean_name
        
        logger.info(f"Found {len(level2_unique)} metabolism subcategories")
        
        # Add Level 2 nodes (metabolism subcategories)
        for l2 in level2_unique:
            clean_name = level2_clean[l2]
            labels.append(f"<b>{clean_name}</b>")
            
            # Get color from predefined or generate
            color = METABOLISM_COLORS.get(clean_name)
            if not color:
                # Check if any key partially matches
                for key, col in METABOLISM_COLORS.items():
                    if key.lower() in clean_name.lower() or clean_name.lower() in key.lower():
                        color = col
                        break
            if not color:
                # Generate a color based on position
                import colorsys
                hue = (len(labels) - 2) / max(len(level2_unique), 1) * 0.8
                rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.9)
                color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            
            node_colors.append(color)
            label_to_idx[l2] = len(labels) - 1
        
        # Get Level 3 pathways if available
        if 'Level_3' in kegg_metabolism.columns:
            level3_by_level2 = {}
            for _, row in kegg_metabolism.iterrows():
                l2 = row.get('Level_2')
                l3 = row.get('Level_3')
                if pd.notna(l2) and pd.notna(l3):
                    if l2 not in level3_by_level2:
                        level3_by_level2[l2] = set()
                    level3_by_level2[l2].add(str(l3))
            
            # Add Level 3 nodes (specific pathways) - limit for clarity
            max_l3_per_l2 = 5  # Maximum pathways per subcategory
            for l2 in level2_unique:
                if l2 in level3_by_level2:
                    pathways = sorted(list(level3_by_level2[l2]))[:max_l3_per_l2]
                    parent_color = node_colors[label_to_idx[l2]]
                    
                    for l3 in pathways:
                        # Clean and truncate pathway names
                        clean_l3 = l3
                        if l3.startswith('0'):
                            parts = l3.split(' ', 1)
                            if len(parts) > 1:
                                clean_l3 = parts[1]
                        
                        # Truncate long names
                        if len(clean_l3) > 40:
                            clean_l3 = clean_l3[:37] + "..."
                        
                        labels.append(clean_l3)
                        
                        # Create lighter version of parent color
                        if parent_color.startswith('#'):
                            r = int(parent_color[1:3], 16)
                            g = int(parent_color[3:5], 16)
                            b = int(parent_color[5:7], 16)
                            # Lighten
                            r = min(255, r + 60)
                            g = min(255, g + 60)
                            b = min(255, b + 60)
                            new_color = f'#{r:02x}{g:02x}{b:02x}'
                        else:
                            new_color = '#E0E0E0'
                        
                        node_colors.append(new_color)
                        label_to_idx[l3] = len(labels) - 1
        
        # Add taxonomy nodes (phyla)
        all_phyla = set()
        mags_subset = []
        
        if self.taxonomy is not None and 'cog' in self.data:
            # Get MAG names from COG data
            mags_subset = list(self.data['cog'].columns)[:40]  # Use more MAGs for metabolism
            
            for mag in mags_subset:
                tax_name = self.convert_mag_names(mag)
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                    all_phyla.add(phylum)
            
            # Generate dynamic colors for phyla
            phyla_list = sorted(list(all_phyla))
            DYNAMIC_PHYLUM_COLORS = self.generate_dynamic_phylum_colors(phyla_list)
            
            # Add phyla nodes
            for phylum in phyla_list:
                labels.append(f"<i><b>{phylum}</b></i>")
                node_colors.append(DYNAMIC_PHYLUM_COLORS[phylum])
                label_to_idx[phylum] = len(labels) - 1
            
            logger.info(f"Added {len(phyla_list)} phyla to metabolism Sankey")
        else:
            phyla_list = []
        
        # Calculate flows
        numeric_cols = kegg_metabolism.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            kegg_metabolism['total'] = kegg_metabolism[numeric_cols].sum(axis=1)
            
            # Flow: Metabolism root ? Level 2 subcategories
            for l2 in level2_unique:
                l2_data = kegg_metabolism[kegg_metabolism['Level_2'] == l2]
                total_flow = l2_data['total'].sum()
                if total_flow > 0:
                    source.append(label_to_idx['METABOLISM'])
                    target.append(label_to_idx[l2])
                    value.append(int(total_flow))
                    link_colors.append(self.hex_to_rgba('#2ECC71', 0.3))
            
            # Flow: Level 2 ? Level 3 (if available)
            if 'Level_3' in kegg_metabolism.columns:
                for _, row in kegg_metabolism.iterrows():
                    l2 = row.get('Level_2')
                    l3 = row.get('Level_3')
                    total = row.get('total', 0)
                    
                    if pd.notna(l2) and pd.notna(l3) and l2 in label_to_idx and l3 in label_to_idx and total > 0:
                        # Check if this link already exists
                        existing = False
                        for i, (s, t) in enumerate(zip(source, target)):
                            if s == label_to_idx[l2] and t == label_to_idx[l3]:
                                value[i] += int(total)
                                existing = True
                                break
                        
                        if not existing:
                            source.append(label_to_idx[l2])
                            target.append(label_to_idx[l3])
                            value.append(int(total))
                            source_color = node_colors[label_to_idx[l2]]
                            link_colors.append(self.hex_to_rgba(source_color, 0.3))
                
                # Flow: Level 3 ? Phyla
                if phyla_list:
                    # Distribute Level 3 flows to phyla based on MAG abundance
                    for l3 in label_to_idx.keys():
                        if l3 in [str(x) for x in level3_by_level2.get(l2, []) for l2 in level2_unique]:
                            # This is a Level 3 pathway
                            # Get total flow into this pathway
                            total_l3_flow = sum(v for s, t, v in zip(source, target, value) if t == label_to_idx.get(l3, -1))
                            
                            if total_l3_flow > 0:
                                # Distribute to phyla (simplified - you could make this more sophisticated)
                                flow_per_phylum = total_l3_flow / len(phyla_list)
                                for phylum in phyla_list:
                                    if phylum in label_to_idx:
                                        source.append(label_to_idx[l3])
                                        target.append(label_to_idx[phylum])
                                        value.append(int(flow_per_phylum * 0.5))  # Scale down for clarity
                                        source_color = node_colors[label_to_idx[l3]]
                                        link_colors.append(self.hex_to_rgba(source_color, 0.25))
            
            else:
                # If no Level 3, connect Level 2 directly to phyla
                if phyla_list:
                    for l2 in level2_unique:
                        if l2 in label_to_idx:
                            l2_data = kegg_metabolism[kegg_metabolism['Level_2'] == l2]
                            total_l2_flow = l2_data['total'].sum()
                            
                            if total_l2_flow > 0:
                                flow_per_phylum = total_l2_flow / len(phyla_list)
                                for phylum in phyla_list:
                                    if phylum in label_to_idx:
                                        source.append(label_to_idx[l2])
                                        target.append(label_to_idx[phylum])
                                        value.append(int(flow_per_phylum * 0.5))
                                        source_color = node_colors[label_to_idx[l2]]
                                        link_colors.append(self.hex_to_rgba(source_color, 0.25))
        
        # Create Sankey
        if len(source) > 0:
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="white", width=1.5),
                    label=labels,
                    color=node_colors,
                    hovertemplate='<b>%{label}</b><br>Total: %{value:,.0f}<extra></extra>',
                    x=[0.1] + [0.3]*len(level2_unique) + [0.6]*(len(labels)-len(level2_unique)-len(phyla_list)-1) + [0.9]*len(phyla_list) if phyla_list else None,
                    y=None  # Let plotly arrange vertically
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_colors,
                    hovertemplate='%{source.label} ? %{target.label}<br>Count: %{value:,.0f}<extra></extra>'
                ),
                arrangement='snap',
                orientation='h'
            )])
            
            # Create title with pathway count
            pathway_count = len([l for l in labels if not l.startswith('<b>') and not l.startswith('<i>')])
            title_text = f"<b>KEGG Metabolism Detailed Flow</b><br>"
            title_text += f"<sup>{len(level2_unique)} subcategories, {pathway_count} pathways"
            if phyla_list:
                title_text += f", {len(phyla_list)} phyla"
            title_text += "</sup>"
            
            fig.update_layout(
                title=dict(
                    text=title_text,
                    font=dict(size=18, family="Arial, sans-serif", color="#27AE60"),
                    x=0.5,
                    xanchor='center'
                ),
                font=dict(size=11, family="Arial, sans-serif"),
                height=1200,
                width=1800,
                paper_bgcolor='#F8F9FA',
                plot_bgcolor='#F8F9FA',
                margin=dict(l=50, r=50, t=100, b=50),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=12,
                    font_family="Arial",
                    bordercolor="#27AE60"
                ),
                annotations=[
                    dict(
                        x=0.05, y=0.5,
                        text="<b>Metabolism</b>",
                        showarrow=False,
                        font=dict(size=12, color="#27AE60"),
                        xref="paper", yref="paper",
                        textangle=-90
                    ),
                    dict(
                        x=0.3, y=1.05,
                        text="<b>Subcategories</b>",
                        showarrow=False,
                        font=dict(size=11, color="#555"),
                        xref="paper", yref="paper"
                    ),
                    dict(
                        x=0.6, y=1.05,
                        text="<b>Pathways</b>",
                        showarrow=False,
                        font=dict(size=11, color="#555"),
                        xref="paper", yref="paper"
                    ),
                    dict(
                        x=0.9, y=1.05,
                        text="<b>Taxonomy</b>",
                        showarrow=False,
                        font=dict(size=11, color="#555"),
                        xref="paper", yref="paper"
                    ) if phyla_list else dict()
                ]
            )
            
            return fig
        
        return None
    
    def _create_metabolism_sankey_without_taxonomy(self, kegg_metabolism, METABOLISM_COLORS):
        """Create metabolism Sankey WITHOUT taxonomy nodes - pathways only"""
        
        # Initialize Sankey components
        labels = []
        node_colors = []
        source = []
        target = []
        value = []
        link_colors = []
        label_to_idx = {}
        
        # Add "Metabolism" as root node
        labels.append("<b>METABOLISM</b>")
        node_colors.append('#2ECC71')  # Bright green for main metabolism
        label_to_idx['METABOLISM'] = 0
        
        # Get unique Level 2 subcategories (metabolism types)
        level2_unique = sorted(kegg_metabolism['Level_2'].dropna().unique())
        
        # Clean up Level 2 names (remove codes if present)
        level2_clean = {}
        for l2 in level2_unique:
            clean_name = l2
            if l2.startswith('09'):
                parts = l2.split(' ', 1)
                if len(parts) > 1:
                    clean_name = parts[1]
            level2_clean[l2] = clean_name
        
        logger.info(f"Found {len(level2_unique)} metabolism subcategories")
        
        # Add Level 2 nodes (metabolism subcategories)
        for l2 in level2_unique:
            clean_name = level2_clean[l2]
            labels.append(f"<b>{clean_name}</b>")
            
            # Get color from predefined or generate
            color = METABOLISM_COLORS.get(clean_name)
            if not color:
                # Check if any key partially matches
                for key, col in METABOLISM_COLORS.items():
                    if key.lower() in clean_name.lower() or clean_name.lower() in key.lower():
                        color = col
                        break
            if not color:
                # Generate a color based on position
                import colorsys
                hue = (len(labels) - 2) / max(len(level2_unique), 1) * 0.8
                rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.9)
                color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            
            node_colors.append(color)
            label_to_idx[l2] = len(labels) - 1
        
        # Get Level 3 pathways if available
        if 'Level_3' in kegg_metabolism.columns:
            level3_by_level2 = {}
            for _, row in kegg_metabolism.iterrows():
                l2 = row.get('Level_2')
                l3 = row.get('Level_3')
                if pd.notna(l2) and pd.notna(l3):
                    if l2 not in level3_by_level2:
                        level3_by_level2[l2] = set()
                    level3_by_level2[l2].add(str(l3))
            
            # Add Level 3 nodes (specific pathways) - increase limit for more detail
            max_l3_per_l2 = 8  # Increased from 5 to show more pathways
            for l2 in level2_unique:
                if l2 in level3_by_level2:
                    pathways = sorted(list(level3_by_level2[l2]))[:max_l3_per_l2]
                    parent_color = node_colors[label_to_idx[l2]]
                    
                    for l3 in pathways:
                        # Clean and truncate pathway names
                        clean_l3 = l3
                        if l3.startswith('0'):
                            parts = l3.split(' ', 1)
                            if len(parts) > 1:
                                clean_l3 = parts[1]
                        
                        # Truncate long names
                        if len(clean_l3) > 45:
                            clean_l3 = clean_l3[:42] + "..."
                        
                        labels.append(clean_l3)
                        
                        # Create lighter version of parent color
                        if parent_color.startswith('#'):
                            r = int(parent_color[1:3], 16)
                            g = int(parent_color[3:5], 16)
                            b = int(parent_color[5:7], 16)
                            # Lighten
                            r = min(255, r + 60)
                            g = min(255, g + 60)
                            b = min(255, b + 60)
                            new_color = f'#{r:02x}{g:02x}{b:02x}'
                        else:
                            new_color = '#E0E0E0'
                        
                        node_colors.append(new_color)
                        label_to_idx[l3] = len(labels) - 1
        
        # Calculate flows WITHOUT taxonomy
        numeric_cols = kegg_metabolism.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            kegg_metabolism['total'] = kegg_metabolism[numeric_cols].sum(axis=1)
            
            # Flow: Metabolism root ? Level 2 subcategories
            for l2 in level2_unique:
                l2_data = kegg_metabolism[kegg_metabolism['Level_2'] == l2]
                total_flow = l2_data['total'].sum()
                if total_flow > 0:
                    source.append(label_to_idx['METABOLISM'])
                    target.append(label_to_idx[l2])
                    value.append(int(total_flow))
                    link_colors.append(self.hex_to_rgba('#2ECC71', 0.3))
            
            # Flow: Level 2 ? Level 3 (if available)
            if 'Level_3' in kegg_metabolism.columns:
                for _, row in kegg_metabolism.iterrows():
                    l2 = row.get('Level_2')
                    l3 = row.get('Level_3')
                    total = row.get('total', 0)
                    
                    if pd.notna(l2) and pd.notna(l3) and l2 in label_to_idx and l3 in label_to_idx and total > 0:
                        # Check if this link already exists
                        existing = False
                        for i, (s, t) in enumerate(zip(source, target)):
                            if s == label_to_idx[l2] and t == label_to_idx[l3]:
                                value[i] += int(total)
                                existing = True
                                break
                        
                        if not existing:
                            source.append(label_to_idx[l2])
                            target.append(label_to_idx[l3])
                            value.append(int(total))
                            source_color = node_colors[label_to_idx[l2]]
                            link_colors.append(self.hex_to_rgba(source_color, 0.3))
        
        # Create Sankey WITHOUT taxonomy nodes
        if len(source) > 0:
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=25,  # Increased padding for better vertical separation
                    thickness=25,
                    line=dict(color="white", width=1.5),
                    label=labels,
                    color=node_colors,
                    hovertemplate='<b>%{label}</b><br>Total: %{value:,.0f}<extra></extra>',
                    # Optimized x positions: closer L1-L2, wider L2-L3
                    x=[0.05, *[0.30]*len(level2_unique), *[0.80]*(len(labels)-len(level2_unique)-1)],
                    y=None  # Let plotly arrange vertically
                ),
                link=dict(
                    source=source,
                    target=target,
                    value=value,
                    color=link_colors,
                    hovertemplate='%{source.label} ? %{target.label}<br>Count: %{value:,.0f}<extra></extra>'
                ),
                arrangement='snap',
                orientation='h'
            )])
            
            # Create title
            pathway_count = len([l for l in labels if not l.startswith('<b>')])
            title_text = f"<b>KEGG Metabolism Pathways Flow</b><br>"
            title_text += f"<sup>{len(level2_unique)} subcategories, {pathway_count} pathways</sup>"
            
            fig.update_layout(
                title=dict(
                    text=title_text,
                    font=dict(size=18, family="Arial, sans-serif", color="#27AE60"),
                    x=0.5,
                    xanchor='center'
                ),
                font=dict(size=11, family="Arial, sans-serif"),
                height=1200,  # Increased height for better vertical spacing
                width=1600,
                paper_bgcolor='#F8F9FA',
                plot_bgcolor='#F8F9FA',
                margin=dict(l=50, r=50, t=100, b=50),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=12,
                    font_family="Arial",
                    bordercolor="#27AE60"
                ),
                annotations=[
                    dict(
                        x=0.05, y=0.5,
                        text="<b>Metabolism</b>",
                        showarrow=False,
                        font=dict(size=12, color="#27AE60"),
                        xref="paper", yref="paper",
                        textangle=-90
                    ),
                    dict(
                        x=0.30, y=1.05,
                        text="<b>Subcategories</b>",
                        showarrow=False,
                        font=dict(size=11, color="#555"),
                        xref="paper", yref="paper"
                    ),
                    dict(
                        x=0.80, y=1.05,
                        text="<b>Pathways</b>",
                        showarrow=False,
                        font=dict(size=11, color="#555"),
                        xref="paper", yref="paper"
                    )
                ]
            )
            
            return fig
        
        return None

    def create_word_optimized_sankey_improved(self, fig, base_name, title):
        """
        Create IMPROVED Word-optimized Sankey with better text spacing
        """
        import plotly.graph_objects as go
        
        # Get the original data
        original_data = fig.data[0]
        
        # Safely get x and y coordinates if they exist
        x_coords = None
        y_coords = None
        if hasattr(original_data.node, 'x'):
            x_coords = original_data.node.x
        if hasattr(original_data.node, 'y'):
            y_coords = original_data.node.y
        
        # Create COMPACT figure with BETTER SPACING
        fig_word = go.Figure(data=[go.Sankey(
            node=dict(
                pad=50,  # Increased from 40 for better separation
                thickness=60,  # Slightly reduced from 70
                line=dict(color="black", width=3),  # Reduced from 4
                label=original_data.node.label,
                color=original_data.node.color,
                hovertemplate=original_data.node.hovertemplate,
                x=x_coords,
                y=y_coords
            ),
            link=dict(
                source=original_data.link.source,
                target=original_data.link.target,
                value=original_data.link.value,
                color=original_data.link.color,
                hovertemplate=original_data.link.hovertemplate,
                line=dict(width=1, color="rgba(0,0,0,0.1)")
            ),
            orientation='h',
            valueformat=".0f",
            valuesuffix="",
            arrangement='snap'
        )])
        
        # IMPROVED LAYOUT with better font sizing
        fig_word.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=48,  # Reduced from 56
                    family="Arial Black, sans-serif",
                    color="#000000"
                ),
                x=0.5,
                xanchor='center',
                y=0.98,
                yanchor='top'
            ),
            font=dict(
                size=36,  # Reduced from 42
                family="Arial, sans-serif",
                color="#000000"
            ),
            # Better proportions
            height=1800,
            width=2800,  # Increased width for better spacing
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=180, r=180, t=200, b=150),  # Increased margins
            showlegend=False
        )
        
        # Update text font for better readability
        fig_word.update_traces(
            textfont=dict(
                size=34,  # Reduced from 40 for better fit
                family="Arial, sans-serif",
                color="black"
            ),
            selector=dict(type='sankey')
        )
        
        return fig_word
    
    def save_publication_figure_enhanced_improved(self, fig, base_name, title):
        """
        Enhanced save with IMPROVED spacing for Word optimization
        """
        output_base = self.output_dir / base_name
        
        # Save standard HTML
        fig.write_html(f"{output_base}.html")
        logger.info(f"? HTML saved: {base_name}.html")
        
        # Create IMPROVED COMPACT version
        fig_compact = self.create_word_optimized_sankey_improved(fig, base_name, title)
        
        # Create ULTRA version with MAXIMUM spacing optimization
        # Calculate x positions for ULTRA version
        node_count = len(fig.data[0].node.label)
        x_coords_ultra = []
        
        if hasattr(fig.data[0].node, 'label'):
            labels = fig.data[0].node.label
            for i, label in enumerate(labels):
                if i == 0:  # METABOLISM (root)
                    x_coords_ultra.append(0.05)
                elif '<b>METABOLISM</b>' in str(label):  # Root node
                    x_coords_ultra.append(0.05)
                elif label.startswith('<b>') and not label.startswith('<i>'):  # Level 2
                    x_coords_ultra.append(0.20)  # Very close to root
                elif label.startswith('<i>'):  # Taxonomy nodes (if present)
                    x_coords_ultra.append(0.90)
                else:  # Level 3 (pathways)
                    x_coords_ultra.append(0.70)  # Far from Level 2 for maximum separation
        else:
            x_coords_ultra = None
        
        fig_ultra = go.Figure(data=[go.Sankey(
            node=dict(
                pad=70,  # Maximum vertical padding for clear separation
                thickness=45,  # Thinner bars to maximize space for text
                line=dict(color="black", width=2),
                label=fig.data[0].node.label,
                color=fig.data[0].node.color,
                x=x_coords_ultra,  # Use calculated positions
                y=None  # Automatic vertical arrangement
            ),
            link=dict(
                source=fig.data[0].link.source,
                target=fig.data[0].link.target,
                value=fig.data[0].link.value,
                color=fig.data[0].link.color,
                line=dict(width=0.5, color="rgba(0,0,0,0.05)")
            ),
            orientation='h',
            arrangement='snap'
        )])
        
        # ULTRA layout with maximum readability
        fig_ultra.update_layout(
            title=dict(
                text=f"<b>{title}</b>",
                font=dict(
                    size=48,  # Clear title
                    family="Arial Black",
                    color="#000000"
                ),
                x=0.5,
                xanchor='center',
                y=0.98
            ),
            font=dict(
                size=34,  # Balanced font size
                family="Arial, sans-serif",
                color="#000000"
            ),
            # Maximum height for best vertical spacing
            height=2400,  # Increased significantly
            width=3400,  # Wider for better horizontal distribution
            paper_bgcolor='white',
            plot_bgcolor='white',
            margin=dict(l=180, r=180, t=160, b=100),
            showlegend=False
        )
        
        fig_ultra.update_traces(
            textfont=dict(
                size=32,  # Clear, readable text
                family="Arial, sans-serif",
                color="black"
            ),
            selector=dict(type='sankey')
        )
        
        # Try to export
        try:
            # Compact Word version
            fig_compact.write_image(
                f"{output_base}_WORD_COMPACT.png", 
                width=2800, 
                height=2200, 
                scale=2
            )
            logger.info(f"? COMPACT PNG saved: {base_name}_WORD_COMPACT.png")
            
            # Ultra version
            fig_ultra.write_image(
                f"{output_base}_WORD_ULTRA.png",
                width=3400,
                height=2400,
                scale=2
            )
            logger.info(f"? ULTRA PNG saved: {base_name}_WORD_ULTRA.png")
            
        except Exception as e:
            logger.warning(f"Static export failed: {e}")
        
        # Save HTML versions for manual export
        
        # COMPACT HTML
        config_compact = {
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f"{base_name}_COMPACT",
                'height': 2200,
                'width': 2800,
                'scale': 2
            },
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToAdd': ['downloadImage']
        }
        
        fig_compact.write_html(
            f"{output_base}_COMPACT.html",
            config=config_compact,
            include_plotlyjs='cdn'
        )
        
        # ULTRA HTML
        config_ultra = {
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f"{base_name}_ULTRA",
                'height': 2400,
                'width': 3400,
                'scale': 2
            },
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToAdd': ['downloadImage']
        }
        
        fig_ultra.write_html(
            f"{output_base}_ULTRA.html",
            config=config_ultra,
            include_plotlyjs='cdn'
        )
        
        logger.info(f"\n? COMPACT version: {base_name}_COMPACT.html")
        logger.info(f"   ? Size: 2800x2200 (taller), Font: 30pt labels")
        logger.info(f"   ? Level spacing: 0.05 ? 0.25 ? 0.65 (closer L1-L2, wider L2-L3)")
        logger.info(f"? ULTRA version: {base_name}_ULTRA.html")
        logger.info(f"   ? Size: 3400x2400 (maximum height), Font: 32pt labels")
        logger.info(f"   ? Level spacing: 0.05 ? 0.20 ? 0.70 (minimal L1-L2, maximum L2-L3)")
        logger.info(f"   ? Optimized vertical padding (70px) for clear text separation!")



    def create_functional_pca_analysis(self):
        """
        Create high-quality PCA plots for functional analysis of MAGs
        Analyzes COG, KO, and CAZyme profiles with proper taxonomy mapping
        """
        logger.info("Creating publication-quality functional PCA analysis...")
        
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        import matplotlib.patches as mpatches
        from matplotlib.patches import Ellipse  # Added this import
        import seaborn as sns
        
        # Set publication quality parameters
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        
        # Check available data
        available_data = {}
        if 'cog' in self.data:
            available_data['COG'] = self.data['cog']
        if 'ko' in self.data:
            available_data['KEGG'] = self.data['ko']
        if 'cazyme' in self.data:
            available_data['CAZyme'] = self.data['cazyme']
        
        if not available_data:
            logger.warning("No functional data available for PCA analysis")
            return
        
        # Process each functional dataset
        for data_type, data_df in available_data.items():
            logger.info(f"Processing PCA for {data_type} data...")
            self._create_pca_for_dataset(data_type, data_df)
        
        # Create combined functional PCA if multiple datasets available
        if len(available_data) > 1:
            self._create_combined_functional_pca(available_data)
        
        logger.info("Functional PCA analysis completed!")
    
    def _create_pca_for_dataset(self, data_type, data_df):
        """Create publication-quality PCA analysis for a specific functional dataset"""
        
        from matplotlib.patches import Ellipse  # Import here too for safety
        
        # Prepare data - handle different data structures
        if data_type == 'CAZyme':
            # CAZyme data: MAGs are index (rows), CAZyme categories are columns
            # Need to transpose so MAGs are rows for PCA
            data_matrix = data_df
            # Remove 'Total' column if exists
            if 'Total' in data_matrix.columns:
                data_matrix = data_matrix.drop('Total', axis=1)
            mag_names = data_matrix.index.tolist()
        else:
            # COG and KO data: MAGs are columns, features are rows
            # Transpose so MAGs are rows
            data_matrix = data_df.T
            mag_names = data_matrix.index.tolist()
        
        logger.info(f"Data matrix shape for {data_type}: {data_matrix.shape}")
        logger.info(f"Sample MAG names: {mag_names[:3]}")
        
        # Map MAGs to taxonomy with proper name conversion
        mag_taxonomy = {}
        mag_phyla = []
        unmatched_count = 0
        
        for mag in mag_names:
            # Use the convert_mag_names function to handle different naming formats
            tax_name = self.convert_mag_names(mag)
            
            # Debug logging for first few MAGs
            if unmatched_count < 3:
                logger.info(f"Converting {data_type} MAG: {mag} -> {tax_name}")
            
            if self.taxonomy is not None:
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                    genus = tax_row.iloc[0].get('genus', 'Unknown')
                    if unmatched_count < 3:
                        logger.info(f"  ? Matched to phylum: {phylum}")
                else:
                    phylum = genus = 'Unknown'
                    unmatched_count += 1
                    if unmatched_count <= 3:
                        logger.warning(f"  ? No match found for: {tax_name}")
            else:
                phylum = genus = 'Unknown'
            
            mag_taxonomy[mag] = {'phylum': phylum, 'genus': genus}
            mag_phyla.append(phylum)
        
        logger.info(f"Taxonomy matching: {len(mag_names) - unmatched_count}/{len(mag_names)} MAGs matched")
        
        # Get unique phyla and generate colors
        unique_phyla = sorted(list(set(mag_phyla)))
        phylum_colors = self._generate_publication_colors(unique_phyla)
        
        # Remove rows/columns with all zeros or low variance
        data_matrix = data_matrix.loc[:, (data_matrix != 0).any(axis=0)]
        data_matrix = data_matrix.loc[(data_matrix != 0).any(axis=1), :]
        
        # Standardize the data
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data_matrix)
        
        # Perform PCA
        pca = PCA()
        data_pca = pca.fit_transform(data_scaled)
        
        # Calculate explained variance
        explained_var = pca.explained_variance_ratio_
        cumulative_var = np.cumsum(explained_var)
        
        # ========== CREATE PUBLICATION QUALITY FIGURE ==========
        fig = plt.figure(figsize=(16, 10))
        
        # Create grid
        gs = plt.GridSpec(2, 3, figure=fig, 
                         width_ratios=[1.2, 1.2, 1],
                         height_ratios=[1, 1],
                         hspace=0.3, wspace=0.3)
        
        # ========== MAIN PCA PLOT (PC1 vs PC2) ==========
        ax1 = fig.add_subplot(gs[0, :2])
        
        # Plot points by phylum with nice colors
        for phylum in unique_phyla:
            mask = [p == phylum for p in mag_phyla]
            indices = [i for i, m in enumerate(mask) if m]
            if indices:
                ax1.scatter(data_pca[indices, 0], data_pca[indices, 1],
                          c=phylum_colors[phylum], label=phylum,
                          s=120, alpha=0.8, edgecolors='white', linewidth=1.5)
        
        # Add confidence ellipses for each phylum
        for phylum in unique_phyla:
            mask = [p == phylum for p in mag_phyla]
            indices = [i for i, m in enumerate(mask) if m]
            if len(indices) > 2:  # Need at least 3 points for ellipse
                points = data_pca[indices, :2]
                mean = np.mean(points, axis=0)
                cov = np.cov(points.T)
                try:
                    eigenvalues, eigenvectors = np.linalg.eig(cov)
                    if np.all(eigenvalues > 0):  # Check for positive eigenvalues
                        angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
                        width, height = 2 * np.sqrt(eigenvalues) * 1.5  # 1.5 sigma
                        ellipse = Ellipse(mean, width, height, angle=angle,
                                         facecolor=phylum_colors[phylum], alpha=0.1,
                                         edgecolor=phylum_colors[phylum], linewidth=1)
                        ax1.add_patch(ellipse)
                except:
                    # Skip ellipse if covariance calculation fails
                    pass
        
        ax1.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}% variance)', fontsize=12, fontweight='bold')
        ax1.set_ylabel(f'PC2 ({explained_var[1]*100:.1f}% variance)', fontsize=12, fontweight='bold')
        ax1.set_title(f'{data_type} Functional Profile - Principal Component Analysis', 
                     fontsize=14, fontweight='bold')
        
        # Remove grid lines but keep axes
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        ax1.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        
        # ========== PC2 vs PC3 ==========
        if data_pca.shape[1] >= 3:
            ax2 = fig.add_subplot(gs[1, 0])
            
            for phylum in unique_phyla:
                mask = [p == phylum for p in mag_phyla]
                indices = [i for i, m in enumerate(mask) if m]
                if indices:
                    ax2.scatter(data_pca[indices, 1], data_pca[indices, 2],
                              c=phylum_colors[phylum],
                              s=100, alpha=0.8, edgecolors='white', linewidth=1)
            
            ax2.set_xlabel(f'PC2 ({explained_var[1]*100:.1f}% variance)', fontsize=11, fontweight='bold')
            ax2.set_ylabel(f'PC3 ({explained_var[2]*100:.1f}% variance)', fontsize=11, fontweight='bold')
            ax2.set_title(f'PC2 vs PC3', fontsize=12, fontweight='bold')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
            ax2.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        
        # ========== PC1 vs PC3 ==========
        if data_pca.shape[1] >= 3:
            ax3 = fig.add_subplot(gs[1, 1])
            
            for phylum in unique_phyla:
                mask = [p == phylum for p in mag_phyla]
                indices = [i for i, m in enumerate(mask) if m]
                if indices:
                    ax3.scatter(data_pca[indices, 0], data_pca[indices, 2],
                              c=phylum_colors[phylum],
                              s=100, alpha=0.8, edgecolors='white', linewidth=1)
            
            ax3.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}% variance)', fontsize=11, fontweight='bold')
            ax3.set_ylabel(f'PC3 ({explained_var[2]*100:.1f}% variance)', fontsize=11, fontweight='bold')
            ax3.set_title(f'PC1 vs PC3', fontsize=12, fontweight='bold')
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            ax3.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
            ax3.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        
        # ========== SCREE PLOT ==========
        ax4 = fig.add_subplot(gs[0, 2])
        
        n_components_to_show = min(15, len(explained_var))
        x_pos = np.arange(1, n_components_to_show + 1)
        
        # Bar plot for individual variance
        bars = ax4.bar(x_pos, explained_var[:n_components_to_show] * 100,
                      alpha=0.7, color='#2E86AB', edgecolor='#1C5577', linewidth=1.5)
        
        # Cumulative line
        ax4_twin = ax4.twinx()
        line = ax4_twin.plot(x_pos, cumulative_var[:n_components_to_show] * 100,
                           color='#A23B72', linewidth=2.5, marker='o', markersize=6,
                           markerfacecolor='white', markeredgecolor='#A23B72', markeredgewidth=2)
        
        # 90% threshold
        ax4_twin.axhline(y=90, color='#C73E1D', linestyle='--', linewidth=2, alpha=0.7)
        
        ax4.set_xlabel('Principal Component', fontsize=11, fontweight='bold')
        ax4.set_ylabel('Variance Explained (%)', fontsize=11, fontweight='bold', color='#2E86AB')
        ax4_twin.set_ylabel('Cumulative Variance (%)', fontsize=11, fontweight='bold', color='#A23B72')
        ax4.set_title('Scree Plot', fontsize=12, fontweight='bold')
        ax4.tick_params(axis='y', labelcolor='#2E86AB')
        ax4_twin.tick_params(axis='y', labelcolor='#A23B72')
        ax4.spines['top'].set_visible(False)
        ax4.set_xticks(x_pos[::2])
        
        # ========== TOP CONTRIBUTING FEATURES ==========
        ax5 = fig.add_subplot(gs[1, 2])
        
        # Get absolute loadings for PC1
        pc1_loadings = np.abs(pca.components_[0])
        top_indices = np.argsort(pc1_loadings)[-12:]
        
        feature_names = data_matrix.columns.tolist()
        top_features = [feature_names[i] for i in top_indices]
        top_values = pc1_loadings[top_indices]
        
        # Clean feature names
        clean_features = []
        for feat in top_features:
            if len(feat) > 20:
                clean_features.append(feat[:17] + '...')
            else:
                clean_features.append(feat)
        
        y_pos = np.arange(len(clean_features))
        bars = ax5.barh(y_pos, top_values, color='#4A7C59', alpha=0.8, 
                       edgecolor='#2F5233', linewidth=1.5)
        
        # Color gradient for bars
        for i, bar in enumerate(bars):
            bar.set_facecolor(plt.cm.Greens(0.4 + i * 0.05))
        
        ax5.set_yticks(y_pos)
        ax5.set_yticklabels(clean_features, fontsize=9)
        ax5.set_xlabel('Absolute Loading', fontsize=11, fontweight='bold')
        ax5.set_title('Top Contributing Features (PC1)', fontsize=12, fontweight='bold')
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        
        # ========== LEGEND ==========
        # Create custom legend with counts
        legend_elements = []
        for phylum in unique_phyla:
            n_mags = sum(1 for p in mag_phyla if p == phylum)
            if n_mags > 0:
                legend_elements.append(plt.scatter([], [], s=120, c=phylum_colors[phylum],
                                                  alpha=0.8, edgecolors='white', linewidth=1.5,
                                                  label=f'{phylum} (n={n_mags})'))
        
        ax1.legend(handles=legend_elements, loc='best', frameon=True, 
                  fancybox=True, shadow=True, fontsize=10,
                  title='Phylum', title_fontsize=11, 
                  edgecolor='gray', facecolor='white', framealpha=0.95)
        
        # Add subtitle with statistics
        fig.text(0.5, 0.96, 
                f'{len(mag_names)} MAGs | {data_matrix.shape[1]} Features | {len(unique_phyla)} Phyla',
                ha='center', fontsize=11, style='italic', color='#333333')
        
        plt.suptitle(f'{data_type} Functional Analysis - PCA', 
                    fontsize=16, fontweight='bold', y=0.99)
        
        # Save figure
        output_name = f'pca_{data_type.lower()}_publication'
        plt.savefig(self.output_dir / f'{output_name}.pdf', dpi=600, bbox_inches='tight', facecolor='white')
        plt.savefig(self.output_dir / f'{output_name}.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(self.output_dir / f'{output_name}.svg', format='svg', bbox_inches='tight', facecolor='white')
        plt.close()
        
        logger.info(f"Publication-quality PCA for {data_type} saved successfully!")   

    def _generate_publication_colors(self, phyla_list):
        """Generate publication-quality colors for phyla - completely dynamic"""
        import colorsys
        
        phylum_colors = {}
        
        # Special case for Unknown
        if 'Unknown' in phyla_list:
            phylum_colors['Unknown'] = '#CCCCCC'  # Light gray
            remaining_phyla = [p for p in phyla_list if p != 'Unknown']
        else:
            remaining_phyla = phyla_list.copy()
        
        # Generate distinct colors for all phyla dynamically
        n_phyla = len(remaining_phyla)
        
        if n_phyla > 0:
            # Use different strategies based on number of phyla
            if n_phyla <= 12:
                # For small number, use maximally distinct colors
                for i, phylum in enumerate(remaining_phyla):
                    hue = i / n_phyla
                    saturation = 0.85
                    value = 0.85
                    rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                    phylum_colors[phylum] = '#{:02x}{:02x}{:02x}'.format(
                        int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            
            elif n_phyla <= 20:
                # For medium number, use golden ratio for better distribution
                golden_ratio = 0.618033988749895
                for i, phylum in enumerate(remaining_phyla):
                    hue = (i * golden_ratio) % 1
                    # Vary saturation and value for more distinction
                    saturation = 0.65 + (i % 3) * 0.15  # 0.65 to 0.95
                    value = 0.75 + (i % 2) * 0.15  # 0.75 to 0.90
                    rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                    phylum_colors[phylum] = '#{:02x}{:02x}{:02x}'.format(
                        int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            
            else:
                # For large number, use multiple strategies
                # Use ColorBrewer-inspired approach
                base_hues = [0, 0.15, 0.3, 0.45, 0.6, 0.75]  # Red, Yellow, Green, Cyan, Blue, Purple
                
                for i, phylum in enumerate(remaining_phyla):
                    base_idx = i % len(base_hues)
                    variation = (i // len(base_hues)) / (n_phyla // len(base_hues) + 1)
                    
                    hue = (base_hues[base_idx] + variation * 0.1) % 1
                    saturation = 0.5 + (1 - variation) * 0.4  # Decrease saturation with variation
                    value = 0.6 + variation * 0.3  # Increase value with variation
                    
                    rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                    phylum_colors[phylum] = '#{:02x}{:02x}{:02x}'.format(
                        int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
        
        return phylum_colors

    def _create_combined_functional_pca(self, available_data):
        """Create combined PCA for all functional data types"""
        logger.info("Creating combined functional PCA...")
        
        from matplotlib.patches import Ellipse
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        
        # Prepare combined data matrix
        combined_features = []
        feature_labels = []
        
        # Get common MAGs across all datasets
        common_mags = None
        for data_type, data_df in available_data.items():
            if data_type == 'CAZyme':
                current_mags = set(data_df.index)
            else:
                current_mags = set(data_df.columns)
            
            if common_mags is None:
                common_mags = current_mags
            else:
                common_mags = common_mags & current_mags
        
        common_mags = list(common_mags)
        logger.info(f"Found {len(common_mags)} common MAGs across all datasets")
        
        if len(common_mags) == 0:
            logger.warning("No common MAGs found across datasets")
            return
        
        # Build combined matrix
        for data_type, data_df in available_data.items():
            if data_type == 'CAZyme':
                data_matrix = data_df.loc[common_mags, :]
                if 'Total' in data_matrix.columns:
                    data_matrix = data_matrix.drop('Total', axis=1)
            else:
                data_matrix = data_df[common_mags].T
            
            # Add prefix to feature names
            for col in data_matrix.columns:
                feature_labels.append(f'{data_type}_{col}')
            
            combined_features.append(data_matrix)
        
        # Combine all features
        combined_df = pd.concat(combined_features, axis=1)
        combined_df.columns = feature_labels
        
        # Get taxonomy for common MAGs
        mag_phyla = []
        for mag in common_mags:
            tax_name = self.convert_mag_names(mag)
            if self.taxonomy is not None:
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                else:
                    phylum = 'Unknown'
            else:
                phylum = 'Unknown'
            mag_phyla.append(phylum)
        
        # Get unique phyla and colors
        unique_phyla = sorted(list(set(mag_phyla)))
        phylum_colors = self._generate_publication_colors(unique_phyla)
        
        # Standardize and perform PCA
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(combined_df)
        
        pca = PCA()
        data_pca = pca.fit_transform(data_scaled)
        
        # Create visualization
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # ========== SUBPLOT 1: PCA PLOT ==========
        ax1 = axes[0]
        
        # Plot by phylum
        for phylum in unique_phyla:
            mask = [p == phylum for p in mag_phyla]
            indices = [i for i, m in enumerate(mask) if m]
            if indices:
                ax1.scatter(data_pca[indices, 0], data_pca[indices, 1],
                          c=phylum_colors[phylum], label=phylum,
                          s=150, alpha=0.8, edgecolors='white', linewidth=2)
        
        # Add confidence ellipses
        for phylum in unique_phyla:
            mask = [p == phylum for p in mag_phyla]
            indices = [i for i, m in enumerate(mask) if m]
            if len(indices) > 2:
                points = data_pca[indices, :2]
                mean = np.mean(points, axis=0)
                cov = np.cov(points.T)
                try:
                    eigenvalues, eigenvectors = np.linalg.eig(cov)
                    if np.all(eigenvalues > 0):
                        angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
                        width, height = 2 * np.sqrt(eigenvalues) * 1.5
                        ellipse = Ellipse(mean, width, height, angle=angle,
                                         facecolor=phylum_colors[phylum], alpha=0.1,
                                         edgecolor=phylum_colors[phylum], linewidth=1)
                        ax1.add_patch(ellipse)
                except:
                    pass
        
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', 
                     fontsize=12, fontweight='bold')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', 
                     fontsize=12, fontweight='bold')
        ax1.set_title('Combined Functional Profile (All Data Types)', 
                    fontsize=14, fontweight='bold')
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        ax1.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
        
        # Legend
        ax1.legend(loc='best', frameon=True, fancybox=True, shadow=True,
                 fontsize=10, title='Phylum', title_fontsize=11)
        
        # ========== SUBPLOT 2: Data Type Contribution ==========
        ax2 = axes[1]
        
        # Calculate contribution of each data type to PC1
        pc1_loadings = np.abs(pca.components_[0])
        data_type_contributions = {}
        
        for data_type in available_data.keys():
            mask = [data_type in label for label in feature_labels]
            contribution = np.sum(pc1_loadings[mask])
            data_type_contributions[data_type] = contribution
        
        # Normalize contributions
        total_contrib = sum(data_type_contributions.values())
        if total_contrib > 0:
            for key in data_type_contributions:
                data_type_contributions[key] = (data_type_contributions[key] / total_contrib) * 100
        
        # Create bar plot
        data_types = list(data_type_contributions.keys())
        contributions = list(data_type_contributions.values())
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(data_types)))
        bars = ax2.bar(data_types, contributions, color=colors, edgecolor='black', linewidth=1.5)
        
        ax2.set_ylabel('Contribution to PC1 (%)', fontsize=12, fontweight='bold')
        ax2.set_title('Data Type Contributions', fontsize=14, fontweight='bold')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        
        # Add value labels on bars
        for bar, val in zip(bars, contributions):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=10)
        
        plt.suptitle('Combined Functional Analysis - PCA', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save
        plt.savefig(self.output_dir / 'pca_combined_publication.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(self.output_dir / 'pca_combined_publication.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / 'pca_combined_publication.svg', format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info("Combined functional PCA saved successfully!")    

    def create_enhanced_network_analysis(self):
        """Create enhanced network with MAG labels and better layout"""
        logger.info("Creating enhanced MAG network analysis...")
        
        if not NETWORKX_AVAILABLE:
            logger.warning("NetworkX not available, skipping network analysis")
            return
        
        if 'ko' not in self.data:
            logger.warning("No KO data for network analysis")
            return
        
        ko_data = self.data['ko']
        mags = list(ko_data.columns)[:40]  # Limit for clarity
        
        # Calculate similarity matrix
        ko_subset = ko_data[mags]
        ko_binary = (ko_subset > 0).astype(int)
        
        # Jaccard similarity
        similarity_matrix = np.zeros((len(mags), len(mags)))
        for i, mag1 in enumerate(mags):
            for j, mag2 in enumerate(mags):
                if i != j:
                    intersection = (ko_binary[mag1] & ko_binary[mag2]).sum()
                    union = (ko_binary[mag1] | ko_binary[mag2]).sum()
                    if union > 0:
                        similarity_matrix[i, j] = intersection / union
        
        # Create network
        G = nx.Graph()
        
        # Add nodes with attributes
        mag_info = {}
        for mag in mags:
            tax_name = self.convert_mag_names(mag)
            
            # Get taxonomy
            if self.taxonomy is not None:
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                    genus = tax_row.iloc[0].get('genus', 'Unknown')
                    species = tax_row.iloc[0].get('species', '')
                else:
                    phylum = genus = 'Unknown'
                    species = ''
            else:
                phylum = genus = 'Unknown'
                species = ''
            
            # Create short label
            if '_' in mag:
                parts = mag.split('_')
                sample = parts[0] if parts[0].startswith('SRR') else parts[0][:5]
                bin_info = '_'.join(parts[-2:]) if len(parts) > 2 else parts[-1]
                short_label = f"{sample}_{bin_info[:10]}"
            else:
                short_label = mag[:15]
            
            mag_info[mag] = {
                'phylum': phylum,
                'genus': genus,
                'species': species,
                'label': short_label,
                'ko_count': ko_binary[mag].sum()
            }
            
            G.add_node(mag, **mag_info[mag])
        
        # Add edges
        threshold = np.percentile(similarity_matrix[similarity_matrix > 0], 75)
        for i, mag1 in enumerate(mags):
            for j, mag2 in enumerate(mags):
                if i < j and similarity_matrix[i, j] > threshold:
                    G.add_edge(mag1, mag2, weight=similarity_matrix[i, j])
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
        
        # Subplot 1: Network colored by phylum
        pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
        
        # Draw edges
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        nx.draw_networkx_edges(G, pos, alpha=0.2, width=[w*2 for w in weights], 
                              edge_color='gray', ax=ax1)
        
        # Draw nodes colored by phylum
        phyla = list(set(nx.get_node_attributes(G, 'phylum').values()))
        phylum_colors = {p: PHYLUM_COLORS.get(p, '#B09C85') for p in phyla}
        
        node_colors = [phylum_colors[G.nodes[node]['phylum']] for node in G.nodes()]
        node_sizes = [G.nodes[node]['ko_count'] * 10 for node in G.nodes()]
        
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                              alpha=0.7, edgecolors='black', linewidths=0.5, ax=ax1)
        
        # Add labels for all nodes
        labels = {node: G.nodes[node]['label'] for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax1)
        
        # Legend for phylum
        legend_elements = [plt.scatter([], [], c=color, s=100, label=phylum, alpha=0.7)
                          for phylum, color in list(phylum_colors.items())[:8]]
        ax1.legend(handles=legend_elements, title='Phylum', loc='upper left', 
                  fontsize=7, title_fontsize=8)
        
        ax1.set_title('MAG Functional Similarity Network (Colored by Phylum)', 
                     fontsize=11, fontweight='bold')
        ax1.axis('off')
        
        # Subplot 2: Network colored by degree centrality
        degree_centrality = nx.degree_centrality(G)
        
        # Color by centrality
        centrality_colors = [degree_centrality[node] for node in G.nodes()]
        
        nx.draw_networkx_edges(G, pos, alpha=0.2, width=[w*2 for w in weights], 
                              edge_color='gray', ax=ax2)
        
        nodes = nx.draw_networkx_nodes(G, pos, node_color=centrality_colors, 
                                       node_size=node_sizes, cmap='YlOrRd',
                                       alpha=0.7, edgecolors='black', 
                                       linewidths=0.5, ax=ax2, vmin=0, vmax=max(centrality_colors))
        
        # Add labels for highly connected nodes
        high_centrality_nodes = [node for node, cent in degree_centrality.items() if cent > 0.3]
        high_labels = {node: G.nodes[node]['label'] for node in high_centrality_nodes}
        nx.draw_networkx_labels(G, pos, high_labels, font_size=7, font_weight='bold', ax=ax2)
        
        # Add labels for other nodes in smaller font
        other_nodes = [node for node in G.nodes() if node not in high_centrality_nodes]
        other_labels = {node: G.nodes[node]['label'] for node in other_nodes}
        nx.draw_networkx_labels(G, pos, other_labels, font_size=5, alpha=0.7, ax=ax2)
        
        ax2.set_title('MAG Network (Colored by Centrality)', fontsize=11, fontweight='bold')
        ax2.axis('off')
        
        # Add colorbar for centrality
        sm = plt.cm.ScalarMappable(cmap='YlOrRd', 
                                   norm=plt.Normalize(vmin=0, vmax=max(centrality_colors)))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax2, fraction=0.046, pad=0.04)
        cbar.set_label('Degree Centrality', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        
        plt.suptitle('MAG Functional Network Analysis', fontsize=13, fontweight='bold')
        plt.tight_layout()
        
        plt.savefig(self.output_dir / 'enhanced_mag_network.pdf', dpi=600)
        plt.savefig(self.output_dir / 'enhanced_mag_network.png', dpi=300)
        plt.close()
        
        logger.info("Enhanced network analysis saved")

    def create_functional_network_visualization(self):
        """
        Create network visualization with pie charts showing functional composition at each node
        Similar to phylogenetic networks with functional annotation distributions
        """
        logger.info("Creating functional network visualization with pie charts...")
        
        if not NETWORKX_AVAILABLE:
            logger.warning("NetworkX not available, skipping network visualization")
            return
        
        import networkx as nx
        from matplotlib.patches import Circle, Wedge
        from matplotlib.collections import PatchCollection
        import matplotlib.patches as mpatches
        
        # Check available data
        available_data = {}
        if 'cog' in self.data:
            available_data['COG'] = self.data['cog']
        if 'cazyme' in self.data:
            available_data['CAZyme'] = self.data['cazyme']
        if 'ko' in self.data:
            available_data['KEGG'] = self.data['ko']
        
        if not available_data:
            logger.warning("No functional data available for network visualization")
            return
        
        # Create visualizations for each data type
        for data_type, data_df in available_data.items():
            self._create_network_with_pies(data_type, data_df)
        
        # Create combined multi-layer network if multiple datasets
        if len(available_data) > 1:
            self._create_multilayer_functional_network(available_data)
        
        logger.info("Functional network visualizations completed!")
    
    def _create_network_with_pies(self, data_type, data_df):
        """Create network with pie charts for specific functional data"""
        
        import networkx as nx
        from matplotlib.patches import Wedge, Circle
        import matplotlib.patches as mpatches
        
        # Prepare data based on type
        if data_type == 'CAZyme':
            # CAZyme data structure
            mag_data = data_df
            if 'Total' in mag_data.columns:
                mag_data = mag_data.drop('Total', axis=1)
            feature_categories = {
                'GH': '#FF6B6B',  # Glycoside Hydrolases - Red
                'GT': '#4ECDC4',  # Glycosyl Transferases - Teal
                'PL': '#45B7D1',  # Polysaccharide Lyases - Blue
                'CE': '#96CEB4',  # Carbohydrate Esterases - Green
                'AA': '#FFEAA7',  # Auxiliary Activities - Yellow
                'CBM': '#DDA0DD'  # Carbohydrate-Binding Modules - Plum
            }
        elif data_type == 'COG':
            # COG data - aggregate by main categories
            mag_data = data_df.T
            # Aggregate COG categories
            cog_aggregated = {}
            for main_cat, cog_list in COG_MAIN_CATEGORIES.items():
                if any(cog in data_df.index for cog in cog_list):
                    cog_aggregated[main_cat] = data_df.loc[
                        data_df.index.isin(cog_list)
                    ].sum(axis=0)
            
            mag_data = pd.DataFrame(cog_aggregated).T
            feature_categories = {
                'INFORMATION STORAGE AND PROCESSING': '#3498DB',
                'CELLULAR PROCESSES AND SIGNALING': '#E74C3C',
                'METABOLISM': '#2ECC71',
                'POORLY CHARACTERIZED': '#9B59B6'
            }
        else:  # KEGG
            # For KEGG, use top pathways or aggregate by level
            mag_data = data_df.T
            # Get top 6 most abundant features
            feature_sums = data_df.sum(axis=1)
            top_features = feature_sums.nlargest(6).index
            mag_data = mag_data[top_features]
            
            # Generate colors for KEGG features
            feature_categories = {}
            colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
            for i, feat in enumerate(top_features):
                feature_categories[feat] = colors[i % len(colors)]
        
        # Get MAG names and taxonomy
        mag_names = mag_data.index.tolist()
        mag_taxonomy = {}
        mag_phyla = []
        
        for mag in mag_names:
            tax_name = self.convert_mag_names(mag)
            if self.taxonomy is not None:
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                    genus = tax_row.iloc[0].get('genus', 'Unknown')
                else:
                    phylum = genus = 'Unknown'
            else:
                phylum = genus = 'Unknown'
            
            mag_taxonomy[mag] = {'phylum': phylum, 'genus': genus}
            mag_phyla.append(phylum)
        
        # Get unique phyla and colors
        unique_phyla = sorted(list(set(mag_phyla)))
        phylum_colors = self._generate_publication_colors(unique_phyla)
        
        # Calculate similarity matrix for network edges
        from scipy.spatial.distance import pdist, squareform
        from scipy.stats import pearsonr
        
        # Normalize data for similarity calculation
        mag_data_norm = mag_data.div(mag_data.sum(axis=1), axis=0).fillna(0)
        
        # Calculate correlation matrix
        similarity_matrix = mag_data_norm.T.corr(method='pearson')
        
        # Create network
        G = nx.Graph()
        
        # Add nodes
        for i, mag in enumerate(mag_names):
            G.add_node(mag, phylum=mag_phyla[i], genus=mag_taxonomy[mag]['genus'])
        
        # Add edges (only strong correlations)
        threshold = 0.7
        for i, mag1 in enumerate(mag_names):
            for j, mag2 in enumerate(mag_names):
                if i < j and similarity_matrix.iloc[i, j] > threshold:
                    G.add_edge(mag1, mag2, weight=similarity_matrix.iloc[i, j])
        
        # ========== CREATE FIGURE ==========
        fig = plt.figure(figsize=(16, 12))
        gs = plt.GridSpec(2, 2, figure=fig, 
                         width_ratios=[3, 1],
                         height_ratios=[4, 1],
                         hspace=0.15, wspace=0.15)
        
        # Main network plot
        ax_main = fig.add_subplot(gs[0, 0])
        
        # Calculate layout - group by phylum
        pos = self._calculate_grouped_layout(G, mag_phyla, mag_names)
        
        # Draw edges with varying thickness
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        
        if edges:
            nx.draw_networkx_edges(G, pos, 
                                  width=[w*2 for w in weights],
                                  alpha=0.3, 
                                  edge_color='gray',
                                  ax=ax_main)
        
        # Draw nodes as pie charts
        for node, (x, y) in pos.items():
            # Get functional composition for this MAG
            if node in mag_data.index:
                values = mag_data.loc[node]
                total = values.sum()
                if total > 0:
                    # Normalize values
                    proportions = values / total
                    
                    # Create pie chart
                    start_angle = 0
                    for category, proportion in proportions.items():
                        if proportion > 0 and category in feature_categories:
                            end_angle = start_angle + proportion * 360
                            
                            wedge = Wedge((x, y), 0.03,  # radius
                                        start_angle, end_angle,
                                        facecolor=feature_categories.get(category, '#808080'),
                                        edgecolor='white', linewidth=0.5)
                            ax_main.add_patch(wedge)
                            start_angle = end_angle
                    
                    # Add black border circle
                    circle = Circle((x, y), 0.03, 
                                  fill=False, 
                                  edgecolor='black', 
                                  linewidth=1.5)
                    ax_main.add_patch(circle)
                    
                    # Add phylum color as outer ring
                    phylum = mag_taxonomy[node]['phylum']
                    outer_circle = Circle((x, y), 0.035,
                                        fill=False,
                                        edgecolor=phylum_colors[phylum],
                                        linewidth=3)
                    ax_main.add_patch(outer_circle)
        
        # Set axis properties
        ax_main.set_xlim(-1.2, 1.2)
        ax_main.set_ylim(-1.2, 1.2)
        ax_main.set_aspect('equal')
        ax_main.axis('off')
        ax_main.set_title(f'{data_type} Functional Network', fontsize=16, fontweight='bold')
        
        # ========== BAR PLOT: Total abundances ==========
        ax_bar = fig.add_subplot(gs[0, 1])
        
        # Calculate total abundance per category
        category_totals = mag_data.sum(axis=0)
        categories = list(feature_categories.keys())
        
        # Filter to only categories present in data
        present_categories = [cat for cat in categories if cat in category_totals.index]
        if present_categories:
            values = [category_totals[cat] for cat in present_categories]
            colors = [feature_categories[cat] for cat in present_categories]
            
            bars = ax_bar.bar(range(len(present_categories)), values, 
                            color=colors, edgecolor='black', linewidth=1.5)
            
            ax_bar.set_xticks(range(len(present_categories)))
            ax_bar.set_xticklabels(present_categories, rotation=45, ha='right', fontsize=9)
            ax_bar.set_ylabel('Total Abundance', fontsize=11, fontweight='bold')
            ax_bar.set_title(f'{data_type} Distribution', fontsize=12, fontweight='bold')
            ax_bar.spines['top'].set_visible(False)
            ax_bar.spines['right'].set_visible(False)
            
            # Add value labels on bars
            for bar, val in zip(bars, values):
                height = bar.get_height()
                ax_bar.text(bar.get_x() + bar.get_width()/2., height,
                          f'{int(val)}', ha='center', va='bottom', fontsize=8)
        
        # ========== LEGEND: Phylum colors ==========
        ax_legend1 = fig.add_subplot(gs[1, 0])
        ax_legend1.axis('off')
        
        # Create phylum legend
        phylum_patches = []
        for phylum in unique_phyla[:10]:  # Limit to 10 for space
            n_mags = sum(1 for p in mag_phyla if p == phylum)
            patch = mpatches.Patch(color=phylum_colors[phylum], 
                                  label=f'{phylum} (n={n_mags})')
            phylum_patches.append(patch)
        
        phylum_legend = ax_legend1.legend(handles=phylum_patches, 
                                         loc='center', ncol=5,
                                         title='Phylum (Outer Ring)', 
                                         title_fontsize=11,
                                         fontsize=9, frameon=True)
        
        # ========== LEGEND: Functional categories ==========
        ax_legend2 = fig.add_subplot(gs[1, 1])
        ax_legend2.axis('off')
        
        # Create functional category legend
        func_patches = []
        for category, color in feature_categories.items():
            if category in present_categories:
                # Shorten long names
                label = category if len(category) <= 20 else category[:17] + '...'
                patch = mpatches.Patch(color=color, label=label)
                func_patches.append(patch)
        
        func_legend = ax_legend2.legend(handles=func_patches,
                                       loc='center', ncol=2,
                                       title=f'{data_type} Categories',
                                       title_fontsize=11,
                                       fontsize=9, frameon=True)
        
        # Add main title with statistics
        n_edges = G.number_of_edges()
        avg_degree = sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0
        
        plt.suptitle(f'{data_type} Functional Network Analysis', fontsize=18, fontweight='bold')
        fig.text(0.5, 0.93, 
                f'{len(mag_names)} MAGs | {n_edges} Connections | Avg Degree: {avg_degree:.1f}',
                ha='center', fontsize=11, style='italic')
        
        # Save figure
        output_name = f'functional_network_{data_type.lower()}'
        plt.savefig(self.output_dir / f'{output_name}.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(self.output_dir / f'{output_name}.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / f'{output_name}.svg', format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info(f"Functional network for {data_type} saved successfully!")
    
    def _calculate_grouped_layout(self, G, mag_phyla, mag_names):
        """Calculate network layout grouped by phylum"""
        import networkx as nx
        
        # Group nodes by phylum
        phylum_groups = {}
        for i, mag in enumerate(mag_names):
            phylum = mag_phyla[i]
            if phylum not in phylum_groups:
                phylum_groups[phylum] = []
            phylum_groups[phylum].append(mag)
        
        # Calculate positions for each group
        pos = {}
        n_groups = len(phylum_groups)
        
        # Arrange groups in a circle
        angle_step = 2 * np.pi / n_groups
        group_radius = 0.7
        
        for i, (phylum, nodes) in enumerate(phylum_groups.items()):
            # Calculate center for this group
            angle = i * angle_step
            center_x = group_radius * np.cos(angle)
            center_y = group_radius * np.sin(angle)
            
            # Create subgraph for this group
            subgraph = G.subgraph(nodes)
            
            # Calculate layout for subgraph
            if len(nodes) == 1:
                sub_pos = {nodes[0]: np.array([0, 0])}
            else:
                sub_pos = nx.spring_layout(subgraph, k=0.3, iterations=50)
            
            # Scale and translate positions
            for node, (x, y) in sub_pos.items():
                pos[node] = (center_x + x * 0.3, center_y + y * 0.3)
        
        return pos
    
    def _create_multilayer_functional_network(self, available_data):
        """Create multi-layer network showing all functional annotations"""
        logger.info("Creating multi-layer functional network...")
        
        import networkx as nx
        from matplotlib.patches import Rectangle, FancyBboxPatch
        import matplotlib.patches as mpatches
        
        # Get common MAGs
        common_mags = None
        for data_type, data_df in available_data.items():
            if data_type == 'CAZyme':
                current_mags = set(data_df.index)
            else:
                current_mags = set(data_df.columns)
            
            if common_mags is None:
                common_mags = current_mags
            else:
                common_mags = common_mags & current_mags
        
        common_mags = list(common_mags)[:20]  # Limit for clarity
        
        if len(common_mags) < 3:
            logger.warning("Too few common MAGs for multi-layer network")
            return
        
        # Get taxonomy for common MAGs
        mag_taxonomy = {}
        mag_phyla = []
        
        for mag in common_mags:
            tax_name = self.convert_mag_names(mag)
            if self.taxonomy is not None:
                tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
                if not tax_row.empty:
                    phylum = tax_row.iloc[0].get('phylum', 'Unknown')
                else:
                    phylum = 'Unknown'
            else:
                phylum = 'Unknown'
            
            mag_taxonomy[mag] = phylum
            mag_phyla.append(phylum)
        
        unique_phyla = sorted(list(set(mag_phyla)))
        phylum_colors = self._generate_publication_colors(unique_phyla)
        
        # Create figure
        fig = plt.figure(figsize=(18, 10))
        ax = fig.add_subplot(111)
        
        # Define layer positions
        layer_height = 0.3
        layer_spacing = 0.35
        n_layers = len(available_data)
        
        # Colors for each data type
        layer_colors = {
            'CAZyme': '#FFE5CC',
            'COG': '#E5F2FF',
            'KEGG': '#E5FFE5'
        }
        
        # Draw layers and nodes
        for layer_idx, (data_type, data_df) in enumerate(available_data.items()):
            y_base = layer_idx * layer_spacing
            
            # Draw layer background
            layer_rect = FancyBboxPatch((-1, y_base - 0.05), 2, layer_height,
                                       boxstyle="round,pad=0.02",
                                       facecolor=layer_colors.get(data_type, '#F0F0F0'),
                                       edgecolor='black', linewidth=2,
                                       alpha=0.3)
            ax.add_patch(layer_rect)
            
            # Add layer label
            ax.text(-1.1, y_base + layer_height/2, data_type,
                   fontsize=14, fontweight='bold', va='center', ha='right')
            
            # Position nodes in this layer
            x_spacing = 2.0 / (len(common_mags) + 1)
            
            for i, mag in enumerate(common_mags):
                x = -0.9 + (i + 1) * x_spacing
                y = y_base + layer_height/2
                
                # Draw node
                phylum = mag_taxonomy[mag]
                node_color = phylum_colors[phylum]
                
                circle = plt.Circle((x, y), 0.03,
                                   facecolor=node_color,
                                   edgecolor='black',
                                   linewidth=2,
                                   zorder=10)
                ax.add_patch(circle)
                
                # Add connections between layers
                if layer_idx > 0:
                    prev_y = (layer_idx - 1) * layer_spacing + layer_height/2
                    ax.plot([x, x], [prev_y + 0.03, y - 0.03],
                           'k-', alpha=0.2, linewidth=1, zorder=1)
        
        # Set axis properties
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-0.1, n_layers * layer_spacing)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add legend
        legend_elements = []
        for phylum in unique_phyla:
            n_mags = sum(1 for p in mag_phyla if p == phylum)
            legend_elements.append(plt.scatter([], [], s=100, c=phylum_colors[phylum],
                                              edgecolors='black', linewidth=2,
                                              label=f'{phylum} (n={n_mags})'))
        
        ax.legend(handles=legend_elements, loc='upper right',
                 title='Phylum', title_fontsize=11, fontsize=9)
        
        plt.title('Multi-Layer Functional Network', fontsize=16, fontweight='bold')
        
        # Save
        plt.savefig(self.output_dir / 'functional_network_multilayer.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(self.output_dir / 'functional_network_multilayer.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Multi-layer functional network saved!")    

    def create_integrated_multiomics_view(self):
        """Create integrated view combining COGs, KOs, and CAZymes with taxonomy"""
        logger.info("Creating integrated multi-omics visualization...")
        
        # Check data availability
        has_cog = 'cog' in self.data
        has_ko = 'ko' in self.data
        has_cazyme = 'cazyme' in self.data
        
        if not any([has_cog, has_ko, has_cazyme]) or self.taxonomy is None:
            logger.warning("Insufficient data for integrated view")
            return
        
        # Get common MAGs across datasets
        common_mags = None
        if has_cog:
            common_mags = set(self.data['cog'].columns)
        if has_ko:
            ko_mags = set(self.data['ko'].columns)
            common_mags = ko_mags if common_mags is None else common_mags & ko_mags
        if has_cazyme:
            cazyme_mags = set(self.data['cazyme'].columns)
            common_mags = cazyme_mags if common_mags is None else common_mags & cazyme_mags
        
        if not common_mags:
            logger.warning("No common MAGs across datasets")
            return
        
        common_mags = list(common_mags)[:30]  # Limit for visualization
        
        # Get taxonomy for common MAGs
        mag_taxonomy = {}
        for mag in common_mags:
            tax_name = self.convert_mag_names(mag)
            tax_row = self.taxonomy[self.taxonomy['user_genome'] == tax_name]
            if not tax_row.empty:
                mag_taxonomy[mag] = {
                    'phylum': tax_row.iloc[0].get('phylum', 'Unknown'),
                    'genus': tax_row.iloc[0].get('genus', 'Unknown')
                }
        
        # Sort by taxonomy
        sorted_mags = sorted(mag_taxonomy.keys(), 
                           key=lambda x: (mag_taxonomy[x]['phylum'], mag_taxonomy[x]['genus']))
        
        # Create figure
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 0.5], 
                     width_ratios=[0.5, 8, 0.5], hspace=0.15, wspace=0.02)
        
        # Shared taxonomy sidebar
        ax_tax = fig.add_subplot(gs[:3, 0])
        phyla = [mag_taxonomy[mag]['phylum'] for mag in sorted_mags]
        
        for i, (mag, phylum) in enumerate(zip(sorted_mags, phyla)):
            color = PHYLUM_COLORS.get(phylum, '#B09C85')
            rect = Rectangle((0, i-0.4), 1, 0.8, facecolor=color, edgecolor='black', linewidth=0.5)
            ax_tax.add_patch(rect)
        
        ax_tax.set_xlim(0, 1)
        ax_tax.set_ylim(-0.5, len(sorted_mags)-0.5)
        ax_tax.set_title('Taxonomy', fontsize=10, rotation=90, x=-0.3, y=0.5)
        ax_tax.axis('off')
        
        # Panel 1: COG categories
        if has_cog:
            ax_cog = fig.add_subplot(gs[0, 1])
            cog_data = self.data['cog'][sorted_mags]
            
            # Group by main categories
            cog_profiles = []
            for main_cat, cogs in COG_MAIN_CATEGORIES.items():
                cat_sum = cog_data.loc[cog_data.index.isin(cogs)].sum(axis=0)
                cog_profiles.append(cat_sum)
            
            cog_matrix = np.array(cog_profiles)
            cog_norm = cog_matrix / cog_matrix.sum(axis=0)
            
            im1 = ax_cog.imshow(cog_norm, aspect='auto', cmap='Blues', vmin=0, vmax=0.5)
            ax_cog.set_yticks(range(len(COG_MAIN_CATEGORIES)))
            ax_cog.set_yticklabels(list(COG_MAIN_CATEGORIES.keys()), fontsize=7)
            ax_cog.set_xticks([])
            ax_cog.set_title('COG Functional Categories', fontsize=10, fontweight='bold')
            
            # Colorbar
            ax_cbar1 = fig.add_subplot(gs[0, 2])
            cbar1 = plt.colorbar(im1, cax=ax_cbar1)
            cbar1.set_label('Proportion', fontsize=8)
            cbar1.ax.tick_params(labelsize=7)
        
        # Panel 2: Top KO pathways
        if has_ko:
            ax_ko = fig.add_subplot(gs[1, 1])
            ko_data = self.data['ko'][sorted_mags]
            
            # Get top variable KOs
            ko_var = ko_data.var(axis=1)
            top_kos = ko_var.nlargest(20).index
            
            ko_matrix = ko_data.loc[top_kos].values
            ko_norm = np.log10(ko_matrix + 1)
            
            im2 = ax_ko.imshow(ko_norm, aspect='auto', cmap='Greens', vmin=0)
            ax_ko.set_yticks(range(len(top_kos)))
            ax_ko.set_yticklabels(top_kos, fontsize=6)
            ax_ko.set_xticks([])
            ax_ko.set_title('KEGG Orthology (Top 20)', fontsize=10, fontweight='bold')
            
            # Colorbar
            ax_cbar2 = fig.add_subplot(gs[1, 2])
            cbar2 = plt.colorbar(im2, cax=ax_cbar2)
            cbar2.set_label('log10(Count+1)', fontsize=8)
            cbar2.ax.tick_params(labelsize=7)
        
        # Panel 3: CAZyme families
        if has_cazyme:
            ax_cazyme = fig.add_subplot(gs[2, 1])
            cazyme_data = self.data['cazyme'][sorted_mags]
            
            # Get top CAZyme families
            cazyme_sum = cazyme_data.sum(axis=1)
            top_cazymes = cazyme_sum.nlargest(15).index
            
            cazyme_matrix = cazyme_data.loc[top_cazymes].values
            cazyme_norm = np.log10(cazyme_matrix + 1)
            
            im3 = ax_cazyme.imshow(cazyme_norm, aspect='auto', cmap='Oranges', vmin=0)
            ax_cazyme.set_yticks(range(len(top_cazymes)))
            ax_cazyme.set_yticklabels(top_cazymes, fontsize=7)
            ax_cazyme.set_xticks(range(len(sorted_mags)))
            ax_cazyme.set_xticklabels([m.split('_')[0] for m in sorted_mags], 
                                      rotation=90, fontsize=6)
            ax_cazyme.set_title('CAZyme Families (Top 15)', fontsize=10, fontweight='bold')
            ax_cazyme.set_xlabel('MAGs', fontsize=9)
            
            # Colorbar
            ax_cbar3 = fig.add_subplot(gs[2, 2])
            cbar3 = plt.colorbar(im3, cax=ax_cbar3)
            cbar3.set_label('log10(Count+1)', fontsize=8)
            cbar3.ax.tick_params(labelsize=7)
        
        # Legend for taxonomy
        ax_legend = fig.add_subplot(gs[3, :])
        ax_legend.axis('off')
        
        unique_phyla = list(set(phyla))[:10]
        patches = [mpatches.Patch(color=PHYLUM_COLORS.get(p, '#B09C85'), label=p) 
                  for p in unique_phyla]
        ax_legend.legend(handles=patches, loc='center', ncol=5, fontsize=8, 
                        title='Phylum', title_fontsize=9)
        
        plt.suptitle('Integrated Multi-Omics Functional Profile', 
                    fontsize=14, fontweight='bold')
        
        plt.savefig(self.output_dir / 'integrated_multiomics_view.pdf', dpi=600)
        plt.savefig(self.output_dir / 'integrated_multiomics_view.png', dpi=300)
        plt.close()
        
        logger.info("Integrated multi-omics view saved")
    
    def run_all_analyses(self):
        """Run all publication-quality analyses"""
        logger.info("="*60)
        logger.info("Starting publication-quality visualization pipeline")
        logger.info("="*60)
        
        # Load all data
        self.load_all_data()
        
        # Create visualizations
        visualizations = [
            ('Taxonomy-CAZyme heatmap', self.create_taxonomy_cazyme_heatmap),
            ('COG Sankey diagram', self.create_cog_sankey_diagram),
            ('KEGG Sankey diagram', self.create_kegg_level_sankey),
            ('KEGG Metabolism Sankey', self.create_kegg_metabolism_sankey),
            ('Functional PCA analysis', self.create_functional_pca_analysis),
            ('Functional network visualization', self.create_functional_network_visualization),
            ('Enhanced network analysis', self.create_enhanced_network_analysis),
            ('Integrated multi-omics view', self.create_integrated_multiomics_view)
        ]
        
        for name, func in visualizations:
            try:
                func()
                logger.info(f"? {name} completed")
            except Exception as e:
                logger.error(f"? Failed to create {name}: {e}")
                import traceback
                traceback.print_exc()
        
        logger.info("="*60)
        logger.info("Publication-quality visualizations complete!")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info("="*60)




# ============================================================================
# PIPELINE INTEGRATION
# ============================================================================

def run(cpus, memory, time, project_input, project_output, 
        mag_lists=None, list_names=None, viz_output_suffix="", viz_types=None):
    """
    Run publication-quality visualizations as a pipeline step
    """
    logger.info("="*60)
    logger.info("Starting Publication-Quality Visualization Pipeline")
    logger.info("="*60)
    
    project_dir = Path(project_output)
    
    try:
        viz = PublicationQualityVisualizations(
            str(project_dir),
            viz_output_suffix
        )
        
        viz.run_all_analyses()
        
        return True
        
    except Exception as e:
        logger.error(f"Error in visualization pipeline: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# ADD THIS SECTION TO THE END OF YOUR EXISTING advanced_visualizations.py
# ============================================================================

# Add novel MAGs filtering and visualization methods to the existing class
def _load_novel_mags_list(self):
    """Load list of novel MAGs from UniqueMags directory"""
    novel_mags_dir = self.project_dir / "Novel_Mags" / "UniqueMags"
    novel_mags_list = []
    novel_mags_functional_names = []
    
    if novel_mags_dir.exists():
        logger.info(f"Loading novel MAGs from: {novel_mags_dir}")
        
        # Get all MAG files in the directory
        mag_files = list(novel_mags_dir.glob("*.fa")) + list(novel_mags_dir.glob("*.fasta"))
        
        for mag_file in mag_files:
            mag_name = mag_file.stem
            novel_mags_list.append(mag_name)
            
            # Create variations for matching
            for suffix in ['.counts', '_fa.counts', '_fa_faa', '_faa', '_fa', '.faa', '.fa']:
                novel_mags_functional_names.append(mag_name + suffix)
            novel_mags_functional_names.append(mag_name)
        
        logger.info(f"Found {len(novel_mags_list)} novel MAGs")
    else:
        logger.warning(f"Novel MAGs directory not found: {novel_mags_dir}")
    
    self.novel_mags_list = novel_mags_list
    self.novel_mags_functional_names = novel_mags_functional_names
    return novel_mags_list

def _filter_data_for_novel_mags(self):
    """Filter all loaded data to include only novel MAGs"""
    if not hasattr(self, 'novel_mags_list') or not self.novel_mags_list:
        return None
    
    filtered_data = {}
    
    # Filter each data type
    if 'cog' in self.data:
        cog_data = self.data['cog']
        novel_cols = [col for col in cog_data.columns 
                     if self.convert_mag_names(col) in self.novel_mags_list 
                     or col in self.novel_mags_functional_names]
        if novel_cols:
            filtered_data['cog'] = cog_data[novel_cols]
            logger.info(f"Filtered COG: {len(novel_cols)}/{len(cog_data.columns)} MAGs")
    
    if 'ko' in self.data:
        ko_data = self.data['ko']
        novel_cols = [col for col in ko_data.columns 
                     if self.convert_mag_names(col) in self.novel_mags_list 
                     or col in self.novel_mags_functional_names]
        if novel_cols:
            filtered_data['ko'] = ko_data[novel_cols]
            logger.info(f"Filtered KO: {len(novel_cols)}/{len(ko_data.columns)} MAGs")
    
    if 'cazyme' in self.data:
        cazyme_data = self.data['cazyme']
        novel_rows = [idx for idx in cazyme_data.index 
                     if self.convert_mag_names(idx) in self.novel_mags_list 
                     or idx in self.novel_mags_functional_names]
        if novel_rows:
            filtered_data['cazyme'] = cazyme_data.loc[novel_rows]
            logger.info(f"Filtered CAZyme: {len(novel_rows)}/{len(cazyme_data.index)} MAGs")
    
    if 'kegg_pathways' in self.data:
        kegg_data = self.data['kegg_pathways']
        numeric_cols = kegg_data.select_dtypes(include=[np.number]).columns
        novel_cols = [col for col in numeric_cols 
                     if self.convert_mag_names(str(col)) in self.novel_mags_list 
                     or str(col) in self.novel_mags_functional_names]
        if novel_cols:
            meta_cols = kegg_data.select_dtypes(exclude=[np.number]).columns.tolist()
            filtered_data['kegg_pathways'] = kegg_data[meta_cols + novel_cols]
            logger.info(f"Filtered KEGG: {len(novel_cols)} novel MAG columns")
    
    # Filter taxonomy
    if self.taxonomy is not None:
        novel_taxonomy_rows = [idx for idx, row in self.taxonomy.iterrows() 
                              if row['user_genome'] in self.novel_mags_list]
        if novel_taxonomy_rows:
            self.novel_taxonomy = self.taxonomy.loc[novel_taxonomy_rows].copy()
            logger.info(f"Filtered taxonomy: {len(novel_taxonomy_rows)} MAGs")
    
    return filtered_data

def _create_novel_mags_visualizations(self):
    """Create all visualizations for novel MAGs only"""
    logger.info("="*60)
    logger.info("Creating Novel MAGs-specific visualizations")
    logger.info("="*60)
    
    # Load novel MAGs
    novel_mags = self._load_novel_mags_list()
    if not novel_mags:
        logger.warning("No novel MAGs found")
        return
    
    # Setup output directory
    novel_output_dir = self.output_dir / "Novel_MAGs_Only"
    novel_output_dir.mkdir(exist_ok=True)
    
    # Backup original data
    original_data = self.data.copy()
    original_taxonomy = self.taxonomy
    original_output_dir = self.output_dir
    
    # Filter for novel MAGs
    filtered_data = self._filter_data_for_novel_mags()
    if not filtered_data:
        logger.warning("No data after filtering")
        return
    
    # Use filtered data
    self.data = filtered_data
    if hasattr(self, 'novel_taxonomy'):
        self.taxonomy = self.novel_taxonomy
    self.output_dir = novel_output_dir
    
    # Save filtered data
    self._save_novel_mags_data()
    
    # Run all visualizations with novel MAGs data
    visualizations = [
        ('CAZyme Heatmap', self.create_taxonomy_cazyme_heatmap),
        ('COG Sankey', self.create_cog_sankey_diagram),
        ('KEGG Sankey', self.create_kegg_level_sankey),
        ('KEGG Metabolism', self.create_kegg_metabolism_sankey),
        ('Functional PCA', self.create_functional_pca_analysis),
        ('Network', self.create_functional_network_visualization),
        ('Enhanced Network', self.create_enhanced_network_analysis),
        ('Multi-omics', self.create_integrated_multiomics_view)
    ]
    
    for name, method in visualizations:
        try:
            logger.info(f"Creating Novel MAGs {name}...")
            method()
            logger.info(f"? {name} completed")
        except Exception as e:
            logger.error(f"? Failed {name}: {e}")
    
    # Create summary
    self._create_novel_mags_summary()
    
    # Restore original data
    self.data = original_data
    self.taxonomy = original_taxonomy
    self.output_dir = original_output_dir
    
    logger.info(f"Novel MAGs visualizations saved in: {novel_output_dir}")

def _save_novel_mags_data(self):
    """Save filtered novel MAGs data"""
    data_dir = self.output_dir / "filtered_data"
    data_dir.mkdir(exist_ok=True)
    
    if 'cog' in self.data:
        self.data['cog'].to_excel(data_dir / "novel_mags_cog.xlsx")
    if 'ko' in self.data:
        self.data['ko'].to_excel(data_dir / "novel_mags_ko.xlsx")
    if 'cazyme' in self.data:
        self.data['cazyme'].to_excel(data_dir / "novel_mags_cazyme.xlsx")
    if 'kegg_pathways' in self.data:
        self.data['kegg_pathways'].to_excel(data_dir / "novel_mags_kegg.xlsx", index=False)
    
    with open(data_dir / "novel_mags_list.txt", 'w') as f:
        for mag in self.novel_mags_list:
            f.write(f"{mag}\n")

def _create_novel_mags_summary(self):
    """Create summary for novel MAGs"""
    with open(self.output_dir / "novel_mags_summary.txt", 'w') as f:
        f.write("NOVEL MAGs ANALYSIS SUMMARY\n")
        f.write("="*60 + "\n")
        f.write(f"Total Novel MAGs: {len(self.novel_mags_list)}\n\n")
        
        if 'cog' in self.data:
            f.write(f"COG: {self.data['cog'].shape}\n")
        if 'ko' in self.data:
            f.write(f"KO: {self.data['ko'].shape}\n")
        if 'cazyme' in self.data:
            f.write(f"CAZyme: {self.data['cazyme'].shape}\n")
        
        if hasattr(self, 'novel_taxonomy') and self.novel_taxonomy is not None:
            f.write(f"\nTaxonomy available for {len(self.novel_taxonomy)} MAGs\n")
            if 'phylum' in self.novel_taxonomy.columns:
                f.write("\nPhylum distribution:\n")
                for phylum, count in self.novel_taxonomy['phylum'].value_counts().items():
                    f.write(f"  {phylum}: {count}\n")

# Add methods to the class
PublicationQualityVisualizations._load_novel_mags_list = _load_novel_mags_list
PublicationQualityVisualizations._filter_data_for_novel_mags = _filter_data_for_novel_mags
PublicationQualityVisualizations._create_novel_mags_visualizations = _create_novel_mags_visualizations
PublicationQualityVisualizations._save_novel_mags_data = _save_novel_mags_data
PublicationQualityVisualizations._create_novel_mags_summary = _create_novel_mags_summary

# Override run_all_analyses to include both versions
_original_run_all = PublicationQualityVisualizations.run_all_analyses

def enhanced_run_all_analyses(self):
    """Run analyses for BOTH all MAGs and novel MAGs"""
    # Phase 1: All MAGs (original)
    logger.info("="*70)
    logger.info("PHASE 1: Creating visualizations for ALL MAGs")
    logger.info("="*70)
    _original_run_all(self)
    
    # Phase 2: Novel MAGs only
    logger.info("\n" + "="*70)
    logger.info("PHASE 2: Creating visualizations for NOVEL MAGs only")
    logger.info("="*70)
    self._create_novel_mags_visualizations()
    
    # Combined summary
    logger.info("\n" + "="*70)
    logger.info("ANALYSIS COMPLETE")
    logger.info(f"All MAGs visualizations: {self.output_dir}")
    logger.info(f"Novel MAGs visualizations: {self.output_dir}/Novel_MAGs_Only")
    logger.info("="*70)

PublicationQualityVisualizations.run_all_analyses = enhanced_run_all_analyses



if __name__ == "__main__":
    """Standalone execution"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Publication-Quality Visualizations for Metagenomics"
    )
    
    parser.add_argument("project_dir", help="Project directory path")
    parser.add_argument("--suffix", default="", help="Output suffix")
    
    args = parser.parse_args()
    
    viz = PublicationQualityVisualizations(args.project_dir, args.suffix)
    viz.run_all_analyses()
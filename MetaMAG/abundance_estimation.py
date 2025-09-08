#!/usr/bin/env python3
"""
Comprehensive abundance estimation module for metagenomic samples.
IMPROVED VERSION: Smart checkpointing - skips already completed steps
"""

import os
import subprocess
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu, entropy
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

from MetaMAG.config import config, get_tool_path, get_tool_command

# Set publication-quality style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'grid.alpha': 0.3
})

def check_kraken_outputs_exist(sample, output_dir):
    """Check if Kraken2 outputs already exist for a sample"""
    kraken_output = os.path.join(output_dir, f"{sample}_kraken.txt")
    kraken_report = os.path.join(output_dir, f"{sample}_kreport.txt")
    
    # Check both files exist and are non-empty
    if os.path.exists(kraken_output) and os.path.exists(kraken_report):
        if os.path.getsize(kraken_output) > 0 and os.path.getsize(kraken_report) > 0:
            return True
    return False

def check_bracken_outputs_exist(sample, output_dir, taxonomic_level):
    """Check if Bracken outputs already exist for a sample at specific level"""
    bracken_output = os.path.join(output_dir, f"{sample}_bracken_{taxonomic_level}.txt")
    
    # Check file exists and is non-empty
    if os.path.exists(bracken_output) and os.path.getsize(bracken_output) > 0:
        return True
    return False

def check_merged_output_exists(output_dir, taxonomic_level):
    """Check if merged Bracken output already exists for a taxonomic level"""
    merged_output = os.path.join(output_dir, f"merged_bracken_{taxonomic_level}.txt")
    
    if os.path.exists(merged_output) and os.path.getsize(merged_output) > 0:
        return True
    return False

def run_comprehensive_abundance(samples, cpus, memory, time, project_input, project_output, 
                               kraken_db_path, taxonomic_levels=['S', 'G'], read_length=150, 
                               threshold=10, metadata_file=None, force_rerun=False):
    """
    Run comprehensive abundance estimation with SMART CHECKPOINTING
    
    Args:
        force_rerun: If True, ignore existing outputs and rerun everything
    """
    print(f"\n[INFO] Starting comprehensive abundance estimation for {len(samples)} samples")
    print(f"[INFO] Smart checkpointing {'DISABLED (force rerun)' if force_rerun else 'ENABLED'}")
    print(f"[INFO] Taxonomic levels: {', '.join(taxonomic_levels)}")
    print(f"[INFO] Kraken database path: {kraken_db_path}")
    
    # Verify kraken_db_path
    if not kraken_db_path or kraken_db_path == "None":
        print("[ERROR] Kraken database path is not provided or is None!")
        return False
    
    if not os.path.exists(kraken_db_path):
        print(f"[ERROR] Kraken database path does not exist: {kraken_db_path}")
        return False
    
    # Check metadata availability
    metadata_available = bool(metadata_file and os.path.exists(metadata_file))
    print(f"[INFO] Metadata file {'provided' if metadata_available else 'not provided'}")
    
    # Step 1: Run Kraken2 abundance estimation (with checkpointing)
    print("\n" + "="*50)
    print("STEP 1: Kraken2 abundance estimation")
    print("="*50)
    
    success = run_kraken2_abundance_smart(samples, cpus, memory, time, project_input, 
                                         project_output, kraken_db_path, force_rerun)
    if not success:
        print("[ERROR] Kraken2 abundance estimation failed")
        return False
    
    # Step 2: Run Bracken for each taxonomic level (with checkpointing)
    print("\n" + "="*50)
    print("STEP 2: Bracken abundance estimation")
    print("="*50)
    
    for level in taxonomic_levels:
        print(f"\n[INFO] Processing taxonomic level: {level}")
        success = run_bracken_abundance_smart(samples, cpus, project_output, kraken_db_path, 
                                             level, read_length, threshold, force_rerun)
        if not success:
            print(f"[WARNING] Bracken estimation failed for level {level}")
            continue
        
        # Step 3: Combine Bracken outputs (with checkpointing)
        print(f"\n[INFO] Combining Bracken outputs for level {level}")
        success = combine_bracken_outputs_smart(samples, project_output, level, force_rerun)
        if not success:
            print(f"[WARNING] Failed to combine Bracken outputs for level {level}")
    
    # Step 4: Generate comprehensive visualizations
    print("\n" + "="*50)
    print("STEP 3: Generating abundance visualizations")
    print("="*50)
    
    success = generate_abundance_plots(project_output, taxonomic_levels, metadata_file, metadata_available)
    if not success:
        print("[ERROR] Failed to generate abundance plots")
        return False
    
    print("\n" + "="*50)
    print("SUCCESS: Comprehensive abundance analysis completed!")
    print("="*50)
    print(f"Results saved in: {project_output}")
    print("- Kraken_Abundance/ - Raw Kraken outputs")
    print("- Bracken_Abundance_*/ - Bracken estimates by taxonomic level")
    print("- Merged_Bracken_Outputs/ - Combined abundance matrices")
    print("- Abundance_Plots/ - Publication-quality visualizations")
    if metadata_available:
        print("- abundance_summary.html - Interactive dashboard with group comparisons")
    else:
        print("- abundance_summary.html - Interactive dashboard with sample-based analysis")
    
    return True


def run_kraken2_abundance_smart(samples, cpus, memory, time, project_input, project_output, 
                               kraken_db_path, force_rerun=False):
    """Run Kraken2 with SMART CHECKPOINTING - skip if outputs exist"""
    # Get paths from config
    perl_path = config["environment"].get("PERL5LIB", "/opt/ghpc/lib64/perl-5.36")
    kraken2_path = get_tool_path("kraken2")
    if not kraken2_path:
        print(f"[ERROR] Kraken2 tool not found in configuration")
        return False
    
    print(f"[INFO] Using PERL5LIB: {perl_path}")
    print(f"[INFO] Using Kraken2 path: {kraken2_path}")
    print(f"[INFO] Using Kraken database: {kraken_db_path}")
    
    # Verify kraken_db_path
    if not kraken_db_path or kraken_db_path.lower() == "none":
        print("[ERROR] Kraken database path is None or undefined!")
        return False
    
    # Set up paths
    host_removal_dir = os.path.join(project_output, "Host_Removal")
    output_dir = os.path.join(project_output, "Kraken_Abundance")
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if host removal directory exists
    if not os.path.exists(host_removal_dir):
        print(f"[ERROR] Host removal directory not found: {host_removal_dir}")
        print("[ERROR] Please run the 'host_removal' step first!")
        return False
    
    success_count = 0
    skipped_count = 0
    failed_count = 0
    
    for sample in samples:
        # Check if outputs already exist
        if not force_rerun and check_kraken_outputs_exist(sample, output_dir):
            print(f"[SKIP] Kraken2 outputs already exist for sample: {sample}")
            skipped_count += 1
            success_count += 1  # Count as success for downstream processing
            continue
        
        print(f"[INFO] Running Kraken2 for sample: {sample}")
        
        # Input files
        r1_file = os.path.join(host_removal_dir, f"{sample}_unmapped_R1.fastq.gz")
        r2_file = os.path.join(host_removal_dir, f"{sample}_unmapped_R2.fastq.gz")
        
        if not os.path.exists(r1_file) or not os.path.exists(r2_file):
            print(f"[WARNING] Input files not found for sample {sample}")
            print(f"         Expected: {r1_file}")
            print(f"         Expected: {r2_file}")
            failed_count += 1
            continue
        
        # Output files
        kraken_output = os.path.join(output_dir, f"{sample}_kraken.txt")
        kraken_report = os.path.join(output_dir, f"{sample}_kreport.txt")
        
        # Construct Kraken2 command
        cmd = f"PERL5LIB={perl_path} {kraken2_path} --paired --gzip-compressed " \
              f"--db \"{kraken_db_path}\" --threads {cpus} --report \"{kraken_report}\" " \
              f"--output \"{kraken_output}\" \"{r1_file}\" \"{r2_file}\""
        
        print(f"[DEBUG] Running command: {cmd}")
        
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            print(f"[SUCCESS] Kraken2 completed for sample {sample}")
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Kraken2 failed for sample {sample}")
            print(f"[ERROR] Exit code: {e.returncode}")
            print(f"[ERROR] stderr: {e.stderr}")
            failed_count += 1
    
    print(f"\n[SUMMARY] Kraken2 Results:")
    print(f"  - Successfully processed: {success_count}/{len(samples)} samples")
    print(f"  - Skipped (already done): {skipped_count}/{len(samples)} samples")
    print(f"  - Failed: {failed_count}/{len(samples)} samples")
    
    return success_count > 0


def run_bracken_abundance_smart(samples, cpus, project_output, kraken_db_path, taxonomic_level, 
                               read_length, threshold, force_rerun=False):
    """Run Bracken with SMART CHECKPOINTING - skip if outputs exist"""
    bracken_path = get_tool_path("bracken")
    if not bracken_path:
        print(f"[ERROR] Bracken tool not found in configuration")
        return False
    
    # Set up paths
    kraken_output_dir = os.path.join(project_output, "Kraken_Abundance")
    bracken_output_dir = os.path.join(project_output, f"Bracken_Abundance_{taxonomic_level}")
    os.makedirs(bracken_output_dir, exist_ok=True)
    
    # Define kmer distribution file
    kmer_distrib_file = os.path.join(kraken_db_path, f"database{read_length}mers.kmer_distrib")
    
    if not os.path.exists(kmer_distrib_file):
        print(f"[ERROR] Kmer distribution file not found: {kmer_distrib_file}")
        print("[ERROR] This file is required for Bracken")
        return False
    
    success_count = 0
    skipped_count = 0
    failed_count = 0
    
    for sample in samples:
        # Check if output already exists
        if not force_rerun and check_bracken_outputs_exist(sample, bracken_output_dir, taxonomic_level):
            print(f"[SKIP] Bracken output already exists for sample {sample} at level {taxonomic_level}")
            skipped_count += 1
            success_count += 1  # Count as success for downstream
            continue
        
        # Input and output files
        kraken_report = os.path.join(kraken_output_dir, f"{sample}_kreport.txt")
        bracken_output = os.path.join(bracken_output_dir, f"{sample}_bracken_{taxonomic_level}.txt")
        
        if not os.path.exists(kraken_report):
            print(f"[WARNING] Kraken report not found for sample {sample}: {kraken_report}")
            failed_count += 1
            continue
        
        print(f"[INFO] Running Bracken for sample {sample} at level {taxonomic_level}")
        
        # Run Bracken
        cmd = f"python3 {bracken_path} -i {kraken_report} -k {kmer_distrib_file} " \
              f"-l {taxonomic_level} -t {threshold} -o {bracken_output}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"[SUCCESS] Bracken completed for sample {sample}")
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"[WARNING] Bracken failed for sample {sample}: {e}")
            failed_count += 1
    
    print(f"\n[SUMMARY] Bracken Results for level {taxonomic_level}:")
    print(f"  - Successfully processed: {success_count}/{len(samples)} samples")
    print(f"  - Skipped (already done): {skipped_count}/{len(samples)} samples")
    print(f"  - Failed: {failed_count}/{len(samples)} samples")
    
    return success_count > 0


def combine_bracken_outputs_smart(samples, project_output, taxonomic_level, force_rerun=False):
    """Combine Bracken outputs with SMART CHECKPOINTING"""
    combine_script = os.path.join(os.path.dirname(__file__), "combine_bracken_outputs.py")
    
    # Set up paths
    bracken_output_dir = os.path.join(project_output, f"Bracken_Abundance_{taxonomic_level}")
    merged_output_dir = os.path.join(project_output, "Merged_Bracken_Outputs")
    os.makedirs(merged_output_dir, exist_ok=True)
    
    # Check if merged output already exists
    if not force_rerun and check_merged_output_exists(merged_output_dir, taxonomic_level):
        print(f"[SKIP] Merged output already exists for level {taxonomic_level}")
        return True
    
    print(f"[INFO] Merging Bracken outputs for level {taxonomic_level}")
    
    # Collect Bracken output files
    bracken_files = []
    valid_samples = []
    
    for sample in samples:
        bracken_file = os.path.join(bracken_output_dir, f"{sample}_bracken_{taxonomic_level}.txt")
        if os.path.exists(bracken_file):
            bracken_files.append(bracken_file)
            valid_samples.append(sample)
    
    if not bracken_files:
        print(f"[WARNING] No valid Bracken output files found for level {taxonomic_level}")
        return False
    
    print(f"[INFO] Found {len(bracken_files)} Bracken files to merge")
    
    # Define output file
    merged_output = os.path.join(merged_output_dir, f"merged_bracken_{taxonomic_level}.txt")
    
    # Construct combine command
    files_arg = " ".join(bracken_files)
    names_arg = ",".join(valid_samples)
    
    # Use the standalone combine script
    cmd = f"python3 {combine_script} --files {files_arg} --names {names_arg} -o {merged_output}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"[SUCCESS] Combined Bracken outputs for level {taxonomic_level}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] Failed to combine Bracken outputs: {e}")
        return False


# Keep all the other functions unchanged (they don't need checkpointing):
# - prepare_sample_metadata
# - create_sample_groups_from_names  
# - generate_abundance_plots
# - load_bracken_data
# - plot_* functions
# - create_summary_dashboard

def prepare_sample_metadata(abundance_matrix, metadata, metadata_available):
    """Prepare sample metadata mapping for both cases: with and without metadata"""
    if metadata_available and metadata is not None:
        if 'Sample' in metadata.columns and 'Treatment' in metadata.columns:
            print("[INFO] Using provided metadata with treatment groups")
            sample_to_treatment = dict(zip(metadata['Sample'], metadata['Treatment']))
            # Fill in any missing samples with "Unknown"
            for sample in abundance_matrix.index:
                if sample not in sample_to_treatment:
                    sample_to_treatment[sample] = "Unknown"
        else:
            print("[WARNING] Metadata file missing 'Sample' or 'Treatment' columns")
            print("[INFO] Creating groups based on sample names")
            sample_to_treatment = create_sample_groups_from_names(abundance_matrix.index)
    else:
        print("[INFO] No metadata provided - creating groups based on sample names")
        sample_to_treatment = create_sample_groups_from_names(abundance_matrix.index)
    
    return sample_to_treatment


def create_sample_groups_from_names(sample_names):
    """Create groups based on sample naming patterns when no metadata is available"""
    sample_to_treatment = {}
    
    # Try to extract patterns from sample names
    for sample in sample_names:
        # Pattern 1: Extract prefix before first number or underscore
        parts = sample.split('_')
        if len(parts) > 1:
            # Use first part as group
            group = parts[0]
        else:
            # Pattern 2: Extract alphabetic prefix from name with numbers
            import re
            match = re.match(r'^([A-Za-z]+)', sample)
            if match:
                group = match.group(1)
            else:
                # Pattern 3: Create groups based on position
                position = len(sample_to_treatment)
                group = f"Group_{(position // 4) + 1}"
        
        sample_to_treatment[sample] = group
    
    # If too many groups (>10), create fewer groups based on alphabetical sorting
    unique_groups = list(set(sample_to_treatment.values()))
    if len(unique_groups) > 10:
        print(f"[INFO] Too many groups detected ({len(unique_groups)}), creating fewer groups")
        sorted_samples = sorted(sample_names)
        samples_per_group = max(3, len(sorted_samples) // 5)  # Max 5 groups, min 3 samples per group
        
        sample_to_treatment = {}
        for i, sample in enumerate(sorted_samples):
            group_num = (i // samples_per_group) + 1
            sample_to_treatment[sample] = f"Group_{group_num}"
    
    return sample_to_treatment


def generate_abundance_plots(project_output, taxonomic_levels, metadata_file=None, metadata_available=False):
    """Generate comprehensive abundance visualizations"""
    output_dir = os.path.join(project_output, "Abundance_Plots")
    os.makedirs(output_dir, exist_ok=True)
    
    # Load metadata if provided
    metadata = None
    if metadata_available and metadata_file:
        try:
            metadata = pd.read_csv(metadata_file)
            print(f"[INFO] Loaded metadata with {len(metadata)} entries")
        except Exception as e:
            print(f"[WARNING] Failed to load metadata file: {e}")
            metadata_available = False
    
    # Generate plots for each taxonomic level
    for level in taxonomic_levels:
        try:
            print(f"[INFO] Generating plots for {level} level")
            
            # Load Bracken data
            abundance_matrix = load_bracken_data(project_output, level)
            level_name = {'S': 'Species', 'G': 'Genus', 'F': 'Family', 
                         'O': 'Order', 'C': 'Class', 'P': 'Phylum'}.get(level, level)
            
            # Prepare sample metadata
            sample_to_treatment = prepare_sample_metadata(abundance_matrix, metadata, metadata_available)
            
            # Generate all plots - ACTUALLY CALL THE FUNCTIONS
            print(f"[INFO] Creating visualizations for {level_name} level...")
            plot_top_taxa_barplot(abundance_matrix, level_name, output_dir)
            plot_abundance_heatmap(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available)
            plot_pca_analysis(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available)
            plot_diversity_metrics(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available)
            plot_stacked_composition(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available)
            
            # Additional plots for when we have groups
            if metadata_available or len(set(sample_to_treatment.values())) > 1:
                plot_group_comparison_boxplots(abundance_matrix, level_name, output_dir, sample_to_treatment)
                plot_beta_diversity_between_groups(abundance_matrix, level_name, output_dir, sample_to_treatment)
            
        except Exception as e:
            print(f"[WARNING] Failed to generate plots for level {level}: {e}")
            import traceback
            traceback.print_exc()
    
    # Create summary dashboard
    create_summary_dashboard(project_output, taxonomic_levels, metadata, metadata_available)
    
    print(f"[INFO] All plots saved to: {output_dir}")
    return True


def load_bracken_data(project_output, taxonomic_level):
    """Load and process Bracken data"""
    merged_dir = os.path.join(project_output, "Merged_Bracken_Outputs")
    merged_file = os.path.join(merged_dir, f"merged_bracken_{taxonomic_level}.txt")
    
    if not os.path.exists(merged_file):
        raise FileNotFoundError(f"Merged Bracken file not found: {merged_file}")
    
    # Load data
    data = pd.read_csv(merged_file, sep='\t')
    
    # Extract fraction columns and reorganize
    frac_cols = [col for col in data.columns if col.endswith('_frac')]
    sample_names = [col.replace('_frac', '') for col in frac_cols]
    
    # Create abundance matrix
    abundance_matrix = data[frac_cols].T
    abundance_matrix.columns = data['name']
    abundance_matrix.index = sample_names
    
    # Clean taxonomic names
    cleaned_columns = []
    for col in abundance_matrix.columns:
        cleaned = col.split('__')[-1] if '__' in col else col
        cleaned_columns.append(cleaned)
    abundance_matrix.columns = cleaned_columns
    
    return abundance_matrix.fillna(0)


def plot_top_taxa_barplot(abundance_matrix, level_name, output_dir, top_n=20):
    """Plot top abundant taxa as barplot"""
    print(f"  - Creating top {top_n} taxa barplot...")
    
    # Get top taxa by mean abundance
    mean_abundance = abundance_matrix.mean(axis=0).sort_values(ascending=False)
    top_taxa = mean_abundance.head(top_n)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create barplot
    bars = ax.bar(range(len(top_taxa)), top_taxa.values * 100)
    
    # Color bars by abundance
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_taxa)))
    for bar, color in zip(bars, colors):
        bar.set_color(color)
    
    # Customize plot
    ax.set_xticks(range(len(top_taxa)))
    ax.set_xticklabels(top_taxa.index, rotation=45, ha='right')
    ax.set_xlabel(f'{level_name} Name', fontsize=12)
    ax.set_ylabel('Mean Relative Abundance (%)', fontsize=12)
    ax.set_title(f'Top {top_n} Most Abundant {level_name}', fontsize=14, fontweight='bold')
    
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, top_taxa.values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{val*100:.1f}%', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'top_{top_n}_{level_name.lower()}_barplot.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_abundance_heatmap(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available):
    """Plot abundance heatmap with hierarchical clustering"""
    print(f"  - Creating abundance heatmap...")
    
    # Select top taxa for visualization
    top_n = min(50, len(abundance_matrix.columns))
    mean_abundance = abundance_matrix.mean(axis=0).sort_values(ascending=False)
    top_taxa = mean_abundance.head(top_n).index
    data_subset = abundance_matrix[top_taxa]
    
    # Transform data (log10 with pseudocount)
    data_log = np.log10(data_subset + 1e-6)
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    
    # Add treatment colors if metadata available
    if metadata_available and sample_to_treatment:
        # Create color map for treatments
        treatments = [sample_to_treatment.get(s, 'Unknown') for s in data_subset.index]
        unique_treatments = list(set(treatments))
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_treatments)))
        treatment_colors = [colors[unique_treatments.index(t)] for t in treatments]
        
        # Create clustermap with row colors
        g = sns.clustermap(data_log.T, 
                          col_colors=treatment_colors,
                          cmap='YlOrRd', 
                          figsize=(14, 10),
                          cbar_kws={'label': 'Log10(Abundance)'},
                          xticklabels=True,
                          yticklabels=True)
        
        # Add legend for treatments
        for i, treatment in enumerate(unique_treatments):
            g.ax_col_dendrogram.bar(0, 0, color=colors[i], label=treatment, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=min(5, len(unique_treatments)), 
                                   bbox_to_anchor=(0.5, 0.9), title="Treatment")
    else:
        # Create clustermap without row colors
        g = sns.clustermap(data_log.T, 
                          cmap='YlOrRd', 
                          figsize=(14, 10),
                          cbar_kws={'label': 'Log10(Abundance)'},
                          xticklabels=True,
                          yticklabels=True)
    
    # Customize
    g.ax_heatmap.set_xlabel('Samples', fontsize=12)
    g.ax_heatmap.set_ylabel(f'{level_name}', fontsize=12)
    plt.suptitle(f'{level_name} Abundance Heatmap (Top {top_n})', fontsize=14, fontweight='bold', y=0.98)
    
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_abundance_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_pca_analysis(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available):
    """Perform and plot PCA analysis"""
    print(f"  - Creating PCA plot...")
    
    # Prepare data for PCA
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(abundance_matrix)
    
    # Perform PCA
    pca = PCA(n_components=min(3, len(abundance_matrix)))
    pca_result = pca.fit_transform(data_scaled)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot PC1 vs PC2
    if metadata_available and sample_to_treatment:
        treatments = [sample_to_treatment.get(s, 'Unknown') for s in abundance_matrix.index]
        unique_treatments = list(set(treatments))
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_treatments)))
        
        for treatment in unique_treatments:
            mask = [t == treatment for t in treatments]
            ax1.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                       label=treatment, s=100, alpha=0.7)
        ax1.legend(title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax1.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)
        # Add sample labels
        for i, sample in enumerate(abundance_matrix.index):
            ax1.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                        fontsize=8, alpha=0.6)
    
    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax1.set_title(f'PCA - {level_name} Level')
    ax1.grid(True, alpha=0.3)
    
    # Plot explained variance
    ax2.bar(range(1, len(pca.explained_variance_ratio_) + 1), 
            pca.explained_variance_ratio_ * 100)
    ax2.set_xlabel('Principal Component')
    ax2.set_ylabel('Explained Variance (%)')
    ax2.set_title('Scree Plot')
    ax2.set_xticks(range(1, len(pca.explained_variance_ratio_) + 1))
    
    plt.suptitle(f'PCA Analysis - {level_name} Level', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_pca_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_diversity_metrics(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available):
    """Calculate and plot alpha diversity metrics"""
    print(f"  - Creating diversity metrics plot...")
    
    # Calculate diversity metrics
    # Shannon diversity
    shannon_diversity = abundance_matrix.apply(lambda x: entropy(x[x > 0]), axis=1)
    
    # Simpson diversity
    simpson_diversity = abundance_matrix.apply(lambda x: 1 - sum(x**2), axis=1)
    
    # Observed richness
    observed_richness = (abundance_matrix > 0).sum(axis=1)
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Prepare data for plotting
    if metadata_available and sample_to_treatment:
        df = pd.DataFrame({
            'Shannon': shannon_diversity,
            'Simpson': simpson_diversity,
            'Richness': observed_richness,
            'Treatment': [sample_to_treatment.get(s, 'Unknown') for s in abundance_matrix.index]
        })
        
        # Shannon diversity
        df.boxplot(column='Shannon', by='Treatment', ax=axes[0])
        axes[0].set_title('Shannon Diversity')
        axes[0].set_xlabel('Treatment')
        axes[0].set_ylabel('Shannon Index')
        axes[0].get_figure().suptitle('')
        
        # Simpson diversity
        df.boxplot(column='Simpson', by='Treatment', ax=axes[1])
        axes[1].set_title('Simpson Diversity')
        axes[1].set_xlabel('Treatment')
        axes[1].set_ylabel('Simpson Index')
        
        # Observed richness
        df.boxplot(column='Richness', by='Treatment', ax=axes[2])
        axes[2].set_title('Observed Richness')
        axes[2].set_xlabel('Treatment')
        axes[2].set_ylabel('Number of Taxa')
    else:
        # Create bar plots for each sample
        x_pos = np.arange(len(abundance_matrix.index))
        
        axes[0].bar(x_pos, shannon_diversity)
        axes[0].set_title('Shannon Diversity')
        axes[0].set_xlabel('Sample')
        axes[0].set_ylabel('Shannon Index')
        axes[0].set_xticks(x_pos)
        axes[0].set_xticklabels(abundance_matrix.index, rotation=45, ha='right')
        
        axes[1].bar(x_pos, simpson_diversity)
        axes[1].set_title('Simpson Diversity')
        axes[1].set_xlabel('Sample')
        axes[1].set_ylabel('Simpson Index')
        axes[1].set_xticks(x_pos)
        axes[1].set_xticklabels(abundance_matrix.index, rotation=45, ha='right')
        
        axes[2].bar(x_pos, observed_richness)
        axes[2].set_title('Observed Richness')
        axes[2].set_xlabel('Sample')
        axes[2].set_ylabel('Number of Taxa')
        axes[2].set_xticks(x_pos)
        axes[2].set_xticklabels(abundance_matrix.index, rotation=45, ha='right')
    
    plt.suptitle(f'Alpha Diversity Metrics - {level_name} Level', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_diversity_metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_stacked_composition(abundance_matrix, level_name, output_dir, sample_to_treatment, metadata_available):
    """Create stacked bar plot of taxonomic composition"""
    print(f"  - Creating stacked composition plot...")
    
    # Select top taxa
    top_n = 15
    mean_abundance = abundance_matrix.mean(axis=0).sort_values(ascending=False)
    top_taxa = mean_abundance.head(top_n).index
    
    # Create "Others" category
    data_top = abundance_matrix[top_taxa].copy()
    data_top['Others'] = abundance_matrix.drop(columns=top_taxa).sum(axis=1)
    
    # Sort samples by treatment if available
    if metadata_available and sample_to_treatment:
        sample_order = sorted(abundance_matrix.index, 
                            key=lambda x: sample_to_treatment.get(x, 'Unknown'))
        data_top = data_top.loc[sample_order]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create stacked bar plot
    data_top.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.8)
    
    # Customize
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Relative Abundance', fontsize=12)
    ax.set_title(f'{level_name} Composition (Top {top_n} + Others)', fontsize=14, fontweight='bold')
    ax.legend(title='Taxa', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=1)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Add treatment labels if available
    if metadata_available and sample_to_treatment:
        # Add colored bars under x-axis to show treatments
        treatments = [sample_to_treatment.get(s, 'Unknown') for s in data_top.index]
        unique_treatments = list(set(treatments))
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_treatments)))
        
        for i, sample in enumerate(data_top.index):
            treatment = sample_to_treatment.get(sample, 'Unknown')
            color = colors[unique_treatments.index(treatment)]
            ax.axhspan(-0.05, -0.02, xmin=i/len(data_top.index), 
                      xmax=(i+1)/len(data_top.index), color=color, 
                      transform=ax.get_xaxis_transform())
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_stacked_composition.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_group_comparison_boxplots(abundance_matrix, level_name, output_dir, sample_to_treatment):
    """Create boxplots comparing top taxa between groups"""
    print(f"  - Creating group comparison boxplots...")
    
    # Select top taxa
    top_n = 10
    mean_abundance = abundance_matrix.mean(axis=0).sort_values(ascending=False)
    top_taxa = mean_abundance.head(top_n).index
    
    # Prepare data
    data_list = []
    for taxon in top_taxa:
        for sample in abundance_matrix.index:
            data_list.append({
                'Taxon': taxon,
                'Abundance': abundance_matrix.loc[sample, taxon],
                'Treatment': sample_to_treatment.get(sample, 'Unknown'),
                'Sample': sample
            })
    
    df = pd.DataFrame(data_list)
    
    # Create figure
    fig, axes = plt.subplots(2, 5, figsize=(18, 8))
    axes = axes.flatten()
    
    for i, taxon in enumerate(top_taxa):
        ax = axes[i]
        data_taxon = df[df['Taxon'] == taxon]
        treatments = data_taxon['Treatment'].unique()
        
        # Create boxplot data
        boxplot_data = [data_taxon[data_taxon['Treatment'] == t]['Abundance'].values 
                       for t in treatments]
        
        bp = ax.boxplot(boxplot_data, labels=treatments)
        ax.set_title(taxon[:20] + '...' if len(taxon) > 20 else taxon, fontsize=10)
        ax.set_xlabel('Treatment', fontsize=9)
        ax.set_ylabel('Abundance', fontsize=9)
        ax.tick_params(axis='x', rotation=45)
    
    plt.suptitle(f'Top {top_n} {level_name} Comparison Between Groups', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_group_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_beta_diversity_between_groups(abundance_matrix, level_name, output_dir, sample_to_treatment):
    """Calculate and visualize beta diversity between groups"""
    print(f"  - Creating beta diversity plot...")
    
    # Calculate Bray-Curtis dissimilarity
    bray_curtis = pairwise_distances(abundance_matrix, metric='braycurtis')
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot distance matrix
    im = ax1.imshow(bray_curtis, cmap='YlOrRd', aspect='auto')
    ax1.set_xticks(range(len(abundance_matrix.index)))
    ax1.set_yticks(range(len(abundance_matrix.index)))
    ax1.set_xticklabels(abundance_matrix.index, rotation=45, ha='right')
    ax1.set_yticklabels(abundance_matrix.index)
    ax1.set_title('Bray-Curtis Dissimilarity Matrix')
    plt.colorbar(im, ax=ax1)
    
    # Plot distribution of within vs between group distances
    treatments = [sample_to_treatment.get(s, 'Unknown') for s in abundance_matrix.index]
    unique_treatments = list(set(treatments))
    
    within_distances = []
    between_distances = []
    
    for i in range(len(abundance_matrix.index)):
        for j in range(i+1, len(abundance_matrix.index)):
            if treatments[i] == treatments[j]:
                within_distances.append(bray_curtis[i, j])
            else:
                between_distances.append(bray_curtis[i, j])
    
    # Create violin plot
    if within_distances and between_distances:
        parts = ax2.violinplot([within_distances, between_distances], 
                               positions=[1, 2], widths=0.7,
                               showmeans=True, showmedians=True)
        ax2.set_xticks([1, 2])
        ax2.set_xticklabels(['Within Groups', 'Between Groups'])
        ax2.set_ylabel('Bray-Curtis Dissimilarity')
        ax2.set_title('Beta Diversity Comparison')
        
        # Add statistical test
        if len(within_distances) > 1 and len(between_distances) > 1:
            statistic, pvalue = mannwhitneyu(within_distances, between_distances)
            ax2.text(1.5, ax2.get_ylim()[1] * 0.95, 
                    f'Mann-Whitney U test\np-value: {pvalue:.3e}',
                    ha='center', va='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle(f'Beta Diversity Analysis - {level_name} Level', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{level_name.lower()}_beta_diversity.png'), dpi=300, bbox_inches='tight')
    plt.close()


def create_summary_dashboard(project_output, taxonomic_levels, metadata, metadata_available):
    """Create an interactive HTML dashboard summarizing all results"""
    print("  - Creating interactive summary dashboard...")
    
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Abundance Analysis Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1 { color: #333; }
            h2 { color: #666; margin-top: 30px; }
            .summary-box { 
                background-color: #f0f0f0; 
                padding: 15px; 
                border-radius: 5px; 
                margin: 10px 0;
            }
            .level-section {
                margin: 20px 0;
                padding: 20px;
                border: 1px solid #ddd;
                border-radius: 5px;
            }
            table { border-collapse: collapse; width: 100%; margin: 10px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #4CAF50; color: white; }
            .plot-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; }
            .plot-item { text-align: center; }
            img { max-width: 100%; height: auto; border: 1px solid #ddd; }
        </style>
    </head>
    <body>
        <h1>Metagenomic Abundance Analysis Summary</h1>
    """
    
    # Add analysis overview
    html_content += f"""
        <div class="summary-box">
            <h2>Analysis Overview</h2>
            <p><strong>Project Output:</strong> {project_output}</p>
            <p><strong>Taxonomic Levels Analyzed:</strong> {', '.join(taxonomic_levels)}</p>
            <p><strong>Metadata Available:</strong> {'Yes' if metadata_available else 'No'}</p>
    """
    
    if metadata_available and metadata is not None:
        html_content += f"""
            <p><strong>Number of Samples:</strong> {len(metadata)}</p>
            <p><strong>Treatment Groups:</strong> {', '.join(metadata['Treatment'].unique()) if 'Treatment' in metadata.columns else 'N/A'}</p>
        """
    
    html_content += "</div>"
    
    # Add results for each taxonomic level
    for level in taxonomic_levels:
        level_name = {'S': 'Species', 'G': 'Genus', 'F': 'Family', 
                     'O': 'Order', 'C': 'Class', 'P': 'Phylum'}.get(level, level)
        
        html_content += f"""
        <div class="level-section">
            <h2>{level_name} Level Analysis</h2>
            
            <div class="plot-grid">
                <div class="plot-item">
                    <h3>Top Taxa</h3>
                    <img src="Abundance_Plots/top_20_{level_name.lower()}_barplot.png" alt="Top Taxa">
                </div>
                <div class="plot-item">
                    <h3>Abundance Heatmap</h3>
                    <img src="Abundance_Plots/{level_name.lower()}_abundance_heatmap.png" alt="Heatmap">
                </div>
                <div class="plot-item">
                    <h3>PCA Analysis</h3>
                    <img src="Abundance_Plots/{level_name.lower()}_pca_analysis.png" alt="PCA">
                </div>
                <div class="plot-item">
                    <h3>Diversity Metrics</h3>
                    <img src="Abundance_Plots/{level_name.lower()}_diversity_metrics.png" alt="Diversity">
                </div>
                <div class="plot-item">
                    <h3>Taxonomic Composition</h3>
                    <img src="Abundance_Plots/{level_name.lower()}_stacked_composition.png" alt="Composition">
                </div>
        """
        
        if metadata_available:
            html_content += """
                <div class="plot-item">
                    <h3>Group Comparisons</h3>
                    <img src="Abundance_Plots/{}_group_comparison.png" alt="Group Comparison">
                </div>
                <div class="plot-item">
                    <h3>Beta Diversity</h3>
                    <img src="Abundance_Plots/{}_beta_diversity.png" alt="Beta Diversity">
                </div>
            """.format(level_name.lower(), level_name.lower())
        
        html_content += """
            </div>
        </div>
        """
    
    html_content += """
    </body>
    </html>
    """
    
    # Save HTML file
    output_file = os.path.join(project_output, 'abundance_summary.html')
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"  - Dashboard saved to: {output_file}")


# Main function with force_rerun option
def run(samples, cpus, memory, time, project_input, project_output, kraken_db_path,
        taxonomic_levels=['S', 'G'], read_length=150, threshold=10, metadata_file=None, 
        force_rerun=False):
    """
    Main function for comprehensive abundance estimation with SMART CHECKPOINTING.
    
    Args:
        force_rerun: If True, ignore existing outputs and rerun everything
    """
    print("\n" + "="*70)
    print("COMPREHENSIVE ABUNDANCE ESTIMATION WITH SMART CHECKPOINTING")
    print("="*70)
    
    if force_rerun:
        print("[WARNING] Force rerun enabled - will ignore existing outputs!")
    else:
        print("[INFO] Smart checkpointing enabled - will skip completed steps")
    
    return run_comprehensive_abundance(samples, cpus, memory, time, project_input, 
                                     project_output, kraken_db_path, taxonomic_levels, 
                                     read_length, threshold, metadata_file, force_rerun)


if __name__ == "__main__":
    print("This module is designed to be called from the MetaG pipeline.")
    print("Use run_step.py with --step abundance_estimation")
    print("\nTo force rerun all steps (ignore existing outputs):")
    print("  Add --force_rerun flag to your command")
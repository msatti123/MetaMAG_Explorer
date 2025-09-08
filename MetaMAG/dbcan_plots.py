#!/usr/bin/env python3
"""
Module for creating visualizations of CAZyme data from dbCAN results.
This module provides functions to generate various plots for CAZyme analysis.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

def create_heatmap(cazyme_dir, output_dir=None, selected_families=None, 
                  min_prevalence=0.2, cmap="YlGnBu", figsize=(20, 12)):
    """
    Create a heatmap of CAZyme family abundance across MAGs.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    selected_families : list, optional
        List of CAZyme families to include (e.g., ['GH5', 'GH10']). If None, selects based on prevalence.
    min_prevalence : float, optional
        Minimum fraction of MAGs that must contain a family to be included (0.0-1.0)
    cmap : str, optional
        Colormap for the heatmap
    figsize : tuple, optional
        Figure size (width, height) in inches
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load GH family data
    gh_file = os.path.join(cazyme_dir, 'gh_cazyme_counts.xlsx')
    if not os.path.exists(gh_file):
        print(f"[ERROR] GH CAZyme counts file not found: {gh_file}")
        return
    
    df = pd.read_excel(gh_file, index_col=0)
    
    # Filter data based on prevalence if no specific families selected
    if selected_families is None:
        # Calculate prevalence (fraction of MAGs containing each family)
        prevalence = (df > 0).mean(axis=1)
        selected_families = prevalence[prevalence >= min_prevalence].index.tolist()
    
    # Filter dataframe to selected families
    if selected_families:
        # Ensure all requested families exist in the dataframe
        valid_families = [fam for fam in selected_families if fam in df.index]
        if len(valid_families) < len(selected_families):
            missing = set(selected_families) - set(valid_families)
            print(f"[WARNING] Some requested families were not found in the data: {missing}")
        
        if not valid_families:
            print(f"[ERROR] None of the requested families were found in the data")
            return
            
        df_selected = df.loc[valid_families]
    else:
        print(f"[WARNING] No families meet the minimum prevalence threshold of {min_prevalence}")
        # Use the top 20 families by average abundance
        avg_abundance = df.mean(axis=1)
        df_selected = df.loc[avg_abundance.nlargest(20).index]
    
    # Save the selected data to Excel
    selected_file = os.path.join(output_dir, 'selected_gh.xlsx')
    df_selected.to_excel(selected_file)
    
    # Create heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(df_selected, annot=False, cmap=cmap)
    
    plt.title('Heatmap of GH Families Across MAGs', fontsize=14)
    plt.xlabel('MAGs', fontsize=12)
    plt.ylabel('GH Family', fontsize=12)
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'gh_families_heatmap.pdf'))
    plt.savefig(os.path.join(output_dir, 'gh_families_heatmap.png'), dpi=300)
    plt.close()
    
    print(f"[INFO] Heatmap created and saved to {output_dir}")
    
    # Optional: Create a heatmap with values shown
    if len(df_selected) <= 40 and df_selected.shape[1] <= 40:  # Limit for readability
        plt.figure(figsize=figsize)
        sns.heatmap(df_selected, annot=True, cmap=cmap, fmt="d")
        
        plt.title('Heatmap of GH Families Across MAGs (with values)', fontsize=14)
        plt.xlabel('MAGs', fontsize=12)
        plt.ylabel('GH Family', fontsize=12)
        
        # Adjust layout for better readability
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(os.path.join(output_dir, 'gh_families_heatmap_with_values.pdf'))
        plt.savefig(os.path.join(output_dir, 'gh_families_heatmap_with_values.png'), dpi=300)
        plt.close()

def create_category_barplot(cazyme_dir, output_dir=None, figsize=(12, 8)):
    """
    Create a barplot showing CAZyme category distribution across MAGs.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    figsize : tuple, optional
        Figure size (width, height) in inches
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load category data
    category_file = os.path.join(cazyme_dir, 'category_counts.xlsx')
    if not os.path.exists(category_file):
        print(f"[ERROR] Category counts file not found: {category_file}")
        return
    
    df = pd.read_excel(category_file, index_col=0)
    
    # Remove Total column if it exists
    if 'Total' in df.columns:
        df = df.drop('Total', axis=1)
    
    # Calculate the average count for each category
    category_avg = df.mean()
    
    # Create stacked bar chart of all MAGs
    plt.figure(figsize=figsize)
    df.plot(kind='bar', stacked=True, figsize=figsize, colormap='viridis')
    
    plt.title('CAZyme Category Distribution by MAG', fontsize=14)
    plt.xlabel('MAG', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend(title='CAZyme Category')
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'cazyme_category_distribution.pdf'))
    plt.savefig(os.path.join(output_dir, 'cazyme_category_distribution.png'), dpi=300)
    plt.close()
    
    # Create average category bar chart
    plt.figure(figsize=(10, 6))
    category_avg.plot(kind='bar', figsize=(10, 6), color=sns.color_palette('viridis', len(category_avg)))
    
    plt.title('Average CAZyme Category Counts Across MAGs', fontsize=14)
    plt.xlabel('CAZyme Category', fontsize=12)
    plt.ylabel('Average Count', fontsize=12)
    
    # Add count values on top of bars
    for i, v in enumerate(category_avg):
        plt.text(i, v + 0.5, f'{v:.1f}', ha='center')
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'average_cazyme_categories.pdf'))
    plt.savefig(os.path.join(output_dir, 'average_cazyme_categories.png'), dpi=300)
    plt.close()
    
    print(f"[INFO] Category distribution plots created and saved to {output_dir}")

def create_gh_family_distribution(cazyme_dir, output_dir=None, top_n=20, figsize=(14, 10)):
    """
    Create plots showing the distribution of GH (Glycoside Hydrolase) families.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    top_n : int, optional
        Number of top GH families to include in the plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load GH family data
    gh_file = os.path.join(cazyme_dir, 'gh_cazyme_counts.xlsx')
    if not os.path.exists(gh_file):
        print(f"[ERROR] GH CAZyme counts file not found: {gh_file}")
        return
    
    df = pd.read_excel(gh_file, index_col=0)
    
    # Calculate total count for each GH family across all MAGs
    gh_totals = df.sum(axis=1).sort_values(ascending=False)
    
    # Get the top N GH families
    top_gh_families = gh_totals.head(top_n)
    
    # Create bar chart of top GH families
    plt.figure(figsize=figsize)
    top_gh_families.plot(kind='bar', figsize=figsize, color=sns.color_palette('viridis', len(top_gh_families)))
    
    plt.title(f'Top {top_n} GH Families Across All MAGs', fontsize=14)
    plt.xlabel('GH Family', fontsize=12)
    plt.ylabel('Total Count', fontsize=12)
    
    # Add count values on top of bars
    for i, v in enumerate(top_gh_families):
        plt.text(i, v + 0.5, str(v), ha='center')
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f'top_{top_n}_gh_families.pdf'))
    plt.savefig(os.path.join(output_dir, f'top_{top_n}_gh_families.png'), dpi=300)
    plt.close()
    
    # Calculate prevalence (fraction of MAGs containing each family)
    prevalence = (df > 0).mean(axis=1).sort_values(ascending=False)
    
    # Get the top N GH families by prevalence
    top_prevalent_gh = prevalence.head(top_n)
    
    # Create bar chart of top prevalent GH families
    plt.figure(figsize=figsize)
    top_prevalent_gh.plot(kind='bar', figsize=figsize, color=sns.color_palette('viridis', len(top_prevalent_gh)))
    
    plt.title(f'Top {top_n} Most Prevalent GH Families', fontsize=14)
    plt.xlabel('GH Family', fontsize=12)
    plt.ylabel('Fraction of MAGs', fontsize=12)
    
    # Add prevalence values on top of bars
    for i, v in enumerate(top_prevalent_gh):
        plt.text(i, v + 0.02, f'{v:.2f}', ha='center')
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f'top_{top_n}_prevalent_gh_families.pdf'))
    plt.savefig(os.path.join(output_dir, f'top_{top_n}_prevalent_gh_families.png'), dpi=300)
    plt.close()
    
    print(f"[INFO] GH family distribution plots created and saved to {output_dir}")

def create_cazyme_dendrogram(cazyme_dir, output_dir=None, method='ward', metric='euclidean', figsize=(18, 10)):
    """
    Create a dendrogram clustering MAGs based on their CAZyme profiles.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    method : str, optional
        Linkage method for hierarchical clustering
    metric : str, optional
        Distance metric for hierarchical clustering
    figsize : tuple, optional
        Figure size (width, height) in inches
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load all CAZyme data
    cazyme_file = os.path.join(cazyme_dir, 'all_cazyme_counts.xlsx')
    if not os.path.exists(cazyme_file):
        print(f"[ERROR] All CAZyme counts file not found: {cazyme_file}")
        return
    
    df = pd.read_excel(cazyme_file, index_col=0)
    
    # Transpose so MAGs are rows and CAZymes are columns
    df_t = df.T
    
    # Create a clustergram/heatmap with dendrogram
    plt.figure(figsize=figsize)
    
    # Using clustermap from seaborn which includes dendrogram
    g = sns.clustermap(df_t, method=method, metric=metric, 
                      cmap="YlGnBu", figsize=figsize, 
                      xticklabels=True, yticklabels=True)
    
    # Adjust the plot
    plt.title('Hierarchical Clustering of MAGs by CAZyme Profile', fontsize=14)
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'cazyme_dendrogram.pdf'))
    plt.savefig(os.path.join(output_dir, 'cazyme_dendrogram.png'), dpi=300)
    plt.close()
    
    print(f"[INFO] CAZyme dendrogram created and saved to {output_dir}")

def create_functional_enrichment_plot(cazyme_dir, output_dir=None, figsize=(14, 12)):
    """
    Create plot showing functional enrichment based on CAZyme family types.
    Group GH families by their known substrate preferences.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    figsize : tuple, optional
        Figure size (width, height) in inches
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load GH family data
    gh_file = os.path.join(cazyme_dir, 'gh_cazyme_counts.xlsx')
    if not os.path.exists(gh_file):
        print(f"[ERROR] GH CAZyme counts file not found: {gh_file}")
        return
    
    df = pd.read_excel(gh_file, index_col=0)
    
    # Define common substrate categories for GH families 
    # This is a simplified mapping - in reality, many GH families are more versatile
    substrate_mapping = {
        'Cellulose': ['GH1', 'GH3', 'GH5', 'GH6', 'GH7', 'GH8', 'GH9', 'GH12', 'GH44', 'GH45', 'GH48', 'GH74'],
        'Hemicellulose': ['GH10', 'GH11', 'GH26', 'GH30', 'GH31', 'GH43', 'GH51', 'GH54', 'GH62', 'GH67', 'GH115'],
        'Starch': ['GH13', 'GH14', 'GH15', 'GH57', 'GH119'],
        'Pectin': ['GH28', 'GH53', 'GH78', 'GH88', 'GH93', 'GH105'],
        'Chitin': ['GH18', 'GH19', 'GH20', 'GH85'],
        'Peptidoglycan': ['GH23', 'GH24', 'GH25', 'GH73'],
        'Oligosaccharides': ['GH2', 'GH29', 'GH35', 'GH38', 'GH39', 'GH42'],
        'Other': []  # Will catch all families not in the above categories
    }
    
    # Group GH families by substrate
    substrate_counts = defaultdict(lambda: defaultdict(int))
    
    for gh_family in df.index:
        # Get base family name (without subfamily)
        base_family = gh_family.split('_')[0]
        
        # Find which substrate category this family belongs to
        assigned = False
        for substrate, families in substrate_mapping.items():
            if base_family in families:
                for mag in df.columns:
                    substrate_counts[substrate][mag] += df.loc[gh_family, mag]
                assigned = True
                break
        
        # If not assigned to any specific category, put in "Other"
        if not assigned:
            for mag in df.columns:
                substrate_counts['Other'][mag] += df.loc[gh_family, mag]
    
    # Convert to DataFrame
    substrate_df = pd.DataFrame(substrate_counts)
    
    # Create a stacked bar chart
    plt.figure(figsize=figsize)
    substrate_df.plot(kind='bar', stacked=True, figsize=figsize, colormap='tab10')
    
    plt.title('Substrate-Related CAZyme Distribution by MAG', fontsize=14)
    plt.xlabel('MAG', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend(title='Substrate', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust layout for better readability
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'substrate_cazyme_distribution.pdf'))
    plt.savefig(os.path.join(output_dir, 'substrate_cazyme_distribution.png'), dpi=300)
    plt.close()
    
    # Create a pie chart of overall substrate distribution
    substrate_totals = substrate_df.sum()
    
    plt.figure(figsize=(10, 10))
    plt.pie(substrate_totals, labels=substrate_totals.index, autopct='%1.1f%%', 
            shadow=True, startangle=90, colors=plt.cm.tab10.colors[:len(substrate_totals)])
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    plt.title('Overall Distribution of Substrate-Related CAZymes', fontsize=14)
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, 'substrate_cazyme_pie.pdf'))
    plt.savefig(os.path.join(output_dir, 'substrate_cazyme_pie.png'), dpi=300)
    plt.close()
    
    print(f"[INFO] Functional enrichment plots created and saved to {output_dir}")

def run_all_visualizations(cazyme_dir, output_dir=None):
    """
    Run all visualization functions to create a comprehensive set of plots.
    
    Parameters:
    -----------
    cazyme_dir : str
        Path to directory containing CAZyme annotation results (with Excel files)
    output_dir : str, optional
        Directory to save output plots (defaults to cazyme_dir/plots)
    """
    if output_dir is None:
        output_dir = os.path.join(cazyme_dir, "plots")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"[INFO] Creating all CAZyme visualization plots")
    
    # Run all visualization functions
    create_heatmap(cazyme_dir, output_dir)
    create_category_barplot(cazyme_dir, output_dir)
    create_gh_family_distribution(cazyme_dir, output_dir)
    create_cazyme_dendrogram(cazyme_dir, output_dir)
    create_functional_enrichment_plot(cazyme_dir, output_dir)
    
    print(f"[INFO] All visualization plots created and saved to {output_dir}")
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from MetaMAG.utils import ensure_directory_exists, run_command
from MetaMAG.config import config, get_tool_path, get_tool_command

def run_checkm(samples, cpus, memory, time, project_input, project_output, input_dir=None):
    """
    Runs CheckM lineage workflow on the dereplicated genomes produced by dRep.
    Skips execution if expected outputs already exist.
    
    Parameters:
    -----------
    samples : list
        List of sample IDs (may not be used in this step but needed for interface consistency)
    cpus : int
        Number of CPUs to use
    memory : str
        Memory allocation
    time : str
        Time allocation
    project_input : str
        Base input directory for the project
    project_output : str
        Base output directory for the project
    input_dir : str, optional
        Custom input directory containing genomes for CheckM analysis.
        If not provided, the default dRep output directory will be used.
    """
    print(f"[INFO] Running CheckM on dRep dereplicated genomes")

    base_dir = project_output
    
    # Use provided input_dir if available, otherwise use default path
    if input_dir and os.path.exists(input_dir):
        input_genomes_dir = input_dir
        print(f"[INFO] Using custom input directory: {input_dir}")
    else:
        input_genomes_dir = os.path.join(base_dir, "Bin_Refinement", "drep", "dRep_output", "dereplicated_genomes")
        print(f"[INFO] Using default input directory: {input_genomes_dir}")

    if not os.path.exists(input_genomes_dir) or not os.listdir(input_genomes_dir):
        print(f"[ERROR] Genomes directory is missing or empty: {input_genomes_dir}")
        return

    # Get CheckM tool path
    checkm_tool = get_tool_path("checkm")
    if not checkm_tool:
        print(f"[ERROR] CheckM tool not found in configuration")
        return

    checkm_output = os.path.join(base_dir, "Bin_Refinement", "checkm", "checkm_output")
    ensure_directory_exists(checkm_output)

    # Define key output files
    bins_folder = os.path.join(checkm_output, "bins")
    storage_folder = os.path.join(checkm_output, "storage")
    bin_stats_file = os.path.join(storage_folder, "bin_stats_ext.tsv")
    checkm_csv_file = os.path.join(checkm_output, "cleaned_checkm_output.csv")
    lineage_file = os.path.join(checkm_output, "lineage.ms")
    plots_folder = os.path.join(checkm_output, "plots")
    plot1 = os.path.join(plots_folder, "completeness_contamination_plot.pdf")
    plot2 = os.path.join(plots_folder, "genome_quality_plot.pdf")
    plot3 = os.path.join(plots_folder, "scaffolds_histogram.pdf")

    # === Skip entire run if all expected outputs already exist ===
    if all([
        os.path.exists(bins_folder),
        os.path.exists(storage_folder),
        os.path.exists(lineage_file),
        os.path.exists(bin_stats_file),
        os.path.exists(checkm_csv_file),
        os.path.exists(plot1),
        os.path.exists(plot2),
        os.path.exists(plot3)
    ]):
        print(f"[INFO] All CheckM outputs and plots already exist. Skipping CheckM, processing, and plotting.")
        return

    # === Run CheckM if needed ===
    if not (os.path.exists(bins_folder) and os.path.exists(storage_folder) and os.path.exists(lineage_file)):
        print(f"[INFO] Running CheckM lineage workflow...")
        print(f"       Tool config: {checkm_tool}")
        print(f"       Input: {input_genomes_dir}")
        print(f"       Output: {checkm_output}")
        print(f"       CPUs: {cpus}")
        
        log_file = os.path.join(checkm_output, "checkm_output.log")
        
        # Check if the tool path looks like an activation command (old config style)
        if "activate" in checkm_tool and ("checkm_env" in checkm_tool or "checkm" in checkm_tool):
            print(f"[INFO] Detected conda activation command in config")
            
            # Extract the activation script path and environment name
            # Typically format: "/path/to/bin/activate env_name"
            parts = checkm_tool.split()
            if len(parts) >= 2:
                activate_script = parts[0]
                env_name = parts[1]
            else:
                # Fallback - try to extract from environment config
                conda_env_cmd = config.get("environment", {}).get("conda_env", "")
                if "source" in conda_env_cmd and "activate" in conda_env_cmd:
                    parts = conda_env_cmd.split()
                    for i, part in enumerate(parts):
                        if "activate" in part and i > 0:
                            activate_script = parts[i]
                            break
                    else:
                        activate_script = "conda activate"
                else:
                    activate_script = "conda activate"
                env_name = "checkm_env"
            
            # Construct command with conda activation
            cmd = (
                f"bash -c 'source {activate_script} {env_name} && "
                f"checkm lineage_wf -x fa {input_genomes_dir} {checkm_output} --threads {cpus}'"
            )
        else:
            # Direct path to checkm executable
            print(f"[INFO] Using direct path to CheckM executable")
            cmd = f"{checkm_tool} lineage_wf -x fa {input_genomes_dir} {checkm_output} --threads {cpus}"
        
        print(f"[DEBUG] Running command: {cmd}")
        
        try:
            # Redirect output to log file
            cmd_with_log = f"{cmd} > {log_file} 2>&1"
            run_command(cmd_with_log)
            print(f"[INFO] CheckM lineage workflow completed.")
            print(f"[INFO] CheckM log saved to: {log_file}")
        except Exception as e:
            print(f"[ERROR] CheckM failed: {e}")
            if os.path.exists(log_file):
                print(f"[INFO] Check CheckM log for details: {log_file}")
                # Try to print last few lines of the log for debugging
                try:
                    with open(log_file, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            print("[DEBUG] Last few lines of CheckM log:")
                            for line in lines[-10:]:
                                print(f"       {line.strip()}")
                except:
                    pass
            return

    # === Generate CSV if needed ===
    if not os.path.exists(checkm_csv_file):
        print(f"[INFO] Processing CheckM bin stats...")
        process_bin_stats(checkm_output)

    # === Generate plots if needed ===
    if not (os.path.exists(plot1) and os.path.exists(plot2) and os.path.exists(plot3)):
        print(f"[INFO] Generating CheckM plots...")
        generate_plots(checkm_output)

    print(f"[INFO] CheckM evaluation completed successfully!")

def process_bin_stats(checkm_output):
    """
    Extracts CheckM data from `bin_stats_ext.tsv` and saves it as a CSV.
    """
    bin_stats_file = os.path.join(checkm_output, 'storage', 'bin_stats_ext.tsv')

    # Check if the file exists
    if not os.path.exists(bin_stats_file):
        print(f"[ERROR] The file bin_stats_ext.tsv was not found at: {bin_stats_file}")
        return

    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(bin_stats_file, sep='\t', header=None, names=["Bin ID", "Data"])

    # Process the bin data
    data_cleaned = []
    for _, row in df.iterrows():
        bin_id = row["Bin ID"]
        try:
            bin_data = eval(row["Data"])  # Convert string to dictionary
            data_cleaned.append([
                bin_id,
                bin_data.get('marker lineage', 'N/A'),
                bin_data.get('Completeness', 'N/A'),
                bin_data.get('Contamination', 'N/A'),
                bin_data.get('GC', 'N/A'),
                bin_data.get('Genome size', 'N/A'),
                bin_data.get('# scaffolds', 'N/A'),
                bin_data.get('# contigs', 'N/A'),
                bin_data.get('Longest scaffold', 'N/A')
            ])
        except Exception as e:
            print(f"[ERROR] Failed to parse bin data for {bin_id}: {e}")

    # Define columns for the CSV
    columns = ['Bin ID', 'Marker Lineage', 'Completeness', 'Contamination', 'GC',
               'Genome Size', '# Scaffolds', '# Contigs', 'Longest Scaffold']
    
    # Create DataFrame from cleaned data
    df_cleaned = pd.DataFrame(data_cleaned, columns=columns)

    # Save the DataFrame as a CSV file
    csv_file = os.path.join(checkm_output, "cleaned_checkm_output.csv")
    df_cleaned.to_csv(csv_file, index=False)
    print(f"[INFO] Cleaned CheckM table saved as CSV: {csv_file}")

def generate_plots(checkm_output):
    """
    Generates three plots based on CheckM output CSV.
    """
    print("[INFO] Generating plots based on CheckM CSV output...")
    
    csv_file = os.path.join(checkm_output, "cleaned_checkm_output.csv")
    
    if not os.path.exists(csv_file):
        print(f"[ERROR] CheckM CSV file not found: {csv_file}")
        return
        
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"[ERROR] Failed to read CheckM CSV: {e}")
        return

    # Clean column names
    df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_").str.replace("#", "")
    print(f"[DEBUG] Processed Column Names: {df.columns.tolist()}")

    # Create a dedicated folder for plots
    plots_folder = os.path.join(checkm_output, "plots")
    ensure_directory_exists(plots_folder)

    # **Generate all three plots**
    plot1 = os.path.join(plots_folder, "completeness_contamination_plot.pdf")
    plot2 = os.path.join(plots_folder, "genome_quality_plot.pdf")
    plot3 = os.path.join(plots_folder, "scaffolds_histogram.pdf")

    try:
        # Generate first plot: Completeness vs Contamination (Bubble Plot)
        generate_completeness_contamination_plot(df, plot1)
        
        # Generate second plot: Quality Categorization
        generate_quality_plot(df, plot2)
        
        # Generate third plot: Histogram of Scaffolds
        generate_scaffolds_histogram(df, plot3)
        
        print("[INFO] All plots generated successfully!")
        
    except Exception as e:
        print(f"[ERROR] Failed to generate plots: {e}")

def generate_completeness_contamination_plot(df, output_path):
    """Generate completeness vs contamination bubble plot."""
    # Convert genome size into MB and categorize into 1-5 MB
    df["genome_size_mb"] = (df["genome_size"] / 1e6).round().astype('Int64')
    df["genome_size_mb"] = df["genome_size_mb"].clip(1, 5)  # Ensure values range from 1MB to 5MB
    
    # Define a high-contrast and vibrant color palette
    high_contrast_palette = {
        1: "#FF0000",  # Pure Red
        2: "#FF7F00",  # Bright Orange
        3: "#00FF00",  # Neon Green
        4: "#0000FF",  # Pure Blue
        5: "#FF00FF"   # Bright Magenta
    }

    # Create a scatter plot with circle markers and vibrant colors
    plt.figure(figsize=(10, 8))

    scatter = sns.scatterplot(
        data=df,
        x="completeness",
        y="contamination",
        size="genome_size_mb",
        sizes=(20, 200),  # Adjust size range for better visualization
        hue="genome_size_mb",
        palette=high_contrast_palette,  # Apply the vibrant color scheme
        alpha=0.8,
        edgecolor="k",
        marker="o"  # Use circle markers
    )

    # Labels and title
    plt.xlabel("Completeness (%)", fontsize=12, fontweight='bold')
    plt.ylabel("Contamination (%)", fontsize=12, fontweight='bold')
    plt.title("Distribution of Completeness and Contamination\n(Bubble Size = Genome Size in MB)", fontsize=14, fontweight='bold')

    # Ensure the legend includes all sizes from 1 to 5 MB with corresponding colors
    unique_sizes = [1, 2, 3, 4, 5]  # Ensure all sizes are present in the legend
    legend_handles = [plt.scatter([], [], s=size * 50, color=high_contrast_palette[size], edgecolor="k", alpha=0.8, marker="o") 
                      for size in unique_sizes]

    # Update legend with correctly labeled genome sizes
    plt.legend(legend_handles, [f"{size} MB" for size in unique_sizes], 
               title="Genome Size (MB)", bbox_to_anchor=(1.05, 1), loc="upper left")

    # Grid for readability
    plt.grid(True, linestyle="--", alpha=0.6)

    # Save the figure as a PDF
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()  # Close to free memory

    print(f"Plot saved successfully as {output_path}")

def generate_quality_plot(df, output_path):
    """Generate genome quality categorization plot."""
    # Ensure Completeness and Contamination are numeric
    df["completeness"] = pd.to_numeric(df["completeness"], errors="coerce")
    df["contamination"] = pd.to_numeric(df["contamination"], errors="coerce")

    # Drop any missing values in required columns
    df = df.dropna(subset=["completeness", "contamination"])

    # Categorize genomes based on Completeness and Contamination
    def categorize_quality(row):
        if row["completeness"] >= 90 and row["contamination"] <= 5:
            return "Near Complete"
        elif row["completeness"] >= 70 and row["contamination"] <= 10:
            return "Medium Quality"
        elif row["completeness"] >= 50 and row["contamination"] <= 4:
            return "Partial"
        else:
            return "Low Quality"

    df["quality"] = df.apply(categorize_quality, axis=1)

    # Remove "Low Quality" entries for a cleaner visualization
    df = df[df["quality"].isin(["Near Complete", "Medium Quality", "Partial"])]

    # Count genomes in each category
    total_genomes = len(df)
    quality_counts = df["quality"].value_counts()
    quality_percentages = (quality_counts / total_genomes * 100).round(1)

    # Define updated legend labels with genome counts
    updated_legend_labels = {
        category: f"{category} ({quality_counts[category]} genomes, {quality_percentages[category]}%)"
        for category in quality_counts.index
    }

    # Define color palette
    quality_palette = {
        "Near Complete": "#D73027",  # Red
        "Medium Quality": "#4575B4",  # Blue
        "Partial": "#4D4D4D"  # Grey
    }

    # Create figure with subplots for scatter plot and stacked histograms
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(4, 4, hspace=0.1, wspace=0.1)

    # Scatter plot for Completeness vs Contamination with a boxed frame
    ax_scatter = fig.add_subplot(gs[1:, :-1])
    sns.scatterplot(
        data=df,
        x="completeness",
        y="contamination",
        hue="quality",
        palette=quality_palette,
        alpha=0.7,
        edgecolor="k",
        marker="o",
        ax=ax_scatter
    )

    # Set x and y labels (no title)
    ax_scatter.set_xlabel("Completeness (%)", fontsize=12, fontweight="bold")
    ax_scatter.set_ylabel("Contamination (%)", fontsize=12, fontweight="bold")

    # Add a box around the scatter plot
    for spine in ax_scatter.spines.values():
        spine.set_visible(True)

    # Stacked histogram for Completeness on top (without title)
    ax_histx = fig.add_subplot(gs[0, :-1], sharex=ax_scatter)
    sns.histplot(data=df, x="completeness", hue="quality", palette=quality_palette, bins=30, multiple="stack", alpha=0.7, ax=ax_histx, legend=False)
    ax_histx.axis("off")  # Remove labels and axes

    # Stacked histogram for Contamination on the right (without title)
    ax_histy = fig.add_subplot(gs[1:, -1], sharey=ax_scatter)
    sns.histplot(data=df, y="contamination", hue="quality", palette=quality_palette, bins=30, multiple="stack", alpha=0.7, ax=ax_histy, legend=False)
    ax_histy.axis("off")  # Remove labels and axes

    # Remove duplicate legend entries and adjust position
    handles, labels = ax_scatter.get_legend_handles_labels()
    new_handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=10, markerfacecolor=quality_palette[label],
                               markeredgecolor="black") for label in labels if label in quality_palette]
    new_labels = [updated_legend_labels[label] for label in labels if label in updated_legend_labels]

    # Place legend outside the plot with proper circle markers and black outline
    ax_scatter.legend(
        handles=new_handles,
        labels=new_labels,
        title="Genome Quality",
        bbox_to_anchor=(1.35, 1),
        loc="upper left",
        frameon=False
    )

    # Save the figure as a PDF
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()  # Close to free memory

    print(f"Plot saved successfully as {output_path}")

def generate_scaffolds_histogram(df, output_path):
    """Generate histogram of scaffolds per genome."""
    # Ensure required columns are numeric
    df["completeness"] = pd.to_numeric(df["completeness"], errors="coerce")
    df["contamination"] = pd.to_numeric(df["contamination"], errors="coerce")
    df["scaffolds"] = pd.to_numeric(df["_scaffolds"], errors="coerce")

    # Drop any missing values in required columns
    df = df.dropna(subset=["completeness", "contamination", "scaffolds"])

    # Categorize genomes based on Completeness and Contamination
    def categorize_quality(row):
        if row["completeness"] >= 90 and row["contamination"] <= 5:
            return "Near Complete"
        elif row["completeness"] >= 70 and row["contamination"] <= 10:
            return "Medium Quality"
        elif row["completeness"] >= 50 and row["contamination"] <= 4:
            return "Partial"
        else:
            return "Low Quality"

    df["quality"] = df.apply(categorize_quality, axis=1)

    # Remove "Low Quality" entries for a cleaner visualization
    df = df[df["quality"].isin(["Near Complete", "Medium Quality", "Partial"])]

    # Define original color palette
    original_palette = {
        "Near Complete": "#D73027",  # Red
        "Medium Quality": "#4575B4",  # Blue
        "Partial": "#4D4D4D"  # Grey
    }    

    # Create histogram for the number of scaffolds stacked by genome quality
    plt.figure(figsize=(8, 6))
    histplot = sns.histplot(
        data=df,
        x="scaffolds",
        hue="quality",
        palette=original_palette,
        bins=30,
        multiple="stack",
        alpha=0.8
    )

    # Set labels and title
    plt.xlabel("Number of Scaffolds", fontsize=12, fontweight="bold")
    plt.ylabel("Genomes (%)", fontsize=12, fontweight="bold")
    plt.title("Number of Scaffolds per Genome", fontsize=14, fontweight="bold")

    # Create a custom legend with proper color labels
    legend_handles = [
        mpatches.Patch(color=original_palette["Near Complete"], label="Near Complete"),
        mpatches.Patch(color=original_palette["Medium Quality"], label="Medium Quality"),
        mpatches.Patch(color=original_palette["Partial"], label="Partial")
    ]

    plt.legend(handles=legend_handles, title="Genome Quality", bbox_to_anchor=(1.05, 1), loc="upper left")

    # Remove grid lines
    plt.grid(False)

    # Save the figure as a PDF
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()  # Close to free memory
    
    print(f"Plot saved successfully as {output_path}")

def run(samples, cpus, memory, time, project_input, project_output, evaluation_input_dir=None):
    """
    Main function to run CheckM evaluation.
    This is called by the pipeline step runner.
    """
    print(f"[INFO] Starting CheckM evaluation")
    run_checkm(samples, cpus, memory, time, project_input, project_output, evaluation_input_dir)
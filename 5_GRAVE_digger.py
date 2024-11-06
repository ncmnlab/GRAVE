import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import argparse
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze and visualize gene abundance data.")
    parser.add_argument('--output_dir', type=str, required=True, help='Directory for combined abundance output.')
    parser.add_argument('--metadata_file', type=str, required=True, help='CSV file containing metadata.')
    parser.add_argument('--compiled_abundance_file', type=str, required=True, help='TSV file with compiled gene abundance.')
    parser.add_argument('--type', type=str, required=True, help='Class/Mechanism')
    parser.add_argument('--colour', type=str, default="tab20", help='Colour Gradient combination options of matplotlib such as viridis ; plasma ; inferno ; magma ; cividis ; Blues ; Greens ; Oranges ; Purples ; Reds ; Greys ; YlOrBr (Yellow-Orange-Brown) ; YlOrRd (Yellow-Orange-Red) ; PuBu (Purple-Blue) ; BuGn (Blue-Green) ; BuPu (Blue-Purple) ; GnBu (Green-Blue) ; OrRd (Orange-Red) ; PuRd (Purple-Red) ; tab10 ; tab20 ; tab20b ; tab20c ; Pastel1 ; Pastel2 ; Paired ; Set1 ; Set2 ; Set3 ; Dark2 ; Accent    etc., (optional, default=tab20).')
    
    return parser.parse_args()

def generate_percentage_colormap(abundance_file, metadata_file, output_dir, category, legend_title, colour):
    # Read TSV
    data = pd.read_csv(abundance_file, sep='\t', index_col=0)

    # Read the metadata file (CSV format)
    metadata = pd.read_csv(metadata_file, index_col=0)

    # metadata's sample_ID matches the abundance data
    metadata = metadata.loc[data.columns]

    # Normalize the abundance data by columns (samples)
    data_normalized = data.div(data.sum(axis=0), axis=1)

    # Grouping the data
    metadata_grouped = metadata.groupby(category).groups

    # Define a stacked bar color
    qualitative_colors = sns.color_palette(colour, len(data_normalized.index))

    # Create a dictionary to map genes to colors
    gene_color_map = dict(zip(data_normalized.index, qualitative_colors))

    # Create a figure and axis for the bar chart
    fig, ax = plt.subplots(figsize=(10, 8))

    # Prepare log file path
    log_file_path = os.path.join(output_dir, f"{category}_percentage.log")
    
    # Create a DataFrame to store percentage data for logging
    log_data = pd.DataFrame(index=data_normalized.index)

    # Plot each group of samples (gene percentage)
    for group, sample_ids in metadata_grouped.items():
        # Sum normalized abundance for each gene in the current group
        group_abundance = data_normalized[sample_ids].sum(axis=1)
        
        # Convert abundance to percentages so each group sums to 100
        group_abundance_percentage = (group_abundance / group_abundance.sum()) * 100
        
        # Add percentage data to the log
        log_data[group] = group_abundance_percentage
        
        bottom = 0  # Initialize bottom for stacking

        # Plot each gene's abundance percentage in the current group
        for gene, color in gene_color_map.items():
            ax.bar(group, group_abundance_percentage[gene], bottom=bottom, color=color)
            bottom += group_abundance_percentage[gene]

    # Save the percentage data to a log file
    log_data.to_csv(log_file_path)
    print(f"Percentage data table saved to {log_file_path}")

    # Customize the chart
    ax.set_title(f'Stacked Bar Plot Grouped by {category}', fontsize=20)
    ax.set_xlabel(category, fontsize=15)
    ax.set_ylabel('Composition Percentage (%)', fontsize=15)

    # Rotate x-axis labels if needed
    plt.xticks(rotation=90)

    # Creating a legend with gene names and corresponding colors
    legend_elements = [Patch(facecolor=color, label=gene) for gene, color in gene_color_map.items()]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, title=legend_title)

    # Save the plot to a file
    output_path = os.path.join(output_dir, f"stacked_bar_{category}.pdf")
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0.1)

    # Display the plot
    plt.close()
    print(f"Stacked bar plot grouped by {category} saved to {output_path}")

def main():
    # Print welcome message
    print("[-------------------------------------------------------------------------------------] Namaskaram !! Welcome here as 'Grave filler' and Get ready to see the Great Bar of Resistants {{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}")
    time.sleep(10)
    
    # Parse arguments
    args = parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Read the metadata to get categories for qualitative colormap
    metadata = pd.read_csv(args.metadata_file, index_col=0)

    # Generate percentage colormap for each category in the metadata
    for category in metadata.columns:
        print(f"Generating percentage colormap for {category}...")
        generate_percentage_colormap(
            abundance_file=args.compiled_abundance_file,
            metadata_file=args.metadata_file,
            output_dir=args.output_dir,
            category=category,
            legend_title=args.type,
            colour=args.colour
        )

if __name__ == "__main__":
    main()


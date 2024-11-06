import argparse
import os
import subprocess
import csv
from collections import Counter
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, DistanceMatrix
from scipy.spatial.distance import pdist, squareform
import scipy.stats as stats
from scipy.stats import pearsonr
from skbio.diversity import alpha_diversity
from scipy.stats import ttest_ind
import time
import sys

def get_user_input(prompt):
    return input(prompt)

def count_abundance(input_file, output_file):
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        abundance_count = Counter(row[1] for row in reader if row)  # Avoid empty rows
    
    with open(output_file, 'w') as outfile:
        for abundance, count in abundance_count.items():
            outfile.write(f"{abundance}\t{count}\n")

def compile_abundance(directory, output_file):
    files = glob.glob(os.path.join(directory, "*_abund.csv"))
    combined_data = {}
    gene_counts = set()
    headers = []

    for file in files:
        sample_name = os.path.basename(file).replace("_abund.csv", "")
        headers.append(sample_name)

        with open(file, 'r') as f:
            next(f)  # Skip header row
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    print(f"Skipping malformed line: {line.strip()}")
                    continue

                gene, abundance = parts
                combined_key = f"{gene},{sample_name}"
                combined_data[combined_key] = abundance
                gene_counts.add(gene)

    genes_sorted = sorted(gene_counts)
    df = pd.DataFrame(index=genes_sorted, columns=headers)

    for gene in genes_sorted:
        for header in headers:
            combined_key = f"{gene},{header}"
            df.at[gene, header] = combined_data.get(combined_key, '0')

    df.to_csv(output_file, sep='\t')
    print(f"Combined data saved to {output_file}")

def generate_heatmap(file_path, folder_name):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Load the data and calculate percentage
    data = pd.read_csv(file_path, sep='\t', index_col=0)
    data_percentage = data.div(data.sum(axis=0), axis=1) * 100

    # Calculate the total abundance across all samples for each gene
    total_abundance = data_percentage.sum(axis=1)

    # Get the top 50 genes with the highest total abundance
    top_50_genes = total_abundance.nlargest(50).index

    # Filter the data to include only the top 50 genes
    filtered_data = data_percentage.loc[top_50_genes]

    # Check if filtered data is empty
    if filtered_data.empty:
        print("No genes available for heatmap generation after selecting top genes.")
        return

    # Determine the appropriate y-axis title based on the file name
    if 'combined_arg.tsv' in file_path:
        yaxis_title = "ARG"
        title = "Heatmap of Dominant ARGs (Top 50)"
    elif 'combined_class.tsv' in file_path:
        yaxis_title = "Class"
        title = "Heatmap of Dominant Class (Top 50)"        
    elif 'combined_mech.tsv' in file_path:
        yaxis_title = "Mechanisms"
        title = "Heatmap of Dominant Mechanisms"        
    else:
        yaxis_title = "ARG/Class/Mechanism"
        title = "ARG/Class/Mechanism"

    # Create the heatmap using Plotly
    fig = px.imshow(
        filtered_data,
        labels=dict(x="Samples", y=yaxis_title, color="Abundance(%)"),
        x=filtered_data.columns,
        y=filtered_data.index,
        color_continuous_scale='mint',
        aspect="auto"
    )
    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title=yaxis_title,
        width=1000,
        height=1200,
        margin=dict(l=200, r=50, t=50, b=100),
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(filtered_data.index))),
            ticktext=filtered_data.index,
            tickfont=dict(size=10)
        )
    )

    # Save the heatmap as a static PDF file
    static_output_path = os.path.join(folder_name, f"heatmap_top_{yaxis_title}.pdf")
    fig.write_image(static_output_path)
    fig.show()
    print(f"Heatmap saved to {static_output_path}")

    # Save the filtered data to a CSV file with the top 50 genes or classes and their percentages
    csv_output_path = os.path.join(folder_name, f"top_50_{yaxis_title}_percentages.csv")
    filtered_data.to_csv(csv_output_path)
    print(f"Top 50 {yaxis_title} data saved to {csv_output_path}")


def analyze_data(file_path, metadata_file, output_dir):
    # Read the TSV file
    data = pd.read_csv(file_path, sep='\t', index_col=0)

    # Load metadata (CSV format)
    metadata = pd.read_csv(metadata_file, index_col=0)

    # Ensure metadata's sample_ID matches the abundance data
    metadata = metadata.loc[data.columns]

    # Normalize the data to get relative abundance for qualitative colormap
    data_normalized = data.div(data.sum(axis=0), axis=1)

    # Merge abundance data with metadata
    data_with_metadata = pd.concat([data_normalized.T, metadata], axis=1)

    # Generate correlation matrix for all numeric data
    numeric_data_with_metadata = data_with_metadata.select_dtypes(include=[np.number])

    # Calculate the Pearson correlation and p-values
    correlation_matrix = numeric_data_with_metadata.corr(method='pearson')
    p_values_matrix = numeric_data_with_metadata.corr(
        method=lambda x, y: stats.pearsonr(x, y)[1]
    ) - np.eye(len(correlation_matrix))

    # Save the correlation matrix to a file (without filtering by p-value)
    correlation_matrix_path = os.path.join(output_dir, 'correlation_matrix.csv')
    correlation_matrix.to_csv(correlation_matrix_path)
    print(f"Correlation matrix saved to {correlation_matrix_path}")

    # Filter the correlation matrix to only include values with p < 0.05
    significant_correlation_matrix = correlation_matrix.where(p_values_matrix < 0.05)

    # Replace NaN values with zero for a complete heatmap
    significant_correlation_matrix = significant_correlation_matrix.fillna(0)

    # Determine dynamic plot size based on the number of genes
    num_genes = len(significant_correlation_matrix)
    plot_size = max(800, num_genes * 15)  # Adjust 15 as needed for better scaling

    # Create heatmap using plotly for the significant correlation matrix
    fig = px.imshow(
        significant_correlation_matrix, 
        color_continuous_scale='RdBu_r', 
        aspect="auto", 
        title="Correlation Plot (p < 0.05)",
        labels=dict(color="Correlation Coefficient")
    )
    fig.update_layout(
        height=plot_size, 
        width=plot_size, 
        title_x=0.5,
        xaxis_tickangle=-45,  # Rotate tick labels for better readability
        font=dict(size=10)  # Adjust font size for better readability
    )

    # Save the figure
    correlation_fig_path = os.path.join(output_dir, 'correlation_matrix_significant.pdf')
    fig.write_image(correlation_fig_path)
    fig.show()
    print(f"Correlation matrix plot saved to {correlation_fig_path}")


    # Perform PCoA and PERMANOVA analysis
    distance_matrix = pdist(data_normalized.T, metric='braycurtis')
    distance_matrix = squareform(distance_matrix)
    dm = DistanceMatrix(distance_matrix, data_with_metadata.index)

    # Loop through each category in metadata to perform PERMANOVA
    for category in metadata.columns:
        permanova_results = permanova(dm, metadata[category])
        print(f"PERMANOVA results for {category}:\n", permanova_results)
        
        # PCoA Plot with sample ID as labels
        pcoa_results = pcoa(dm)
        pcoa_samples = pcoa_results.samples.copy()
        pcoa_samples['Sample_ID'] = data_with_metadata.index
        fig = px.scatter(
            pcoa_samples, 
            x='PC1', 
            y='PC2',
            color=metadata[category], 
            hover_name='Sample_ID',  # Show Sample_ID on hover
            title=f"PCoA (PERMANOVA) - {category}",
            labels={
               'PC1': f"PC1 ({pcoa_results.proportion_explained[0]:.2%})",
               'PC2': f"PC2 ({pcoa_results.proportion_explained[1]:.2%})"
            },
            color_discrete_sequence=px.colors.qualitative.Vivid
        )

        fig.update_traces(marker=dict(size=12, line=dict(width=2, color='DarkSlateGrey')))
        pcoa_fig_path = os.path.join(output_dir, f'pcoa_{category}.pdf')
        fig.write_image(pcoa_fig_path)
        fig.show()
        print(f"PCoA plot saved for {category} to {pcoa_fig_path}")

    # Core genes figure
    
    if 'combined_class.tsv' in file_path:
        title = "Core Class (#Hits)"
    else:
        title = "Core ARGs (#Hits)"
        
    core_genes = data[(data > 0).sum(axis=1) == len(data.columns)]
    if core_genes.empty:
        print("No core genes found.")
    else:
        fig = ff.create_annotated_heatmap(
            z=core_genes.values, x=core_genes.columns.tolist(),
            y=core_genes.index.tolist(), colorscale='aggrnyl',
            annotation_text=core_genes.round(2).values, hoverinfo="z"
        )
        fig.update_layout(
            title_text=title, title_x=0.5,
            height=800, width=800
        )
        core_genes_fig_path = os.path.join(output_dir, 'core_genes_heatmap.pdf')
        fig.write_image(core_genes_fig_path)
        fig.show()
        print(f"Core genes heatmap saved to {core_genes_fig_path}")

    # Alpha diversity calculation and plotting
    diversity_metrics = alpha_diversity('shannon', data.T, ids=data.columns)
    diversity_metadata = pd.DataFrame(diversity_metrics, index=data.columns, columns=['Shannon_Index'])

    # Add metadata to diversity results for plotting
    for category in metadata.columns:
        diversity_metadata[category] = metadata[category]

    # Determine the appropriate title based on the file name
    if 'combined_class.tsv' in file_path:
        title = "Alpha Diversity (Class)"
    else:
        title = "Alpha Diversity (ARG)"

    # Iterate through each category for plotting
    for category in metadata.columns:
        fig = px.box(
           diversity_metadata, 
           x=category, 
           y='Shannon_Index', 
           points="all",
           hover_data=['Shannon_Index', diversity_metadata.index],  # Include Shannon index and sample IDs in hover data
           title=title,
           color=category, 
           color_discrete_sequence=px.colors.qualitative.Set1
    )
        
        fig.update_traces(quartilemethod="linear", jitter=0.3, marker=dict(size=8, opacity=0.6))
        fig.update_layout(width=1000, height=600)
        fig.write_image(os.path.join(output_dir, f'alpha_diversity_{category}.pdf'))
        fig.show()
        print(f"Alpha diversity plot saved for {category}.")

def run_blast_commands(input_dir, output_dir, class_output_dir, mech_output_dir, gene_db, class_db, mech_db, threads, manifest_file):
    
    # Use the manifest file provided via the command line argument
    with open(manifest_file, 'r') as manifest:
        for line in manifest:
            input_file, output_file = line.strip().split(',')
            input_path = os.path.join(input_dir, input_file)
            output_path = os.path.join(output_dir, output_file)
            class_output_path = os.path.join(class_output_dir, output_file)
            mech_output_path = os.path.join(mech_output_dir, output_file)

            gene_blast_cmd = f"diamond blastp --query {input_path} --db {gene_db} --outfmt 6 --id 60 --threads {threads} --max-target-seqs 1 --out {output_path} --evalue 0.00001"
            class_blast_cmd = f"diamond blastp --query {input_path} --db {class_db} --outfmt 6 --id 60 --threads {threads} --max-target-seqs 1 --out {class_output_path} --evalue 0.00001"
            mech_blast_cmd = f"diamond blastp --query {input_path} --db {mech_db} --outfmt 6 --id 60 --threads {threads} --max-target-seqs 1 --out {mech_output_path} --evalue 0.00001"

            print(f"Running BLAST for gene database: {gene_blast_cmd}")
            subprocess.run(gene_blast_cmd, shell=True, check=True)
            print(f"Running BLAST for class database: {class_blast_cmd}")
            subprocess.run(class_blast_cmd, shell=True, check=True)
            print(f"Running BLAST for mechanism database: {mech_blast_cmd}")
            subprocess.run(mech_blast_cmd, shell=True, check=True)

            count_abundance(output_path, output_path.replace("blast.csv", "abund.csv"))
            count_abundance(class_output_path, class_output_path.replace("blast.csv", "abund.csv"))
            count_abundance(mech_output_path, mech_output_path.replace("blast.csv", "abund.csv"))
    
def parse_args():
    parser = argparse.ArgumentParser(description="Genomic Resistance Annotation and Visualization Environment (GRAVE), GRAVE_yard will help you to annotate Resistance in Peace (RIP)")
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input files.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory for Resistance gene-specific output.')
    parser.add_argument('--class_output_dir', type=str, required=True, help='Directory for Antibiotic class-specific output.')
    parser.add_argument('--mech_output_dir', type=str, required=True, help='Directory for Resistance mechanism-specific output.')
    parser.add_argument('--manifest_file', type=str, required=True, help='CSV file containing input and output file names (manifest3.csv)')
    parser.add_argument('--metadata_file', type=str, required=True, help='CSV file containing metadata.')
    parser.add_argument('--gene_db', type=str, required=True, help='Path to gene database (card_gene_mmyyyy.dmnd) for BLAST.')
    parser.add_argument('--class_db', type=str, required=True, help='Path to class database (card_class_mmyyyy.dmnd) for BLAST.')
    parser.add_argument('--mech_db', type=str, required=True, help='Path to mechanism database (card_mech_mmyyyy.dmnd) for BLAST.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use in BLAST (optional, default=4).')
    parser.add_argument('--skip_blast', action='store_true', help='Skip running BLAST and abundance counting (optional).')
    parser.add_argument('--skip_abundance', action='store_true', help='Skip abundance counting after BLAST (optional).')
    return parser.parse_args()
    
def main():
    # Check if the help argument is being called before proceeding
    if '-h' not in sys.argv and '--help' not in sys.argv:
        print("______________________________________________________________________________________Hello comrades...!! Welcome to the GRAVE yard, now you may have your meal if your sample numbers are more than 30 or you have just upto 20 threads to blast the protein grenades, let it enjoy the data war and in rewards you will get victorious graphical representations and tabulars with war history. ____________________________________________________________________________________________________")
        time.sleep(10)  # Pausing for 10 seconds before continuing

    args = parse_args()

    # Run BLAST and abundance counting if not skipped
    if not args.skip_blast:
        print("Running BLAST and abundance counting...")
        run_blast_commands(args.input_dir, args.output_dir, args.class_output_dir, args.mech_output_dir,
                           args.gene_db, args.class_db, args.mech_db, args.threads, args.manifest_file)

    # If skipping BLAST, ensure abundance is still counted from existing BLAST output files
    if args.skip_blast and not args.skip_abundance:
        print("Skipping BLAST, counting abundance from existing BLAST output...")
        manifest_file = args.manifest_file
        with open(manifest_file, 'r') as manifest:
            for line in manifest:
                input_file, output_file = line.strip().split(',')
                output_path = os.path.join(args.output_dir, output_file)
                class_output_path = os.path.join(args.class_output_dir, output_file)
                mech_output_path = os.path.join(args.mech_output_dir, output_file)

                # Count abundance for the gene, class, and mechanism outputs
                count_abundance(output_path, output_path.replace("blast.csv", "abund.csv"))
                count_abundance(class_output_path, class_output_path.replace("blast.csv", "abund.csv"))
                count_abundance(mech_output_path, mech_output_path.replace("blast.csv", "abund.csv"))

    # Compile abundance, generate heatmaps, qualitative colormaps, and perform analyses
    compile_abundance(args.output_dir, os.path.join(args.output_dir, "combined_arg.tsv"))
    generate_heatmap(os.path.join(args.output_dir, "combined_arg.tsv"), args.output_dir)

    compile_abundance(args.class_output_dir, os.path.join(args.class_output_dir, "combined_class.tsv"))
    generate_heatmap(os.path.join(args.class_output_dir, "combined_class.tsv"), args.class_output_dir)

    compile_abundance(args.mech_output_dir, os.path.join(args.mech_output_dir, "combined_mech.tsv"))
    generate_heatmap(os.path.join(args.mech_output_dir, "combined_mech.tsv"), args.mech_output_dir)

    # Perform statistical analysis and visualizations
    analyze_data(os.path.join(args.output_dir, "combined_arg.tsv"), args.metadata_file, args.output_dir)
    analyze_data(os.path.join(args.class_output_dir, "combined_class.tsv"), args.metadata_file, args.class_output_dir)

if __name__ == "__main__":
    main()

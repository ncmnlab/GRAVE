import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def plot_spider_web_chart(data, output_dir):
    # The columns are now the samples (set as index)
    samples = data.columns.tolist()
    
    # Number of samples (this will be the number of axes on the spider web chart)
    num_samples = len(samples)

    # Angles for each sample (equal spacing around the circle)
    angles = [n / float(num_samples) * 2 * pi for n in range(num_samples)]
    angles += angles[:1]  # close the spider web chart circle

    # Initialize spider web chart
    fig, ax = plt.subplots(figsize=(16, 16), subplot_kw={'polar': True})

    # Define a constant offset distance for uniform label placement
    label_offset_distance = 1  # Slightly outside the circle's edge (multiplier for uniform distance)

    # Draw one axis per sample and add labels at the end of each radial axis
    for i, angle in enumerate(angles[:-1]):
        angle_deg = np.degrees(angle)  # Convert radian to degrees

        # Determine the label rotation based on angle for radial alignment
        rotation = angle_deg if angle_deg <= 90 or angle_deg >= 270 else angle_deg + 180
        horizontal_alignment = 'left' if angle_deg <= 90 or angle_deg >= 270 else 'right'

        # Determine vertical alignment logically based on the angle
        if 0 < angle_deg < 180:
            vertical_alignment = 'bottom'  # Top part of the circle
        elif angle_deg == 0 or angle_deg == 180:
            vertical_alignment = 'center'  # Center for top and bottom
        else:
            vertical_alignment = 'top'  # Bottom part of the circle

        # Place the label just outside the radial axis, exactly where the y-axis ends (100% mark on the y-axis)
        ax.text(angle, label_offset_distance * 100, samples[i], size=9, 
                horizontalalignment=horizontal_alignment, 
                verticalalignment=vertical_alignment, 
                rotation=rotation)

    # Set up color cycle for the plot using the updated colormap method
    colors = plt.colormaps['tab10']  # Updated colormap access

    # Plot each mechanism's data
    for i, mechanism in enumerate(data.index):
        # Values for the mechanism (convert to a list and close the circle by adding the first value to the end)
        values = data.loc[mechanism].tolist()
        values += values[:1]  # Close the spider web chart circle

        # Plot data with lines connecting the points, including connection from last to first
        ax.plot(angles, values, linewidth=2, linestyle='-', marker='o', markersize=8, label=mechanism, color=colors(i))

    # Remove radial degree labels (yticklabels and xticklabels)
    ax.set_xticklabels([])  # Remove angle labels (0, 45, etc.)

    # Make the grid circular and more prominent
    ax.yaxis.grid(True, color='grey', linestyle='--', linewidth=0.5)
    ax.xaxis.grid(True, color='grey', linestyle='--', linewidth=0.5)

    # Set the radial ticks (one per sample) and the y-axis ticks dynamically
    ax.set_xticks(angles[:-1])  # Add x-axis ticks for each sample
    ax.set_yticks(np.linspace(0, 100, 10))  # Keep the radial labels based on data scale

    # Set grid to circular shape
    ax.set_rlabel_position(0)
    plt.yticks(np.linspace(20, 100, 5), color="black", size=10)  # Adjust the radial labels (optional)
    plt.ylim(0, 100)

    # Adjust the legend position so it does not overlap the chart
    plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.05), fontsize=10)

    # Add title
    plt.title(f'Spider Web plot for All Mechanisms', size=15, color='blue', y=1.1)

    # Save figure
    plt.tight_layout()
    output_path = f"{output_dir}combined_spider_web.pdf"
    plt.savefig(output_path)
    plt.close()
    print(f"Combined spider web chart saved at {output_path}")

def main():
    # Command-line arguments
    parser = argparse.ArgumentParser(description="Welcome into the GRAVE (Genomic Resistance Annotation and Visualization Environmnet), This is GRAVE_Web for becoming spiderman to generate spider web")
    parser.add_argument('--compiled_abundance_file', type=str, required=True, help='TSV file with compiled gene abundance.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the combined spider web chart.')
    args = parser.parse_args()

    # Load the abundance file
    data = pd.read_csv(args.compiled_abundance_file, sep='\t', index_col=0)

    # Normalize the data (optional: if you need to convert abundance to percentage)
    data = data.div(data.sum(axis=0), axis=1) * 100

    # Plot combined spider web chart
    plot_spider_web_chart(data, args.output_dir)

if __name__ == "__main__":
    main()


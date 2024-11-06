import os
import pandas as pd
import subprocess
import argparse
import logging

def parse_args():
    parser = argparse.ArgumentParser(description="Assemble WGS/metagenome raw reads using Unicycler or metaSPAdes.")
    parser.add_argument('--manifest_file', type=str, required=True, help='Path to the manifest file containing R1, R2, and output paths.')
    parser.add_argument('--threads', type=int, required=True, help='Number of threads to use.')
    parser.add_argument('--type', type=str, choices=['WGS', 'metagenome'], required=True, help='Type of assembly: WGS or metagenome.')
    parser.add_argument('--log_file', type=str, required=True, help='Path to the log file to store execution details.')
    return parser.parse_args()

def get_assembly_size(output_dir):
    """Get the size of the assembly output in MB."""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(output_dir):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size / (1024 * 1024)  # Convert bytes to MB

def main():
    # Parse command-line arguments
    args = parse_args()

    # Set up logging
    logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s - %(message)s')

    # Welcome message
    logging.info("Namaste! AMR Assemblers")
    logging.info("This script is for assembling WGS/metagenome raw reads using Unicycler/metaSPAdes")

    # List to store final assembly results
    assembly_results = []

    # Read the manifest file and iterate over each line
    manifest = pd.read_csv(args.manifest_file, header=None)
    for index, row in manifest.iterrows():
        input_file_R1 = row[0]
        input_file_R2 = row[1]
        output_dir = row[2]

        # Determine which assembler to use based on input type
        if args.type == 'WGS':
            command = f"unicycler -1 {input_file_R1} -2 {input_file_R2} -o {output_dir} -t {args.threads}"
            logging.info(f"Running Unicycler for {input_file_R1} and {input_file_R2}...")
        elif args.type == 'metagenome':
            command = f"metaspades.py -1 {input_file_R1} -2 {input_file_R2} -o {output_dir} -t {args.threads}"
            logging.info(f"Running metaSPAdes for {input_file_R1} and {input_file_R2}...")

        # Execute the command and handle errors
        try:
            subprocess.run(command, shell=True, check=True)
            logging.info(f"Assembly completed for output: {output_dir}")

            # Calculate the assembly size
            assembly_size = get_assembly_size(output_dir)
            logging.info(f"Directory size for {output_dir}: {assembly_size:.2f} MB")

            # Append result to the assembly results list
            assembly_results.append({
                "sample": os.path.basename(output_dir),
                "status": "Success",
                "size_mb": assembly_size
            })

        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred while running the command for {output_dir}: {e}")
            # Append failure to the assembly results list
            assembly_results.append({
                "sample": os.path.basename(output_dir),
                "status": "Failed",
                "size_mb": 0
            })

    # Final summary of all assemblies
    logging.info("Final assembly summary:")
    for result in assembly_results:
        logging.info(f"Sample: {result['sample']}, Status: {result['status']}, Assembly Size: {result['size_mb']:.2f} MB")

    logging.info("Avengers, Assembled")

if __name__ == "__main__":
    main()


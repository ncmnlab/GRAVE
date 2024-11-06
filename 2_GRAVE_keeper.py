import argparse
import os
import subprocess
import logging

def setup_logging(output_dir):
    """Set up logging to a file in the output directory."""
    log_file = os.path.join(output_dir, 'predation.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Logging initialized.")
    return log_file

def log_file_info(input_path, output_path):
    """Log the size of input and output files."""
    try:
        input_size = os.path.getsize(input_path) / (1024 * 1024)  # in MB
        output_size = os.path.getsize(output_path) / (1024 * 1024)  # in MB
        logging.info(f"Input file size: {input_size:.2f} MB, Output file size: {output_size:.2f} MB")
        return input_size, output_size
    except FileNotFoundError:
        logging.error(f"File not found: {output_path}")
        return None, None

def run_prodigal(input_dir, output_dir, manifest, seq_type):
    print("hello GRAVE keepers")
    print("This is the tool for predating proteins, enjoy your meal.....:)")
    logging.info("GRAVE keepers tool started. This is the tool for predating proteins, enjoy your meal.....:)")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")

    with open(manifest, 'r') as manifest_file, open(os.path.join(output_dir, 'predation_status.txt'), 'w') as status_file:
        for line in manifest_file:
            input_file, output_file = line.strip().split(',')
            
            input_path = os.path.join(input_dir, input_file)
            output_path = os.path.join(output_dir, output_file)
            
            if seq_type == 'WGS':
                command = ['prodigal', '-i', input_path, '-o', os.path.join(output_dir, 'out.txt'), '-a', output_path]
            elif seq_type == 'metagenome':
                command = ['prodigal', '-i', input_path, '-o', os.path.join(output_dir, 'out.txt'), '-a', output_path, '-p', 'meta']
            else:
                logging.error(f"Unknown type: {seq_type}")
                return

            # Run the prodigal command and capture output
            try:
                subprocess.run(command, check=True)
                logging.info(f"Successfully processed: {input_file}")
                status_file.write(f"{input_file},{output_file},Success\n")
                
                # Log file sizes
                input_size, output_size = log_file_info(input_path, output_path)
                if input_size is not None and output_size is not None:
                    status_file.write(f"Input size: {input_size:.2f} MB, Output size: {output_size:.2f} MB\n")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running prodigal on {input_file}: {e}")
                status_file.write(f"{input_file},{output_file},Failed\n")

    logging.info("Tremendous predation, bro....")
    print("Tremendous predation, bro....")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some files using Prodigal.")
    
    parser.add_argument('--input_dir', required=True, help="Input directory containing the files.")
    parser.add_argument('--output_dir', required=True, help="Output directory for the results.")
    parser.add_argument('--manifest', required=True, help="Path to the manifest file.")
    parser.add_argument('--type', required=True, choices=['WGS', 'metagenome'], help="Specify 'WGS' or 'metagenome'.")
    
    args = parser.parse_args()
    
    log_file = setup_logging(args.output_dir)
    
    run_prodigal(args.input_dir, args.output_dir, args.manifest, args.type)


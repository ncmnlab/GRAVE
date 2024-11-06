import argparse
import os
import csv

def process_files(input_dir, output_dir, manifest_file, replace_option):
    # Check if the manifest file exists
    if not os.path.isfile(manifest_file):
        print(f"The manifest file '{manifest_file}' does not exist.")
        return

    # Read each line from the manifest and process only files
    with open(manifest_file, 'r') as manifest:
        reader = csv.reader(manifest)
        for line in reader:
            if len(line) != 2:
                print(f"Skipping invalid line in manifest: {line}")
                continue

            input_file = line[0].strip()
            output_file = line[1].strip()

            input_path = os.path.join(input_dir, input_file)
            output_path = os.path.join(output_dir, output_file)

            print(f"Input path: {input_path}")
            print(f"Output path: {output_path}")

            # Check if the input is a valid file and not a directory
            if os.path.isfile(input_path):
                print(f"Processing file: {input_path}")
                with open(input_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
                    reader = csv.reader(infile, delimiter='\t')
                    writer = csv.writer(outfile, delimiter='\t')

                    for row in reader:
                        # Replace semicolon-separated classes in the second column with the chosen option
                        if ';' in row[1]:
                            row[1] = replace_option
                        writer.writerow(row)

                print(f"Processing of {input_path} completed. Output saved to {output_path}.")
            else:
                print(f"Skipping {input_path}, not a valid file.")

    print("All files processed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CSV files and replace semicolon-separated classes.")
    parser.add_argument("--input_dir", required=True, help="Input directory containing the CSV files.")
    parser.add_argument("--output_dir", required=True, help="Output directory to save the processed files.")
    parser.add_argument("--manifest_file", required=True, help="Manifest file (format: input_file,output_file).")
    parser.add_argument("--replace_option", required=True, choices=["multiclass_resistant", "multi_mode"],
                        help="Replacement text for semicolon-separated classes.")

    args = parser.parse_args()

    process_files(args.input_dir, args.output_dir, args.manifest_file, args.replace_option)

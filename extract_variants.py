import os
import re
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse

def sanitize_column(value):
    """
    Sanitize column value to remove any tabs or newlines

    Args:
        value (str): Input string value

    Returns:
        str: Cleaned string value
    """
    return str(value).replace('\t', ' ').replace('\n', ' ').strip()

def process_file(filepath):
    """
    Process a single file to extract rows with 'rs' and transform path/column separators

    Args:
        filepath (str): Path to the input file

    Returns:
        list: Processed rows with transformed separators and added patient/sample columns
    """
    try:
        # Extract patient ID from filepath (assumes filepath follows the pattern JH-2-XXX/...)
        patient_match = re.search(r'(JH-2-\d+)', filepath)
        patient = sanitize_column(patient_match.group(1) if patient_match else 'Unknown')

        # Extract sample name (assumes it's the filename before .filtered...)
        sample_match = re.search(r'([^/]+)\.filtered', filepath)
        sample = sanitize_column(sample_match.group(1) if sample_match else 'Unknown')

        with open(filepath, 'r') as f:
            # Filter rows containing 'rs' and transform separators
            processed_rows = []
            for line in f:
                if 'rs' in line:
                    # Split the line by tabs first
                    parts = line.strip().split('\t')

                    # Ensure we have enough parts and 'rs' is present
                    if len(parts) >= 11 and any('rs' in part for part in parts):
                        # Find the index of the 'rs' variant
                        rs_index = next((i for i, part in enumerate(parts) if 'rs' in part), -1)

                        if rs_index != -1:
                            # Reconstruct the row with specific column mapping
                            reconstructed_row = (
                                f"{patient}\t{sample}\t" +  # patient and sample
                                f"{parts[0]}\t" +           # chromosome
                                f"{parts[1]}\t{parts[2]}\t" +  # start and end
                                f"{parts[3]}\t{parts[4]}\t" +  # ref and alt
                                f"{parts[5]}\t{parts[6]}\t" +  # type and gene
                                f".\t.\t.\t" +             # placeholder for d1, d2, d3
                                f"{parts[rs_index]}"       # variant (rs number)
                            )
                            processed_rows.append(reconstructed_row)
            return processed_rows
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return []

def main(input_pattern, output_file, num_threads=16):
    """
    Main function to process multiple files using thread pool

    Args:
        input_pattern (str): File path pattern to match (e.g., 'JH-2-*/*.filtered.vcf*.txt')
        output_file (str): Path to save the output
        num_threads (int): Number of threads to use
    """
    # Find all matching files
    input_files = glob.glob(input_pattern)

    if not input_files:
        print(f"No files found matching pattern: {input_pattern}")
        return

    print(f"Found {len(input_files)} files to process")

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit file processing tasks
        futures = {executor.submit(process_file, filepath): filepath for filepath in input_files}

        # Collect results
        all_rows = []
        for future in as_completed(futures):
            rows = future.result()
            all_rows.extend(rows)

    # Write results to output file
    # Add header with column names
    header = "patient\tsample\tchr\tstart\tend\tref\talt\ttype\tgene\td1\td2\td3\tvariant"

    with open(output_file, 'w') as outf:
        outf.write(header + '\n')
        # Remove any rows that don't have exactly 13 columns after processing
        valid_rows = [row for row in all_rows if len(row.split('\t')) == 13]
        outf.write('\n'.join(valid_rows))

    print(f"Processed {len(all_rows)} rows. Wrote {len(valid_rows)} valid rows. Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files to extract 'rs' variants")
    parser.add_argument('--input', default='JH-2-*/*.filtered.vcf*.txt',
                        help='Input file pattern to match')
    parser.add_argument('--output', default='complete_variants.tsv',
                        help='Output file path')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of threads to use')

    args = parser.parse_args()

    main(args.input, args.output, args.threads)

#!/usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Filter out specific gene IDs from a count matrix TSV.")
    parser.add_argument("tsv_file", help="Input count matrix (TSV format)")
    parser.add_argument("filter_list", help="File containing list of gene IDs to filter out (one per line)")
    
    args = parser.parse_args()
    
    try:
        # Load gene IDs to filter
        with open(args.filter_list, 'r') as f:
            ids_to_filter = {line.strip() for line in f if line.strip()}
            
        # Read and filter the TSV line by line to preserve exact formatting.
        # This avoids pandas numeric parsing which can alter values (e.g., adding .0 to integers).
        with open(args.tsv_file, 'r') as f:
            header = f.readline()
            if not header:
                return
            
            # Find the index of 'gene_id' column
            cols = header.rstrip('\r\n').split('\t')
            if 'gene_id' not in cols:
                print(f"Error: Column 'gene_id' not found in {args.tsv_file}", file=sys.stderr)
                sys.exit(1)
            
            gene_id_idx = cols.index('gene_id')
            
            # Output the header exactly as it was
            sys.stdout.write(header)
            
            # Filter and output records
            for line in f:
                if not line.strip():
                    continue
                # Split line to check the gene_id
                row = line.rstrip('\r\n').split('\t')
                if len(row) > gene_id_idx:
                    if row[gene_id_idx] not in ids_to_filter:
                        # Write the original line string to preserve exact formatting
                        sys.stdout.write(line)
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

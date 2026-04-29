#!/usr/bin/env python3

import pandas as pd
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
            
        # Read the TSV matrix
        # Using engine='python' or low_memory=False can sometimes help with large files, 
        # but standard read_csv is usually fine for these matrices.
        df = pd.read_csv(args.tsv_file, sep='\t')
        
        # Check if 'gene_id' column exists
        if 'gene_id' not in df.columns:
            print(f"Error: Column 'gene_id' not found in {args.tsv_file}", file=sys.stderr)
            sys.exit(1)
            
        # Filter out records
        # ~ is the logical NOT operator in pandas
        df_filtered = df[~df['gene_id'].isin(ids_to_filter)]
        
        # Output results to stdout (standard tab-separated)
        df_filtered.to_csv(sys.stdout, sep='\t', index=False)
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

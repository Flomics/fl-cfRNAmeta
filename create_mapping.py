#!/usr/bin/env python3

import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Create sample_name to dataset_batch mapping by joining two metadata files.")
    parser.add_argument("sampleinfo", help="Path to sampleinfo_external_and_internal_datasets.tsv")
    parser.add_argument("cfRNA_meta", help="Path to cfRNA-meta_per_sample_metadata.tsv")
    
    args = parser.parse_args()

    # Dictionary to store metadata: 
    # Key: 'sample_name' from cfRNA_meta (contains SRR IDs)
    # Value: 'dataset_batch' (contains study/batch names)
    metadata_map = {}
    
    try:
        with open(args.cfRNA_meta, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                srr_id = row.get('sample_name')
                batch = row.get('dataset_batch')
                if srr_id:
                    metadata_map[srr_id] = batch

        # Output to stdout
        # Columns requested: sample_name (from sampleinfo) and dataset_batch (from cfRNA_meta)
        print("sample_name\tdataset_batch")

        with open(args.sampleinfo, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Join key in sampleinfo is 'sample_id' (contains SRR IDs)
                sample_id = row.get('sample_id')
                if sample_id in metadata_map:
                    name_in_info = row.get('sample_name') # e.g., block_1, chen_1
                    batch = metadata_map[sample_id]       # e.g., block_150bp, chen
                    print(f"{name_in_info}\t{batch}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import sys
import operator
import csv

def main():
    parser = argparse.ArgumentParser(description="Filter rows of a TSV file based on a column value and operator.")
    parser.add_argument("file", help="Input tab-separated file (use '-' for stdin)")
    parser.add_argument("column", help="Column header to filter on")
    parser.add_argument("operator", choices=["==", ">", ">=", "<", "<=", "!="], help="Comparison operator")
    parser.add_argument("value", help="Value to compare against")

    args = parser.parse_args()

    # Map string operators to operator functions
    ops = {
        "==": operator.eq,
        "!=": operator.ne,
        ">":  operator.gt,
        ">=": operator.ge,
        "<":  operator.lt,
        "<=": operator.le
    }
    
    op_func = ops[args.operator]

    try:
        # Open file or use stdin
        input_handle = sys.stdin if args.file == '-' else open(args.file, 'r', newline='')
        
        with input_handle as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            if args.column not in reader.fieldnames:
                print(f"Error: Column '{args.column}' not found.", file=sys.stderr)
                sys.exit(1)

            # Output the header
            print('\t'.join(reader.fieldnames))

            for row in reader:
                raw_val = row[args.column]
                
                # Attempt to convert to float for numeric comparison
                try:
                    val = float(raw_val)
                    comp_val = float(args.value)
                except (ValueError, TypeError):
                    # Fallback to string comparison
                    val = raw_val
                    comp_val = args.value

                if op_func(val, comp_val):
                    print('\t'.join(row[col] for col in reader.fieldnames))

    except FileNotFoundError:
        print(f"Error: File '{args.file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

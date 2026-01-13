import pandas as pd
import sys
import argparse


def merge_svs(input_file, output_file, distant=800):
    try:
        df = pd.read_csv(input_file, sep='\t')
        df['SVType'] = df['SV_support'].str.split(';').str[-1]
        df = df.sort_values(by=['SVType', 'Target_start'])

        merged_rows = []
        i = 0
        while i < len(df):
            current_row = df.iloc[i]
            current_chrom = current_row['#Target_name']
            current_start = current_row['Target_start']
            current_size = current_row['Target_size']
            current_svtype = current_row['SVType']

            j = i + 1
            while j < len(df):
                next_row = df.iloc[j]
                next_chrom = next_row['#Target_name']
                next_start = next_row['Target_start']
                next_size = next_row['Target_size']
                next_svtype = next_row['SVType']
                # Use the 'distant' parameter instead of hardcoded 800
                if (next_chrom == current_chrom and 
                    next_svtype == current_svtype and 
                    abs(next_start - current_start) < distant):
                    current_size += next_size
                    j += 1
                else:
                    break

            current_row['Target_size'] = current_size
            merged_rows.append(current_row)
            i = j

        merged_df = pd.DataFrame(merged_rows)
        merged_df.drop(columns=['SVType'], inplace=True)
        merged_df.to_csv(output_file, sep='\t', na_rep='nan', index=False)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # Initialize argument parser with a description
    parser = argparse.ArgumentParser(
        description='Merge structural variants (SVs) based on proximity and type.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # Shows default values in help
    )
    
    # Required positional arguments
    parser.add_argument('input_file', 
                        help='Path to the input TSV file containing SV data')
    parser.add_argument('output_file', 
                        help='Path to save the merged output TSV file')
    
    # Optional argument with default value
    parser.add_argument('--distant', type=int, default=800,
                        help='Maximum distance (bp) between SVs to be merged')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the merge function with parsed arguments
    merge_svs(args.input_file, args.output_file, args.distant)

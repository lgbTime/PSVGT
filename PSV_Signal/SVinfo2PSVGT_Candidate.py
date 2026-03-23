import argparse
import sys
import os

def transform_tables(input_file, output_file):
    header = "#Target_name\tTarget_start\tTarget_end\tSVlen\tSVID\tSVType\tseq\tmaq\tcluster_size_prevalent\tsv_rate_prevalent"
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            outfile.write(header + "\n")
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) < 4:
                    continue
                
                # Corrected logic: check the 4th column for SV types
                sv_keywords = ["DEL", "INS", "INV", "DUP"]
                
                if parts[3] in sv_keywords:
                    chrom = parts[0]
                    start = int(parts[1])
                    length = int(parts[2])
                    sv_type = parts[3]
                    
                    if sv_type == "INS":
                        end = start + 1
                    else:
                        # DEL, DUP, INV use start + length
                        end = start + length
                
                    sv_id = f"{chrom}:{start}-{end}_{sv_type}={length}"
                    row = [
                        chrom, str(start), str(end), str(length), 
                        sv_id, sv_type, "*", "60", "1", "1.0"
                    ]
                    outfile.write("\t".join(row) + "\n")
                
                else:
                    # Logic for Translocation (TRA)
                    chrom1 = parts[0]
                    chr1_start = parts[1]
                    chrom2 = parts[2]
                    chr2_start = parts[3]
                    sv_type = "TRA"
                    
                    sv_id = f"{chrom1}:{chr1_start}_{chrom2}:{chr2_start}"
                    row = [
                        chrom1, str(chr1_start), str(chr2_start), "0", 
                        sv_id, sv_type, "*", "60", "1", "1.0"
                    ]
                    outfile.write("\t".join(row) + "\n")
        print(f"Successfully generated {output_file}")
    except Exception as e:
        print(f"Error: {e}")

def main():
    # Using RawTextHelpFormatter to keep our custom formatting in the help menu
    parser = argparse.ArgumentParser(
        description="Reformat SV Table 1 into PSVGT-compatible Table 2.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Input Format Requirements (Table 1):
------------------------------------
The script expects a 4-column tab or space-delimited file.

1. For Standard SVs (DEL, INS, DUP, INV):
   Format: [Chr]  [SV Start]  [Length]  [Type]
   Example: Chr1   10025486  2129      DEL

2. For Translocations (TRA):
   Format: [ChrA] [StartA] [ChrB]    [StartB]
   Example: Chr1   2000      Chr4      3030320
"""
    )
    parser.add_argument("-i", "--input", required=True, help="Path to SV information table file")
    parser.add_argument("-o", "--output", required=True, help="Path for the formatted output file")
    args = parser.parse_class() if hasattr(parser, 'parse_class') else parser.parse_args()
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
    transform_tables(args.input, args.output)
if __name__ == "__main__":
    main()

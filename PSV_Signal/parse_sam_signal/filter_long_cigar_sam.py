import re
import sys
insam =sys.argv[1]
outsam = sys.argv[2]
MAX_LEN = 265000000
cigar_pattern = re.compile(r'(\d+)([MIDNSHPX=])')

with open(insam, 'r') as fin, open(outsam, 'w') as fout:
    for line in fin:
        if line.startswith('@'):  
            fout.write(line)
            continue
        parts = line.strip().split('\t')
        if len(parts) < 6:  
            continue
        cigar = parts[5]
        valid = True
        for length_str, op in cigar_pattern.findall(cigar):
            if int(length_str) > MAX_LEN:
                valid = False
                break
        if valid:
            fout.write(line)

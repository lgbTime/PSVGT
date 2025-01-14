import sys
import re
### samtools mpileup bam.list |python $1  > CapsMarker.input
def count_bases(sequence):
    # Initialize a dictionary to hold counts
    base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    # Loop through each character in the sequence
    for base in sequence.upper():  # Convert to uppercase to handle lowercase inputs
        if base in base_counts:
            base_counts[base] += 1
    return base_counts

for line in sys.stdin:
    line = re.sub(r'\^.',"", line)
    line = line.replace("$", "")
    info = line.strip().split("\t")
    num_base = count_bases(",".join(info[3:]).upper())
    num_samples = int(len(info[3:])/3)
    chroms, pos=info[0],info[1]
    most_freq, second_freq = sorted(num_base.values())[-1], sorted(num_base.values())[-2]
    marker_most_freq =  round(most_freq / num_samples,3)
    second_most_freq =  round(second_freq / num_samples,3)
    if most_freq != second_freq and second_freq >0 and marker_most_freq < 0.95 :
        most_freq_base =   [key for key, value in num_base.items() if value == most_freq][0]
        second_freq_base = [key for key, value in num_base.items() if value == second_freq][0]
        print(chroms, pos, most_freq_base, second_freq_base,most_freq, second_freq ,marker_most_freq, second_most_freq,sep="\t")
    elif most_freq == second_freq and most_freq >0 and marker_most_freq < 0.95:
        most_freq_base =   [key for key, value in num_base.items() if value == most_freq][0]
        second_freq_base = [key for key, value in num_base.items() if value == second_freq][1]
        print(chroms, pos, most_freq_base, second_freq_base,most_freq, second_freq,marker_most_freq ,second_most_freq,sep="\t") 

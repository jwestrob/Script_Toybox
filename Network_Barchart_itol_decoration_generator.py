
import argparse
import re
from Bio import SeqIO

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate an iTOL decoration file from a FASTA file.")
parser.add_argument('--input_fasta', type=str, required=True, help='Input protein FASTA file.')
parser.add_argument('--output_itol', type=str, required=True, help='Output iTOL decoration file.')
args = parser.parse_args()

# Write the header for the iTOL decoration file
with open(args.output_itol, 'w') as file:
    file.write("DATASET_SIMPLEBAR\n")
    file.write("SEPARATOR SPACE\n")
    file.write("DATASET_LABEL Length\n")
    file.write("COLOR #ff0000\n")
    file.write("FIELD_LABELS Length\n")
    file.write("DATA\n")

# Extract the length from the sequence headers and write the data lines for the iTOL decoration file
for record in SeqIO.parse(args.input_fasta, "fasta"):
    length_match = re.search('__len(\d+)', record.id)
    if length_match:
        length = length_match.group(1)
        with open(args.output_itol, 'a') as file:
            file.write(f"{record.id.replace('|','_')} {length}\n")

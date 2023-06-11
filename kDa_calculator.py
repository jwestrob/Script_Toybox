import argparse
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight

def calculate_molecular_weight(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        mw = float(molecular_weight(record.seq, seq_type='protein'))
        mw_kDa = mw/1000.0  # Convert Da to kDa
        print(f"Protein Name: {record.id}, Molecular Weight: {mw_kDa:.2f} kDa")

def main():
    parser = argparse.ArgumentParser(description='Calculate molecular weights of proteins in a FASTA file.')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file')

    args = parser.parse_args()
    
    calculate_molecular_weight(args.fasta_file)

if __name__ == "__main__":
    main()

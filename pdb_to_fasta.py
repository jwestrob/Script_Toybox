import sys
from collections import OrderedDict

def pdb_to_fasta(pdb_file, fasta_file):
    sequence_dict = OrderedDict()
    with open(pdb_file, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM'):
                chain = line[21]
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()

                if chain not in sequence_dict:
                    sequence_dict[chain] = []

                aa_code = three_to_one(res_name)
                if aa_code and (len(sequence_dict[chain]) == 0 or sequence_dict[chain][-1][0] != res_seq):
                    sequence_dict[chain].append((res_seq, aa_code))

    with open(fasta_file, 'w') as outfile:
        for chain, residues in sequence_dict.items():
            outfile.write(f">{chain}\n")
            outfile.write("".join([res[1] for res in residues]) + "\n")

def three_to_one(res_name):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(res_name, None)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdb_to_fasta.py <input_pdb_file> <output_fasta_file>")
    else:
        pdb_file = sys.argv[1]
        fasta_file = sys.argv[2]
        pdb_to_fasta(pdb_file, fasta_file)

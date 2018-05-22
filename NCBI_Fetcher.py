import sys, os, csv
import argparse
from Bio import Entrez, SeqIO
Entrez.email = 'jwestrob@andrew.cmu.edu'

parser=argparse.ArgumentParser(description='Input a comma-delimited .csv file of genbank accession IDs for FASTA download.')

parser.add_argument('-infile', metavar='Genbank IDs infile,', help="Path to comma-delimited .csv file")

args=parser.parse_args()

infile = str(args.infile)

def import_data():
    with open(infile, 'r') as csvfile:
        id_reader = csv.reader(csvfile, delimiter=",")
        comp_ids = list(id_reader) #Reads in as list of lists
        untreated_ids = [val for sublist in comp_ids for val in sublist] #Flatten that bad boy
        ids = [GB.lstrip(' ') for GB in untreated_ids] #Remove leading whitespace, if any
    return ids


def fetch(singleID):
    handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'fasta', retmode= 'text')
    f = open('%s.fasta' % singleID, 'w')
    f.write(handle.read())
    handle.close()
    f.close()



def main():
    ids = import_data()
    print(ids)
    for GB in ids:
        print(GB)
        fetch(GB)
    print("boogie")

if __name__ == "__main__":
    main()

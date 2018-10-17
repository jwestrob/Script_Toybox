import os, sys
from Bio import SeqIO


"""
Assumes that your input sequences are in standard PATRIC format
(species name and scaffold ID in the description of the FASTA header);
gets rid of the record ID and replaces it with the description
so that when you make a tree you can see the species name in iTOL.
"""


#FASTA file to reformat
PATRIC = sys.argv[1]

new_recs = []

for rec in SeqIO.parse(PATRIC, 'fasta'):
    rec.id = rec.description.split('  ')[1]
    rec.description = ''
    new_recs.append(rec)

if PATRIC.split('.')[-1] == 'fna':
    print("Writing reformatted seqs to: " + '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.fna")
    SeqIO.write(new_recs, '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.fna", 'fasta')
elif PATRIC.split('.')[-1] == 'faa':
    print("Writing reformatted seqs to: " + '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.faa")
    SeqIO.write(new_recs, '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.faa", 'fasta')
else:
    print("Wat kind of file extension is this dawg are u sure u got this from PATRIC")
    print("Check urself")
    print("Writing reformatted QUESTIONABLE seqs to: " + '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.fa")
    SeqIO.write(new_recs, '.'.join(PATRIC.split('.')[:-1]) + "_reformatted.faa", 'fasta')

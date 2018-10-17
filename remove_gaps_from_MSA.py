from Bio import SeqIO
import os, sys, pandas as pd

input_fasta = sys.argv[1]

output_fasta = sys.argv[2]

out_recs = []

for rec in SeqIO.parse(input_fasta, 'fasta'):
    new_rec = rec
    new_rec.seq = new_rec.seq.ungap('-')
    out_recs.append(new_rec)

SeqIO.write(out_recs, output_fasta, 'fasta')

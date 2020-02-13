import os, sys, pandas as pd
from Bio import SeqIO

fasta_with_dupes = sys.argv[1]


new_recs = []
ids_seen = []

for rec in SeqIO.parse(fasta_with_dupes, 'fasta'):
    if rec.id not in ids_seen:
        ids_seen.append(rec.id)
        new_recs.append(rec)

SeqIO.write(new_recs, fasta_with_dupes, 'fasta')

import sys, os, pandas as pd
from Bio import SeqIO

infile = sys.argv[1]

seen_headers = []
good_recs = []


for rec in SeqIO.parse(infile, 'fasta'):
	if rec.id not in seen_headers:
		good_recs.append(rec)
		seen_headers.append(rec.id)
	else:
		new_rec = rec
		new_rec.id = rec.id + '_1'
		seen_headers.append(new_rec.id)
		good_recs.append(new_rec)

SeqIO.write(good_recs, infile, 'fasta')

print("Complete!")

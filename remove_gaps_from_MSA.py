from Bio import SeqIO
import os, sys, pandas as pd
#from pathos.multiprocessing import ProcessingPool as Pool

input_fasta = sys.argv[1]

#threads = sys.argv[2]

#output_fasta = sys.argv[2]

out_recs = []

#p = Pool(int(threads))


def remove_gap(rec):
	new_rec = rec
	new_rec.seq = new_rec.seq.ungap('-')
	return new_rec

out_recs = list(map(remove_gap, SeqIO.parse(input_fasta, 'fasta')))

"""
for rec in SeqIO.parse(input_fasta, 'fasta'):
    new_rec = rec
    new_rec.seq = new_rec.seq.ungap('-')
    out_recs.append(new_rec)
"""
SeqIO.write(out_recs, input_fasta, 'fasta')
print("Complete!")
sys.exit(420)

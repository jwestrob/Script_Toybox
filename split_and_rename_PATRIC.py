import os, sys, pandas as pd
from Bio import SeqIO

input_fasta = sys.argv[1]

#every sequence is going to start with the format 'acc|XXXX00000.1' or whatever.
#the XXXX is a four-letter alphabetic identifier; that's how we'll split
#name information is in the description

seen_organism_ids = []

cur_recs = []

for rec in SeqIO.parse(input_fasta, 'fasta'):
    organism_id = '_'.join(rec.description.split(' | ')[1].rstrip(']').split(' '))

    if organism_id not in seen_organism_ids:
        #End of previous genome! Write that if applicable
        if len(cur_recs) > 0:
            prev_organism_id = cur_recs[0].id.split('|')[0]
            SeqIO.write(cur_recs, prev_organism_id + '.fna', 'fasta')
            cur_recs = []
        seen_organism_ids.append(organism_id)
    new_rec = rec
    organism_name = '_'.join(rec.description.split('[')[1].split(' | ')[0].split(' '))

    new_rec.id = organism_name + '|' + rec.id.split('|')[1]
    cur_recs.append(new_rec)

#Write final fasta
prev_organism_id = cur_recs[0].id.split('|')[0]
SeqIO.write(cur_recs, prev_organism_id + '.fna', 'fasta')

print("Complete!")

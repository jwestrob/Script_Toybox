from Bio import SeqIO
import os, sys

fastafile = sys.argv[1]

def labeler(fastafile):
	new_recs = []
	genome_id = fastafile.split('.fa')[0].split('.fna')[0]
	for rec in SeqIO.parse(fastafile, 'fasta'):
		if rec.id.split('|')[0] == genome_id:
			sys.exit()
		new_rec = rec
		new_rec.id = genome_id + '|' + rec.id
		new_recs.append(new_rec)
	SeqIO.write(new_recs, fastafile, 'fasta')
	return

labeler(fastafile)

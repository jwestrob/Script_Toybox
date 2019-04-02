from Bio import SeqIO
import sys

infile = sys.argv[1]

recs = list(SeqIO.parse(sys.argv[1], 'fasta'))

nonempty = list(filter(lambda x: x.seq != '' and "#" not in x.seq, recs))

SeqIO.write(nonempty, sys.argv[1], 'fasta')

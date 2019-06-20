from Bio import SeqIO
import sys

infile = sys.argv[1]

#recs = list(SeqIO.parse(sys.argv[1], 'fasta'))


#I don't remember why I'm checking to make sure the pound character isn't in the seq, but I trust my past self
nonempty = list(filter(lambda x: x.seq != '' and "#" not in x.seq, SeqIO.parse(sys.argv[1], "fasta")))

SeqIO.write(nonempty, sys.argv[1], 'fasta')

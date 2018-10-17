from Bio import AlignIO
import sys

input = sys.argv[1]

stupid = list(AlignIO.parse(input, 'phylip'))

AlignIO.write(stupid, input.split('.phy')[0] + '.mfaa', 'fasta')

print('boogie')

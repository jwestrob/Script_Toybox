import os, sys
from Bio import SeqIO

file_to_rename = sys.argv[1]

extension = file_to_rename.split('.')[-1]

for rec in SeqIO.parse(file_to_rename, 'fasta'):
	new_id = '_'.join(rec.description.split('[')[1].split(' |')[0].split()).replace("/", '' ).replace( "(", '').replace(')', '')
	break

new_filename = new_id + '.' + extension

os.system('mv ' + file_to_rename + ' ' + new_filename)


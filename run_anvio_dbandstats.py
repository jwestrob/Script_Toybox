import os, sys
import argparse

parser = argparse.ArgumentParser(description="Run anvi'o stuff on either all files in a directory (-a) or sets of files (-s). Remember to activate Anvi'o first!")

parser.add_argument('-a', action='store_true', default=False, help="Run db and stats on all .fa/.fasta files in directory.")
parser.add_argument('-f', metavar='fastas to process', default=None, nargs='?', help="COMMA-SEPARATED list of fastas to operate on.")

args = parser.parse_args()

a = args.a
if args.a is False:
	fastas = args.f
	fastas = fastas.split(',')

if args.a:
	for file in os.listdir('.'):
		if file.split('.')[-1] == 'fa' or file.split('.')[-1] == 'fasta':
			os.system('anvi-script-reformat-fasta ' + file + ' -o ' + file.split('.')[0] + '-fixed.fa -l 1000 --simplify-names')
			os.system('anvi-gen-contigs-database -f ' + file.split('.')[0] + '-fixed.fa -o ' + file.split('.')[0] + ".db -n 'contigs database'")
			os.system('anvi-run-hmms -c ' + file.split('.')[0] + '.db --num-threads 4')
			os.system('anvi-run-ncbi-cogs -c ' + file.split('.')[0] + '.db --num-threads 4')
			#os.system('rm ' + file)
			#os.system('mv ' + file.split('.')[0] + '-fixed.fa ' + file)
		else:
			continue
else:
	for file in fastas:
		os.system('anvi-script-reformat-fasta ' + file + ' -o ' + file.split('.')[0] + '-fixed.fa -l 1000 --simplify-names')
		os.system('anvi-gen-contigs-database -f ' + file.split('.')[0] + '-fixed.fa -o ' + file.split('.')[0] + ".db -n 'contigs database'")
		os.system('anvi-run-hmms -c ' + file.split('.')[0] + '.db --num-threads 4')
		os.system('anvi-run-ncbi-cogs -c ' + file.split('.')[0] + '.db --num-threads 4')
		#os.system('rm ' + file)
		#os.system('mv ' + file.split('.')[0] + '-fixed.fa ' + file)

print('boogie')

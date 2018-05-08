import os, sys

for file in os.listdir('.'):
	if file.split('.')[-1] == 'fa':
		os.system('anvi-script-reformat-fasta ' + file + ' -o ' + file.split('-')[0] + '.fa -l 1000 --simplify-names')
		os.system('anvi-gen-contigs-database -f ' + file.split('-')[0] + '.fa -o ' + file.split('-')[0] + ".db -n 'contigs database'")
		os.system('anvi-run-hmms -c ' + file.split('-')[0] + '.db --num-threads 4')
		os.system('anvi-run-ncbi-cogs -c ' + file.split('-')[0] + '.db --num-threads 4')
		#os.system('rm ' + file)
	else:
		continue

print('boogie')

import os, sys

for filename in os.listdir('.'):
	if filename.split('.')[1] == 'db':
		os.system('anvi-run-ncbi-cogs -c ' + filename + ' -T 11 --search-with blastp')

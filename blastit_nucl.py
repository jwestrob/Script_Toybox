import os,sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='BLASTp without looking up the damn outfmt every time.')

parser.add_argument('-db', metavar='[DB NAME]', help='Name of blastp database')
parser.add_argument('-threads', metavar='[NUM THREADS]', help="Number of threads to use.")
parser.add_argument('-outfile', metavar='[OUTPUT]', help='Name of output log file.')
parser.add_argument('-query', metavar='[Query fasta]', help="Query fasta file to search with.")
parser.add_argument('-pident', metavar='[Percent Identity Threshold]', help='You know what this is')
args = parser.parse_args()

db_name = args.db
query = args.query
output = args.outfile
threads = args.threads
pident = float(args.pident)

print(db_name, query, output, threads, pident)
#sys.exit()



os.system('blastn -perc_identity ' + str(pident) + ' -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send sseq evalue bitscore" -db ' + db_name + ' -num_threads ' + threads + \
		' -query ' + query + ' -out ' + output)

header = ['qseqid', 'sseqid', 'pident', 'qlen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sseq', 'evalue', 'bitscore']

hits = pd.read_csv(output, sep='\t', names=header)

hits.to_csv(output, index=False, sep='\t')

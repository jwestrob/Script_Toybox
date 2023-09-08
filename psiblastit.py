import os, sys
import argparse
import pandas as pd
import subprocess

parser = argparse.ArgumentParser(description='PSI-BLAST without looking up the damn outfmt every time.')

parser.add_argument('-db', metavar='[DB NAME]', help='Name of blastp database', required=True)
parser.add_argument('-threads', metavar='[NUM THREADS]', help="Number of threads to use.", default='1')
parser.add_argument('-outfile', metavar='[OUTPUT]', help='Name of output log file.', required=True)
parser.add_argument('-query', metavar='[Query fasta]', help="Query fasta file to search with.", required=True)
parser.add_argument('-evalue', metavar='[E-VALUE]', help="E-value cutoff.", default='0.01')
parser.add_argument('-iterations', metavar='[NUM ITERATIONS]', help="Number of PSI-BLAST iterations.", default='3')

args = parser.parse_args()

command = [
    'psiblast',
    '-outfmt', '6 sseqid qseqid pident qlen slen length mismatch gapopen qstart qend sstart send sseq evalue bitscore',
    '-db', args.db,
    '-num_threads', args.threads,
    '-query', args.query,
    '-out', args.outfile,
    '-num_iterations', args.iterations,
    '-evalue', args.evalue
]

result = subprocess.run(command, capture_output=True, text=True)

# Check for errors
if result.returncode != 0:
    print("Error running psiblast:")
    print(result.stderr)
    sys.exit(1)

header = ['sseqid', 'qseqid', 'pident', 'qlen', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sseq', 'evalue', 'bitscore']

hits = pd.read_csv(args.outfile, sep='\t', names=header)

hits.to_csv(args.outfile, index=False, sep='\t')

import os, sys, csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

parser=argparse.ArgumentParser(description='Takes parsed/filtered blastn output and pulls out the subject sequences, placing them into a multifasta (.mfa) ready to be aligned with MAFFT.')
parser.add_argument('-logfile', metavar='blast log infile,', help="Path to blastn log file")
args=parser.parse_args()
logfile = str(args.logfile)
with open(logfile, 'r') as infile:
	log_reader = csv.reader(infile, delimiter="\t")
	log_list = list(log_reader)
sseqs = [] #Empty list into which we will be placing the subject sequences
ids = []

for row in log_list:
	if row[0].split('#')[0] != '':
		ids.append(row[1].split("_")[-1])
		sseqs.append(row[-2])

ids_and_seqs = zip(ids, sseqs)

records = []
for id, seq in ids_and_seqs:
	records.append(SeqRecord(Seq(seq), id))

SeqIO.write(records, 'HERES_UR_MFA_LOSER.mfa', "fasta")

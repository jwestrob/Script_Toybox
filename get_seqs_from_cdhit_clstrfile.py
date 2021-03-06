from Bio import SeqIO
import os, sys, csv
import argparse
from pathos.multiprocessing import ProcessingPool as Pool

parser=argparse.ArgumentParser(description='Get IDs from clusters generated by CDHIT.')
parser.add_argument('-i', metavar='[INFILE]', help=".clstr file from CDHIT")
parser.add_argument('-f', metavar='[FASTAFILE]', help="fasta file used in clustering (NOT CDHIT'S OUTPUT FASTA)")
parser.add_argument('-o', metavar='[OUTFILE]', help="Desired prefix of fasta outfiles", default='')
parser.add_argument('-t', metavar='[THREADS]', help="Number of threads to use. Default 1.", default=1)
parser.add_argument('-threshold', metavar='[THRESHOLD]', default=2, help="[OPTIONAL] Minimum size of cluster to print out. Default 2.")

args = parser.parse_args()

infile = str(args.i)

fastafile = str(args.f)

out_prefix = str(args.o)

threads = int(args.t)

thresh = int(args.threshold)

def main():
	p = Pool(threads)
	with open(infile) as f:
		lines = f.readlines()
	#Create list comprehension to store your stuff
	clusters_list = []
	for line in lines:
		if ">Cluster" in line:
			clusters_list.append([])
		else:
			#Get contig ID
			ID = line.split(">")[1].split("...")[0]
			clusters_list[-1].append(ID)

	#for index, cluster in enumerate(clusters_list):
	#	print("Length of cluster " + str(index) + " is " + str(len(cluster)))

	clstr_recs = list(SeqIO.parse(fastafile, 'fasta'))
	clstr_fa_ids = [rec.id for rec in clstr_recs]

	# For each id in clusters_list, get a list of the same dimension composed of the indices in original FASTA
	# for those ids

	print("Generating indices...")
	indices = list(p.map(lambda cluster: list(map(lambda x: clstr_fa_ids.index(x), cluster)), clusters_list))

	#Get recs for each cluster
	print("Getting records...")
	cluster_recs_list = list(p.map(lambda cluster_indices: list(map(lambda x: clstr_recs[x], cluster_indices)), \
														 indices))
	#Filter out clusters below length threshold
	cluster_recs_list_filtered = list(filter(lambda x: len(x) >= thresh, cluster_recs_list))

	#Print that junk to fastas
	for index, cluster in enumerate(cluster_recs_list_filtered):
		for rec in cluster:
			if 'cluster' not in rec.id:
				rec.id = rec.id + '_cluster_' + str(index)
		SeqIO.write(cluster, out_prefix + 'cluster_' + str(index) + '.faa', 'fasta')

	print('boogie')

if __name__ == "__main__":
	main()

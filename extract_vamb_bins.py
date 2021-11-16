import pandas as pd
import os, sys
from Bio import SeqIO
from pathos.multiprocessing import ProcessingPool as Pool
import argparse

parser = argparse.ArgumentParser(
    description='Given a directory of VAMB results, extract large (>500kb) bins as FASTA files')
parser.add_argument('-vambdir', metavar='[PATH TO VAMB OUTPUT DIRECTORY]',
                    help='Directory containing VAMB output.')
parser.add_argument('-fastafile', metavar='[PATH TO SCAFFOLDS FASTA]',
                    help='Path to binned fasta file.')
parser.add_argument('-threads', metavar='[NUM THREADS]', help='Give me lots')


args = parser.parse_args()

vambdir = args.vambdir
scaffolds = args.fastafile
threads = int(args.threads)

p = Pool(threads)

#Load VAMB clusters file
clusters = os.path.join(vambdir, 'clusters.tsv')
clusters_df = pd.read_csv(clusters, sep='\t', header=None)
clusters_df.columns = ['cluster', 'seqid']
clusters_df.seqid = clusters_df.seqid.apply(lambda x: x.split(' ')[0])


seqids_and_lengths = list(p.map(lambda rec: (rec.id, len(rec.seq)), SeqIO.parse(scaffolds, 'fasta')))
seqids = [seqid_and_len[0] for seqid_and_len in seqids_and_lengths]
lengths = [seqid_and_len[1] for seqid_and_len in seqids_and_lengths]

#VAMB fucks up fasta headers; fix with annoyingly complicated split/join operation
#clusters_df.seqid = clusters_df.seqid.apply(lambda x: '_'.join(x.split('_')[0:5]) + ' ' \
#                        + '_'.join(x.split('_')[5:8]) \
#                        + ' ' + '_'.join(x.split('_')[8:])  )

lendict = dict(zip(seqids, lengths))

clusters_df['seqlen'] = clusters_df.seqid.apply(lambda x: lendict[x])

clusters = clusters_df.cluster.unique().tolist()
cluster_lengths = []
clusters_df['cluster_id'] = clusters_df.cluster
clusters_df = clusters_df.set_index('cluster')

def grab_cluster_length(cluster_id, clusters_df=clusters_df):
    return clusters_df.loc[cluster_id].seqlen.sum()

cluster_lengths = list(p.map(grab_cluster_length, clusters))

cluster_length_df = pd.DataFrame({'cluster':clusters, 'length':cluster_lengths})
cluster_length_df = cluster_length_df.sort_values(by='length', ascending=False)

sufficient_size = cluster_length_df[cluster_length_df.length >= 500000]

#os.chdir(vambdir)
os.system('mkdir ' + os.path.join(vambdir, 'fastas'))
os.system('mkdir ' + os.path.join(vambdir, 'idfiles'))

good_clusters = sufficient_size.cluster.tolist()
good_clusters_df = clusters_df[clusters_df.cluster_id.isin(good_clusters)]

good_clusters_df.seqid = good_clusters_df.seqid.apply(lambda x: x.split('_read_length')[0])
idfiles_dir = os.path.join(vambdir, 'idfiles')

for cluster in good_clusters:
    try:
        cluster_id = int(cluster.split('_')[-1])
    except:
        cluster_id = int(cluster)
    this_cluster_df = good_clusters_df[good_clusters_df.cluster_id == cluster].copy()
    with open(os.path.join(idfiles_dir, 'vamb_bin_' + str(cluster) + '.seqids.txt'), 'w') as outfile:
        [outfile.writelines(element + '\n') for element in this_cluster_df.seqid.unique().tolist()]

idfiles = list(map(lambda x: os.path.join(idfiles_dir, x), os.listdir(idfiles_dir)))

for idfile in idfiles:
    bin_name = idfile.split('.seqids')[0].split('/')[-1]
    os.system('cat ' + idfile + ' | pullseq -i ' + scaffolds + ' -N > ' + os.path.join(vambdir, 'fastas') + '/' +  bin_name + '.fna')

p.terminate()
sys.exit(420)

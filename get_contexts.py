from pathos.multiprocessing import ProcessingPool as Pool
import os, sys, csv, pandas as pd
from Bio import SeqIO
import numpy as np
import argparse

parser=argparse.ArgumentParser('Get genomic contexts for hits in a fasta file (scaffolds with hits)')

parser.add_argument('-infile', metavar='[INPUT FASTAFILE]', help="File of sequences to get contexts for")
parser.add_argument('-od', metavar='[OUTPUT DIRECTORY]', help="Directory to put contexts in (will be separated by sample)")
parser.add_argument('-threads', metavar='[NUMBER OF THREADS]', default=1, help="Number of threads to use. Default 1.")
args = parser.parse_args()

infile = args.infile
outdir = args.od
threads = int(args.threads)

def get_context(lanm_id, contexts, hits_ids):
    #print('_id: ', lanm_id)
    sample_id = lanm_id.split('_scaffold')[0]
    #print('sample_id: ', sample_id)
    if contexts.endswith('/'):
        if not os.path.exists(contexts + sample_id):
            os.system('mkdir ' + contexts + sample_id)
    else:
        if not os.path.exists(contexts + '/' + sample_id):
            os.system('mkdir ' + contexts + '/' + sample_id)

    #Get feature ID (contig number, protein number, sep='_')
    contig_feat = '_'.join((lanm_id.split('scaffold_')[1].split('_')[0], lanm_id.split('scaffold_')[1].split('_')[1]))
    print('contig_feat: ', contig_feat)
    #get contig number
    contig_num = contig_feat.split('_')[0]
    print('contig_num: ', contig_num)
    #Get bin ID; check for pesky '_sub' substring that has no obvious meaning
    bin_id = '_'.join(lanm_id.split('_scaffold_')[1].split('_')[2:])

    print('bin_id: ' + bin_id)
    #Directory to find individual bin proteins
    sample_hits_dir = '/home/jwestrob/2016_Angelo/Protein_Fastafiles/' + sample_id + '/individual_bin_fastas/'

    out_recs = []

    #Why exactly do we need a try/catch here?
    try:
        #print(sample_hits_dir + bin_id + '.faa', 'fasta')
        for rec in SeqIO.parse(sample_hits_dir + bin_id + '.faa', 'fasta'):
            if contig_num in rec.id:
                if rec.id in hits_ids:
                    rec.id = rec.id + '_LanM'
                    out_recs.append(rec)
                else:
                    out_recs.append(rec)
    except:
        #print(sample_hits_dir + sample_id + '_' + bin_id + '.faa')
        for rec in SeqIO.parse(list(filter(lambda x: '.faa' in x, os.listdir(sample_hits_dir)))[0], 'fasta'):
            if contig_num in rec.id:
                if rec.id in hits_ids:
                    rec.id = rec.id + '_LanM'
                    out_recs.append(rec)
                else:
                    out_recs.append(rec)

    #There are occasionally multiple hits per bin; this block of code looks to see if there's already
    #a file with that bin's name, and prevents duplicate filename creation (i.e. overwriting)
    if not os.path.exists(contexts + sample_id + '/' + bin_id + '.faa'):
        print(contexts + sample_id + '/' + bin_id + '.faa')
        SeqIO.write(out_recs, contexts + sample_id + '/' + bin_id + '.faa', 'fasta')
    else:
        num_unks = len(list(filter(lambda x: bin_id in x, os.listdir(contexts + sample_id))))
        print(contexts + sample_id + '/' + bin_id + '.' + str(num_unks) + '.faa')
        SeqIO.write(out_recs, contexts + sample_id + '/' + bin_id + '.' + str(num_unks) + '.faa', 'fasta')
    return

def main():

        context = os.path.join(os.getcwd(), outdir)

        if not os.path.exists(context):
            os.system('mkdir ' + context)
        hits_recs = list(SeqIO.parse(os.path.join(os.getcwd(),infile), 'fasta'))
        hits_ids = [rec.id for rec in hits_recs]

        p = Pool(threads)
        p.map(lambda hit_id: get_context(hit_id, context, hits_ids), hits_ids)

        print('boogie')

if __name__ == "__main__":
        main()

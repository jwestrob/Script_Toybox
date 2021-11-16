#!/usr/bin/python

import os, sys, csv, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

t1 = time.time()

parser=argparse.ArgumentParser(description='Takes a single s2b file and its corresponding contigs file and separates it out into genomes.')

parser.add_argument('-s2b', metavar='s2bfile', nargs=1, help="Path to s2b file (SCAFFOLDS2BIN FORMAT ONLY)")
parser.add_argument('-c', metavar='contigs', nargs='?', default=None, help="Contig file (FASTA FORMAT ONLY)")
parser.add_argument('-outdir', nargs='?', default='.', help="Name of directory to throw your output FASTAs into.")

args = parser.parse_args()

#Parse arguments to strings
s2b = str(args.s2b[0])
contig_file = str(args.c)

outdir = args.outdir

def main():

    binfile_df = pd.read_csv(s2b, sep='\t', names=['Contig', 'Bin'])
    binfile_name = s2b

    #Parse the right contigs file
    contigs = list(SeqIO.parse(contig_file, 'fasta'))
    print("Werkin on : ", binfile_name)
    unique_bins = binfile_df.Bin.unique()
    for index, bin in enumerate(unique_bins):
        #Get a reduced dataframe with only the rows corresponding to the bin in question
        bin_red_df = binfile_df[binfile_df["Bin"] == bin]
        #Make an empty list to store SeqRecord objects to put in a FASTA
        bin_records = []
        bin_red_df_contigs = bin_red_df.Contig.unique().tolist()
        for record in contigs:
            if record.id in bin_red_df_contigs:
                bin_records.append(record)
        if len(bin_records) == 0:
            print(binfile_name, bin)
            print(len(bin_red_df))
            print(bin_red_df.head())
            sys.exit()

        SeqIO.write(bin_records, os.path.join(outdir, binfile_name.split('/')[-1] + '-' + bin + '.fasta'), 'fasta')

        print('Wrote ' + binfile_name + ' bin ' + str(index) + ' to ' + os.path.join(outdir, binfile_name.split('/')[-1] + '-' + bin + '.fasta'))















if __name__ == "__main__":
    main()

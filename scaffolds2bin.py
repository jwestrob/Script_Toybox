#!/usr/bin/python

import os, sys, csv, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

t1 = time.time()

parser=argparse.ArgumentParser(description='Takes a series of binning files (fasta format) and reads (fasta format), creates a tab-separated scaffolds2bin file for input to DAStool.')

parser.add_argument('-bd', metavar='bindir', help="Path to directory containing binning files (FASTA FORMAT ONLY)")
#Is this necessary?
#parser.add_argument('-c', metavar='contigs', help="Contig file (FASTA FORMAT ONLY)")
parser.add_argument('-o', metavar='outfile', help="Name of outfile (.scaffolds2bin.tsv will be appended to whatever you put here)")

args = parser.parse_args()

#Parse arguments to strings
bindir = str(args.bd)
#contigs = str(args.c)
outfile = str(args.o)

def main():
    #Declare empty list for SeqRecord objects

    binfile_list = []
    for filename in os.listdir(bindir):
        #Check filename endings to make sure everything is a fasta format-compliant
        #Maybe I'll allow for .mfa files sometime soon; probably not
        if filename.split('.')[-1] == 'fa' or filename.split('.')[-1] == 'fasta':
            binfile_list.append(SeqIO.parse(bindir + filename, 'fasta'))
        else:
            print("Your binfile directory has something in it that isn't a FASTA. Please take care of this. XOXO")
            sys.exit()

    #Now let's make a separate list for the record objects.
    contig_list = []
    #Don't think I needed to use enumerate... better safe than sorry, I suppose.
    for index, generator in enumerate(binfile_list):
        #Add names of each contig in bin file to contig_list[index]
        for seq_record in generator:
            #Make sure to name bins so they play nicely with DAStool specified format
            if index < 10:
                contig_list.append([seq_record.id, 'bin_0'+str(index+1)])
            else:
                contig_list.append([seq_record.id, 'bin_'+str(index+1)])

    with open(outfile + '.scaffolds2bin.tsv', 'w') as tsvout:
        writer = csv.writer(tsvout, delimiter='\t')
        for element in contig_list:
            writer.writerow(element)
    print("Complete! Time: ", time.time()-t1)





if __name__ == "__main__":
    main()

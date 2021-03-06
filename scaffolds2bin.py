#!/usr/bin/python

import os, sys, csv, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

t1 = time.time()

parser=argparse.ArgumentParser(description='Takes a series of binning files (fasta format) and reads (fasta format), creates a tab-separated scaffolds2bin file for input to DAStool.')

parser.add_argument('-fbd', metavar='bindir', nargs='?', default=None, help="Path to directory containing binning files (FASTA FORMAT ONLY)")
parser.add_argument('-sbd', metavar='bindir', nargs='?', default=None, help="Path to directory containing binning files (SCAFFOLDS2BIN FORMAT ONLY)")
parser.add_argument('-c', metavar='contigs', nargs='?', default=None, help="Contig file (FASTA FORMAT ONLY)")
parser.add_argument('-cd', metavar='contigs_dir', nargs='?', default=None, help="Contig file directory (ALSO ONLY FASTA FORMAT PLS)")
parser.add_argument('-r', action='store_true', default=False, help="Wanna turn a .scaffolds2bin file into a .fasta bin file?")
parser.add_argument('-o', metavar='outfile', nargs='?', default=None, help="S2B ONLY: Name of outfile (.scaffolds2bin.tsv will be appended to whatever you put here)")
parser.add_argument('-outdir', nargs='?', default='.', help="B2S only: Name of directory to throw your output FASTAs into.")

args = parser.parse_args()

#Parse arguments to strings
if args.fbd is not None:
    fbindir = str(args.fbd)
if args.sbd is not None:
    sbindir = str(args.sbd)
if args.c is not None:
    contig_file = str(args.c)
if args.cd is not None:
    contig_dir = str(args.cd)
if args.o is not None:
    outfile = str(args.o)
if args.outdir != '.':
    outdir = str(args.outdir).split('/')[0]
reverse = args.r

def main():
    if reverse == False:
        #Declare empty list for SeqRecord objects

        binfile_list = []
        for filename in os.listdir(fbindir):
            #Check filename endings to make sure everything is a fasta format-compliant
            #Maybe I'll allow for .mfa files sometime soon; probably not
            if filename.split('.')[-1] == 'fa' or filename.split('.')[-1] == 'fasta':
                binfile_list.append(SeqIO.parse(bindir + filename, 'fasta'))
            #else:
                #print("Your binfile directory has something in it that isn't a FASTA. Please investigate (this will not impede file conversion)")

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

    elif reverse:
        binfile_list = []
        #Get a list of binfiles (as pandas dfs) with corresponding filename
        for filename in os.listdir(sbindir):
            if filename.split('_')[-1] == 'scaffolds2bin.txt':
                binfile_list.append([pd.read_csv(sbindir + '/' + filename, sep='\t', names=["Contig", "Bin"]), filename])
            #else:
                #print("Found some nonsense. Please evaluate: ", filename)

        for binfile in binfile_list:
            binfile_df = binfile[0]
            binfile_name = binfile[1]

            #Parse the right contigs file
            contigs = list(SeqIO.parse(contig_file, 'fasta'))
            print("Werkin on : ", binfile_name)
            unique_bins = binfile_df.Bin.unique()
            for index, bin in enumerate(unique_bins):
                #Get a reduced dataframe with only the rows corresponding to the bin in question
                bin_red_df = binfile_df[binfile_df["Bin"] == bin]
                #Make an empty list to store SeqRecord objects to put in a FASTA
                bin_records = []
                for record in contigs:
                    if record.id in list(bin_red_df["Contig"]):
                        bin_records.append(record)
                if len(bin_records) == 0:
                    print(binfile_name, bin)
                    print(len(bin_red_df))
                    print(bin_red_df.head())
                    sys.exit()
                SeqIO.write(bin_records, outdir + '/' + binfile_name + '-' + bin + '.fasta', 'fasta')
                print('Wrote ' + binfile_name + ' bin ' + str(index) + ' to ' + binfile_name + '-' + bin + '.fasta')








    else:
        print("You didn't select a mode fool!")
        sys.exit(420)







if __name__ == "__main__":
    main()

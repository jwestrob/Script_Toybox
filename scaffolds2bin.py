import os, sys, csv, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
from pathos.multiprocessing import ProcessingPool as Pool

t1 = time.time()

parser=argparse.ArgumentParser(description='Takes a series of binning files (fasta format) and scaffolds (fasta format), creates a tab-separated scaffolds2bin file for input to DAStool.')

parser.add_argument('-fbd', metavar='bindir', nargs='?', default=None, help="Path to directory containing binning files (FASTA FORMAT ONLY)")
parser.add_argument('-sbd', metavar='s2bfile', nargs='?', default=None, help="Path to a single binning file (SCAFFOLDS2BIN FORMAT ONLY)")
parser.add_argument('-c', metavar='contigs', nargs='?', default=None, help="Contig file (FASTA FORMAT ONLY)")
parser.add_argument('-cd', metavar='contigs_dir', nargs='?', default=None, help="Contig file directory (ALSO ONLY FASTA FORMAT PLS)")
parser.add_argument('-r', action='store_true', default=False, help="Wanna turn a .scaffolds2bin file into a .fasta bin file?")
parser.add_argument('-o', metavar='outfile', nargs='?', default=None, help="S2B ONLY: Name of outfile (.scaffolds2bin.tsv will be appended to whatever you put here)")
parser.add_argument('-outdir', nargs='?', default='.', help="B2S only: Name of directory to throw your output FASTAs into.")
parser.add_argument('-threads', nargs=1, default=1, help="Threads to use for sequence extraction (large input sequences only).")

args = parser.parse_args()

#Parse arguments to strings
if args.fbd is not None:
    fbindir = str(args.fbd)
if args.sbd is not None:
    sbinfile = str(args.sbd)
if args.c is not None:
    contig_file = str(args.c)
if args.cd is not None:
    contig_dir = str(args.cd)
if args.o is not None:
    outfile = str(args.o)
if args.outdir != '.':
    outdir = str(args.outdir).split('/')[0]

threads = int(args.threads[0])
reverse = args.r


def extract_hits(bins_to_contig_lists, outdir, contig_file, threads):
    p = Pool(threads)

    pullseq_tmp = os.path.join(outdir, 'pullseq_ids_tmp')
    if not os.path.exists(pullseq_tmp):
        os.system('mkdir ' + pullseq_tmp)


    def pullseq_by_bin(bin_name, contig_list, contig_file):
        #Generates a file with the names of all the contigs to pull out
        #then provides that to pullseq;
        #parses the resulting fasta output from pullseq and then
        #passes it back.
        with open(os.path.join(pullseq_tmp, bin_name + '.txt'), 'w') as outfile:
            for element in contig_list:
                outfile.writelines(element + '\n')

        os.system('pullseq -i ' + contig_file + ' -n ' + os.path.join(pullseq_tmp, bin_name + '.txt') + ' > ' + os.path.join(outdir, bin_name + '.fasta'))


        return


    p.map(lambda x: pullseq_by_bin(x, bins_to_contig_lists[x], contig_file), bins_to_contig_lists)
    #for bin in bins_to_contig_lists:
    #    pullseq_by_bin(bin, bins_to_contig_lists[bin], contig_file)

    os.system('rm -rf ' + pullseq_tmp)
    p.terminate()
    return
def main():
    if reverse == False:
        #Declare empty list for SeqRecord objects

        binfile_list = []
        for filename in os.listdir(fbindir):
            #Check filename endings to make sure everything is a fasta format-compliant
            #Maybe I'll allow for .mfa files sometime soon; probably not
            
            if filename.split('.')[-1] == 'fa' or filename.split('.')[-1] == 'fasta' or filename.split('.')[-1] == 'fna':
                binfile_list.append(SeqIO.parse(os.path.join(os.path.abspath(fbindir), filename), 'fasta'))
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
        # Read the single binfile as a pandas df
        binfile_df = pd.read_csv(sbinfile, sep='\t', header=None)
        if len(binfile_df.columns) > 2:
            binfile_df = binfile_df[binfile_df.columns[0:2]]
            binfile_df.columns = ["Contig", "Bin"]
        else:
            binfile_df.columns = ["Contig", "Bin"]
        binfile_list.append([binfile_df, sbinfile])

        if not os.path.exists(outdir):
            os.system('mkdir ' + outdir)
        for binfile in binfile_list:
            binfile_df = binfile[0]
            binfile_name = binfile[1]

            #Parse the right contigs file
            contigs = SeqIO.parse(contig_file, 'fasta')
            print("Werkin on : ", binfile_name)
            bins_to_contig_lists = {}
            unique_bins = binfile_df.Bin.unique().tolist()
            for index, bin in enumerate(unique_bins):
                #Get a reduced dataframe with only the rows corresponding to the bin in question
                bin_contigs = binfile_df[binfile_df["Bin"] == bin].Contig.tolist()
                #Make an empty list to store SeqRecord objects to put in a FASTA
                bins_to_contig_lists[bin] = bin_contigs

            extract_hits(bins_to_contig_lists, outdir, contig_file, threads)
            











if __name__ == "__main__":
    main()

import os, sys, csv
import pandas as pd
from Bio import SeqIO, SearchIO, AlignIO
import argparse
from pathos.multiprocessing import ProcessingPool as Pool

parser = argparse.ArgumentParser(description='Scan a given protein FASTA file for a protein of interest using an HMM; \
            return FASTA format hits and a table of hit frequencies.' + '\n' + \
            "If you want to align your fasta output, include the '-align' flag." + '\n' + \
            "If you want to align them using the more stringent MAFFT parameters, include -macc as well.")


parser.add_argument('-hmm', metavar='hmm file', help='PATH to HMMER3-compatible .hmm file. (pls kindly include extension)')
parser.add_argument('-p', metavar='protein fastafile', help='PATH to protein FASTA to search')
parser.add_argument('-fo', metavar='fasta output', help='Name of FASTA output file with HMM hits')
parser.add_argument('-ids', metavar='Contig hits ID file name', default=None, \
                    help='Output contig IDs of HMM hits to file at specified PATH')
parser.add_argument('-align', metavar='Output alignment of retrieved sequences', default=None, \
                    help='Align output sequences and place in specified PATH (uses MAFFT; optionally increase accuracy (-macc))')
parser.add_argument('-macc', action='store_true', default=False, \
                    help='Use this flag to use increased accuracy when aligning with MAFFT.')
parser.add_argument('-ht', metavar='hits table', help='Name of hits table to write.')
parser.add_argument('-PAT',  help='PATRIC format data', \
                    action='store_true', default=False)
parser.add_argument('-gg', help='ggKbase format data', \
                    action='store_true', default=False)
parser.add_argument('-threads', metavar='# of threads', default=1, \
                    help='Threads to use for HMMsearch/MAFFT.')

args = parser.parse_args()

hmmfile = str(args.hmm)
protfile = str(args.p)
fastaout = str(args.fo)

PATRIC = args.PAT
ggkbase = args.gg

if args.ids is not None:
    idfile = str(args.ids)
else:
    idfile = None
if args.align is not None:
    align = str(args.align)
else:
    align = None
threads = int(args.threads)
macc = args.macc
hits_table_out = args.ht

def flatten(inlist):
    return [item for sublist in inlist for item in sublist]

def run_hmmsearch(hmmfile, cwd):
    print('------------------------------------------------------------')
    print("Beginning HMMsearch...")
    print('hmmsearch -o ' + cwd + '/' + hmmfile.split('/')[-1].split('.hmm')[0] + \
            '_hmmsearch.out --notextw --cpu ' + str(threads) + ' ' + hmmfile + \
            ' ' + protfile)
    os.system('hmmsearch -o ' + cwd + '/' + hmmfile.split('/')[-1].split('.hmm')[0] + \
            '_hmmsearch.out  --notextw --cpu ' + str(threads) + ' ' + hmmfile + \
            ' ' + protfile)
    print('------------------------------------------------------------')
    return hmmfile.split('/')[-1].split('.hmm')[0] + '_hmmsearch.out'

def get_hits(infile):
    hits = []
    with open(infile, 'rU') as handle:
        for record in SearchIO.parse(handle, 'hmmer3-text'):
            hits.append(list(record.hit_keys))
    return hits

def write_hits(hits):
    print("Writing HMM hit FASTA contig IDs to " + idfile)
    with open(idfile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for row in hits:
            writer.writerow(row)
    print('------------------------------------------------------------')
    return

def pull_out_seqs(hits):
    print("Pullin out seqs...")
    in_recs = list(SeqIO.parse(protfile, 'fasta'))
    out_recs = []
    for id in hits:
        for rec in in_recs:
            if id in rec.id or id in rec.description:
                out_recs.append(rec)
    print('------------------------------------------------------------')
    return out_recs

def finder(rec, hits):
    for id in hits:
        if id in rec.id or id in rec.description:
            return rec
    return None

def pull_out_seqs_ALT(hits):
    print("Pullin out seqs...")
    out_recs = []
    #print(hits)
    #p = Pool(threads)
    recs = list(SeqIO.parse(protfile, 'fasta'))
    p = Pool(threads)
    out_recs_wNone = list(p.map(lambda r: finder(r, hits), recs))
    out_recs = [x for x in out_recs_wNone if x is not None]
    print('------------------------------------------------------------')
    return out_recs

def align_seqs():
    #Tells MAFFT to align your sequences
    if macc:
        print("Aligning sequences using MAFFT with high accuracy.")
        print("Outputting alignment to " + align)
        os.system('mafft --maxiterate 1000 --localpair --thread ' + str(threads) + ' ' + fastaout + ' > ' + align)
    else:
        print("Aligning sequences using MAFFT with low accuracy.")
        print("Outputting alignment to " + align)
        os.system('mafft --thread ' + str(threads) + ' ' + fastaout + ' > ' + align)
    print('------------------------------------------------------------')
    return

def main():
    cwd = os.getcwd()

    #First things first. Run hmmsearch to get the file to parse
    hmm_output_file = run_hmmsearch(hmmfile, cwd)

    #Cool. Now let's get the IDs for the hits.
    hits = get_hits(hmm_output_file)

    if idfile is not None:
        #Write hits names to file
        write_hits(hits)

    #Eliminate potential nested list bc i don't want to deal with it
    hits = flatten(hits)

    #Now let's get the proteins that correspond to the hits.
    hit_recs = pull_out_seqs_ALT(hits)

    if len(hit_recs) == 0:
        print("No hits detected. Examine outfile to check if anything went wrong.")
        sys.exit(420)
    else:
        print("Writing hit sequences to: " + fastaout)
        SeqIO.write(hit_recs, fastaout, 'fasta')

    if align is not None:
        align_seqs()

    if PATRIC:
        hits_ids = list(map(lambda x: x.split('.peg')[0].split('|')[1], hits))
    elif ggkbase:
        hits_ids = list(map(lambda x: x.split('_scaffold')[0] + '_' +\
                            '_'.join(x.split('_scaffold_')[1].split('_')[2:]), hits))
    else:
        hits_ids = hits

    if hits_ids is not None:
        hits_table = [[i] for i in hits_ids]
    else:
        print("No hits. Exiting...")
        sys.exit(420)

    for element in hits_table:
        element.append(0)



    hits_ids = pd.Series(hits_ids)
    hits_ids_counts = hits_ids.value_counts()



    header = ['Organism_ID', 'num_hits']

    df = pd.DataFrame(hits_table, columns=header)

    df.num_hits = df.apply(lambda row: hits_ids_counts[row['Organism_ID']], axis=1)

    print("Writing hits table to: " + hits_table_out)
    df.to_csv(hits_table_out, sep='\t', index=False)
    print("Complete! Check to make sure things are OK.")


if __name__ == '__main__':
    main()

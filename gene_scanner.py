import os, sys, csv
import pandas as pd
from Bio import SeqIO, SearchIO, AlignIO
import argparse

parser = argparse.ArgumentParser(description='Scan a given set of genomes (in ~/concat) for a protein of interest using an HMM.')



parser.add_argument('-hmm', metavar='hmm file', help='PATH to HMMER3-compatible .hmm file. (pls kindly include extension)')
parser.add_argument('-fd', metavar='fasta directory', help='PATH to directory with protein FASTA files')
parser.add_argument('-c', metavar='concat fastafile', help='PATH to concatenated protein FASTA')
parser.add_argument('-fo', metavar='fasta output', help='Name of FASTA output file with HMM hits')
parser.add_argument('-ids', metavar='Contig hits ID file name', default=None, \
                    help='Output contig IDs of HMM hits to file at specified PATH')
parser.add_argument('-align', metavar='Output alignment of retrieved sequences', default=None, \
                    help='Align output sequences and place in specified PATH (uses MAFFT; optionally increase accuracy (-macc))')
parser.add_argument('-macc', action='store_true', default=False, \
                    help='Use this flag to use increased accuracy when aligning with MAFFT.')
parser.add_argument('-ht', metavar='hits table', help='Name of hits table to write.')
parser.add_argument('-threads', metavar='# of threads', default=1, \
                    help='Threads to use for HMMsearch/MAFFT.')

args = parser.parse_args()

hmmfile = str(args.hmm)
fastadir = str(args.fd)
concat = str(args.c)
fastaout = str(args.fo)
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
            '_hmmsearch.out  --cpu ' + str(threads) + ' ' + hmmfile + \
            ' ' + concat)
    os.system('hmmsearch -o ' + cwd + '/' + hmmfile.split('/')[-1].split('.hmm')[0] + \
            '_hmmsearch.out  --cpu ' + str(threads) + ' ' + hmmfile + \
            ' ' + concat)
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
    in_recs = list(SeqIO.parse(concat, 'fasta'))
    out_recs = []
    for id in hits:
        for rec in in_recs:
            if id in rec.id or id in rec.description:
                out_recs.append(rec)
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
    hit_recs = pull_out_seqs(hits)

    if len(hit_recs) == 0:
        print("No hits detected. Examine outfile to check if anything went wrong.")
        sys.exit(420)
    else:
        print("Writing hit sequences to: " + fastaout)
        SeqIO.write(hit_recs, fastaout, 'fasta')

    if align is not None:
        align_seqs()


    hits_ids = list(map(lambda x: x.split('.peg')[0].split('|')[1], hits))

    hits_table = [[i] for i in hits_ids]

    for element in hits_table:
        element.append(0)

    for filename in os.listdir(fastadir):
        if filename.endswith('.faa'):
            fasta_id = filename.split('.PATRIC.faa')[0]
            if fasta_id in hits_ids:
                idx = hits_ids.index(fasta_id)
                hits_table[idx][1] = hits_ids.count(fasta_id)
                recs = list(SeqIO.parse(filename, 'fasta'))

                if ']' in recs[0].description.split('   ')[2]:
                    organism_name = recs[0].description.split('   ')[2].strip('[').strip(']')
                elif ']' in recs[0].description.split('   ')[1]:
                    organism_name = recs[0].description.split('   ')[1].strip('[').strip(']')
                hits_table[idx].append(organism_name)

    headers = ['Fasta ID', '# of XoxF family HMM hits', 'Organism name']
    df = pd.DataFrame(hits_table, columns=headers)
    print("Writing hits table to: " + hits_table_out)
    df.to_csv(hits_table_out, sep='\t', index=False)
    print("Complete! Check to make sure things are OK.")


if __name__ == '__main__':
    main()

from pathos.multiprocessing import ProcessingPool as Pool
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, SearchIO
from Bio.Alphabet import IUPAC
import os, sys, pandas as pd
from pathlib import Path
from Bio.Seq import Seq
import numpy as np
import subprocess
import argparse

parser = argparse.ArgumentParser(
    description='Given a directory of protein fasta files, extract the best hit for each HMM for each bin; get multifasta of each protein')
parser.add_argument('-fastadir', metavar='[PROTEIN FASTA DIRECTORY]',
                    help="Directory of protein fasta files to scan for ribosomal proteins")
parser.add_argument('-hmmdir', metavar='[HMM DIRECTORY]', help="Directory containing HMM files to scan with")
parser.add_argument('-outdir', metavar='[Output Directory]', default='output', help="Directory to store output files")
parser.add_argument('-evalue', metavar='[EVALUE THRESHOLD]', default=0.01,
        help="Evalue threshold for HMMsearch. Default: 0.01. (I pick the hit with the best e-value; \
        wait until alignment stage to curate)")
parser.add_argument('-threads', metavar='[NUM THREADS]', default=1, help="Number of threads to use")
parser.add_argument('-already_scanned', default=False, action='store_true', help='For if you already ran the HMMs')
parser.add_argument('-no_seqs', default=False, action='store_true', help='Dont pull out sequences to fasta')

def run_hmmsearch(protfile, hmmfile, wd, threshold):
    protein_id = protfile.split('/')[-1].split('.faa')[0]

    print('------------------------------------------------------------')
    print("Beginning HMMsearch...")
    print(protein_id, hmmfile)
    cmd = 'hmmsearch -o ' + wd + '/' + protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + \
          '_hmmsearch.out  --notextw -E ' + str(threshold) + ' --cpu ' + str(1) + ' ' + hmmfile + \
          ' ' + protfile
    print(cmd)
    result = subprocess.getstatusoutput(cmd)
    if result[0] != 0:
        print('HMMsearch error (check for empty sequences in your protein FASTAs)')
        print('protein_id: ', protein_id)
        print('hmmfile: ', hmmfile)
        sys.exit()
    print(result)
    print('------------------------------------------------------------')
    return protein_id + '_' + hmmfile.split('/')[-1].split('.hmm')[0] + '_hmmsearch.out'


def extract_hits_by_outfile(dir, infile):
    hits = []
    e_values = []
    with open(dir + '/' + infile[0], 'r') as handle:
        for record in SearchIO.parse(handle, 'hmmer3-text'):
            hits.append(list(record))

    try:
        hits = hits[0]
    except:
        return

    good_hits = [hit._id for hit in hits]
    e_values = [hit.evalue for hit in hits]
    # If you have more than one hit, go with the hit that has the best e-value
    if len(good_hits) > 1:
        return good_hits[e_values.index(max(e_values))]
    else:
        try:
            return good_hits[0]
        except:
            return


def get_recs_for_hits(hits_ids, hmm, fastadict, fastalist_wpath, fastalist, outdir):
    # flat_hits = [item for sublist in hits_ids for item in sublist]
    hit_recs = []
    for hit in hits_ids:
        if hit is None:
            continue
        else:
            # print(hit)
            if '.peg' in hit:
                genome = fastadict[hit.split('.peg')[0]]
            elif '# ID=' in hit:
                genome = fastadict['_'.join(hit.split(' #')[0].split('_')[0:-2])]
            else:
                num_split = len(hit.split('_'))
                genome = fastadict['_'.join(hit.split('_')[0:num_split-1])]
            recs = list(SeqIO.parse(fastalist_wpath[fastalist.index(genome)], 'fasta'))
            hit_rec = list(filter(lambda x: x.id == hit, recs))[0]
            hit_rec.id = genome + '|' + hit
            hit_recs.append(hit_rec)
    print("Writing hits for: ", hmm)
    SeqIO.write(hit_recs, outdir + '/fastas/' + hmm + '.faa', 'fasta')
    return


def main():
    args = parser.parse_args()
    fastadir = str(Path(args.fastadir).absolute())
    hmmdir = str(Path(args.hmmdir).absolute())
    outdir = str(args.outdir)
    threshold = float(args.evalue)
    threads = int(args.threads)
    already_scanned = args.already_scanned
    no_seqs = args.no_seqs

    p = Pool(threads)

    # Make output directory
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
        outdir = str(Path(outdir).absolute())
    else:
        outdir = str(Path(outdir).absolute())

    # Get list of paths of all fastas
    fastalist_wpath = list(map(lambda file: os.path.join(fastadir, file), os.listdir(fastadir)))

    # Get list of all fastas
    fastalist = list(map(lambda file: file.split('.faa')[0], os.listdir(fastadir)))

    # Get list of paths of all HMM files
    hmmlist_wpath = list(map(lambda file: os.path.join(hmmdir, file), os.listdir(hmmdir)))

    # Get list of all HMMs
    hmmlist = list(map(lambda file: file.split('.hmm')[0], os.listdir(hmmdir)))

    hmm_outfiles = []

    def get_fastaheader_id(fasta):
        for rec in SeqIO.parse(fasta, 'fasta'):
            id = rec.id.split('.peg')[0]
            break
        return id

    fasta_header_ids = list(map(get_fastaheader_id, fastalist_wpath))
    fastadict = dict(zip(fasta_header_ids, fastalist))

    # For each fasta, run all hmms
    if not already_scanned:
        for fastafile in fastalist_wpath:
            fastaoutdir = outdir + '/' + fastafile.split('/')[-1].split('.faa')[0]
            # Make outdir for HMMs
            if not os.path.exists(fastaoutdir):
                os.system('mkdir ' + outdir + '/' + fastafile.split('/')[-1].split('.faa')[0])
            hmm_outfiles.append([])

            # Run all HMMs for fastafile
            hmm_outfiles[-1] = list(p.map(lambda hmmfile: run_hmmsearch(fastafile, hmmfile, outdir, threshold), \
                                          hmmlist_wpath))

            # Move all outfiles to corresponding output directory
            for outfile in hmm_outfiles[-1]:
                os.system('mv ' + outdir + '/' + outfile + ' ' + fastaoutdir)

    # Make directory to store fastas
    if not os.path.exists(outdir + '/' + 'fastas'):
        os.system('mkdir ' + outdir + '/' + 'fastas')

    # Make matrix of zeros to store hits

    hits_by_hmm = []

    # test = extract_hits_by_outfile('/home/jacob/Documents/Berkeley/test_ribosomal_nonsense/943347.4.PATRIC', ['943347.4.PATRIC_RecR_hmmsearch.out'])
    # Declare local version of fastalist because global variables are unavailable to lambda functions (??? python...)
    for hmm in hmmlist:
        print("Extracting hits for: ", hmm)
        relevant_outfiles = []
        for fasta in fastalist:
            relevant_outfiles.append(list(filter(lambda x: hmm in x, os.listdir(outdir + '/' + fasta))))
        # Add HMM outfiles to a list; find these with
        hits_by_hmm.append(list(p.map(lambda relevant_outfile: extract_hits_by_outfile( \
            outdir + '/' + fastalist[relevant_outfiles.index(relevant_outfile)],
            relevant_outfile), relevant_outfiles)))


    print("Making hits matrix...")
    hitstable = np.zeros((len(hmmlist), len(fastalist)))
    # Mark hits in table
    for hmm_idx, hmm in enumerate(hits_by_hmm):
        for genome_idx, genome_hits in enumerate(hmm):
            if type(genome_hits) is list:
                hits = len(genome_hits)
            elif type(genome_hits) is str:
                hits = 1
            if genome_hits is None:
                hitstable[hmm_idx][genome_idx] = 0
            else:
                hitstable[hmm_idx][genome_idx] = hits


    # Grab sequences
    if not no_seqs:
        list(p.map(lambda hits:
                   get_recs_for_hits(hits, hmmlist[hits_by_hmm.index(hits)], fastadict, fastalist_wpath, fastalist,
                                     outdir),
                   hits_by_hmm))

    hits = pd.DataFrame(hitstable).T
    hits.columns = hmmlist
    hits['id'] = fastalist

    cols = list(hits.columns.values)
    cols.pop(cols.index('id'))
    hits = hits[['id'] + cols]
    hits.to_csv(outdir + '/HITSTABLE.tsv', sep='\t', index=False)
    print(hits)
    # recs_by_hmm = list(map(lambda hits: get_recs_for_hits(hits), hits_by_hmm))
    print('boogie')


if __name__ == "__main__":
    main()

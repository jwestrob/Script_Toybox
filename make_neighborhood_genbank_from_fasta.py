import warnings
import os, sys
from Bio import SeqIO, SeqFeature, BiopythonWarning
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import argparse
import pyhmmer
import pandas as pd
import collections
from tqdm import tqdm

warnings.filterwarnings("ignore", category=UserWarning, module="Bio")
warnings.simplefilter('ignore', BiopythonWarning)


parser = argparse.ArgumentParser(description='Makes a genbank file out of (one scaffold of) a fasta and its proteins for use with clinker. If you do not specify a start and end ORF, it will default to the entire scaffold. Only works with prodigal output XOXO')

requiredNamed = parser.add_argument_group('required named arguments')

parser.add_argument('-protfile', metavar='[PRODIGAL PROTEIN FASTA]', help='FASTA protein sequences (PRODIGAL) for the provided nucfile.', default=None, required=True)

parser.add_argument('-prot_id', metavar='[PROTEIN ID]', help='ID of the protein of interest.', required=False)

parser.add_argument('-num_orfs', metavar='[NUMBER OF ORFS]', help='Number of ORFs upstream and downstream.', required=False, type=int)

parser.add_argument('-up', metavar='[NUMBER OF UPSTREAM ORFS]', help='Number of upstream ORFs.', required=False, type=int)

parser.add_argument('-down', metavar='[NUMBER OF DOWNSTREAM ORFS]', help='Number of downstream ORFs.', required=False, type=int)

requiredNamed.add_argument('-outfile', metavar='[OUTPUT FILENAME]', help='Output file')

parser.add_argument('-fasta', metavar='[FASTA FILENAME]', help='Optional output FASTA filename', default=None, required=False)

parser.add_argument('--all', action='store_true', help='Use the entire input file.')

parser.add_argument('-pfam', metavar='[PFAM OUTPUT FILE]', help='pfam_scan.pl output file', default=None, required=False)

parser.add_argument('-tsv', metavar='[TSV ANNOTATION FILE]', help='TSV file with HMM annotations', default=None, required=False)

parser.add_argument('-threads', metavar='[NUMBER OF THREADS]', help='Number of threads for PFAM annotation', default=1, type=int)




args = parser.parse_args()

protfile = os.path.abspath(args.protfile)
prot_id = args.prot_id  
num_orfs = args.num_orfs
num_orfs_up = args.up if args.up else num_orfs
num_orfs_down = args.down if args.down else num_orfs
outfile = args.outfile


def grab_scaffold_proteins(protein_recs, scaffold_id):
    # Retrieves proteins from desired scaffold
    def scaffold_filter(rec):
        return '_'.join(rec.name.decode().split('_')[:-1]) == scaffold_id

    return list(filter(scaffold_filter, protein_recs))



def grab_neighborhood(scaffold_proteins, prot_id, num_orfs_up, num_orfs_down):
    # find the index of the protein with the provided ID
    for i, rec in enumerate(scaffold_proteins):
        if rec.name.decode() == prot_id:
            target_index = i
            break
    else:
        raise ValueError("Protein ID not found: {}".format(prot_id))

    # find the start and end indices for the range of proteins we are interested in
    start = max(0, target_index - num_orfs_up)
    end = min(len(scaffold_proteins), target_index + num_orfs_down + 1)  # add 1 because range is exclusive

    # select the proteins in this range
    subset_recs = scaffold_proteins[start:end]

    return subset_recs


    return subset_recs
def parse_single_hmm(hmm_path):
    #Single-file parser for parallelization
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        return hmm_file.read()
    
def parse_hmms(hmm_in):
    #Checks first whether HMMs are provided as a single file or as a directory.

    hmms = []  # Initialize an empty list to store parsed HMMs
    print("Parsing HMMs...")
    # Check if hmm_in is a directory or a single file
    if os.path.isdir(hmm_in):
        if not os.listdir(hmm_in):
            print("hmm_in directory is empty.")
            logging.info('hmm_in directory is empty.')
            sys.exit(1)

        num_files = len(os.listdir(hmm_in))
        if num_files == 1:
            #Only one HMM file in input directory
            #Get full path to file
            hmm_path = os.path.join(hmm_in, os.listdir(hmm_in)[0])
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                #Works in case of single-model or multi-model HMM file
                hmms = list(hmm_file)

        else:
            hmm_files = list(filter(lambda x: x.endswith(('.hmm', '.HMM')), os.listdir(hmm_in)))
            hmm_paths = [os.path.join(hmm_in, hmm_file) for hmm_file in hmm_files]
            
            #I have tried!! Every possible method! To parallelize this!
            #It does not work. SINGLE THREADED IT IS!
            hmms = list(tqdm(map(parse_single_hmm, hmm_paths)))



    elif os.path.isfile(hmm_in):
        if os.path.getsize(hmm_in) == 0:
            print("hmm_in file is empty.")
            logging.info('hmm_in file is empty.')
            sys.exit(1)
        # Parse the single HMM file; handles multi-model files
        with pyhmmer.plan7.HMMFile(hmm_in) as hmm_file:
            hmms = list(hmm_file)
    else:
        print("Invalid HMM input.")
        logging.info("Invalid HMM input.")
        print("If you used pre-installed HMMs, check hmm_databases.json")
        logging.info("If you used pre-installed HMMs, check hmm_databases.json")
        print("Which is located in the databases directory.")
        logging.info("Which is located in the databases directory.")

        print("Thing that threw the error: {}".format(hmm_in))
        sys.exit(1)

    print("HMMs parsed.")

    return list(hmms)

def get_results_attributes(result):
    # turns named object into an ordered list
    # down with OOP
    bitscore = result.bitscore
    evalue = result.evalue
    cog = result.hmm_name
    c_evalue = result.c_evalue
    i_evalue = result.i_evalue
    query = result.sequence_id
    env_from = result.env_from
    env_to = result.env_to
    dom_bitscore = result.dom_bitscore
    return [query, cog, bitscore, evalue, c_evalue, i_evalue, env_from, env_to, bitscore]

def run_hmmsearch(hmms, pyhmmer_sequences, threads):

    results = []
    #Store as a global so we don't have to define it multiple times
    Result = collections.namedtuple("Result", ["sequence_id", "hmm_name", "bitscore", "evalue","c_evalue", "i_evalue", 
                                          "env_from", "env_to", "dom_bitscore"])
    # Run PyHMMER HMMsearch functionality
    for hits in pyhmmer.hmmsearch(hmms, pyhmmer_sequences, cpus=threads, bit_cutoffs='gathering'):
        # Get sequence ID
        cog = hits.query_name.decode()

        for hit in hits:
            if hit.included:
                # Get HMM name
                hit_name = hit.name.decode()
                full_bitscore = hit.score 
                full_evalue = hit.evalue
                # Domain specific hit information; generally very similar but can differ depending on sequence length
                for domain in hit.domains.reported:
                    results.append(Result(hit_name, cog, full_bitscore, full_evalue, domain.c_evalue, 
                            domain.i_evalue, domain.env_from, domain.env_to, domain.score))

    result_df = pd.DataFrame(list(map(get_results_attributes, results)), columns=["sequence_id", "hmm_name", "bitscore", "evalue","c_evalue", "i_evalue", "env_from", "env_to", "dom_bitscore"])
    return result_df

# Controller function; this is called from the main pipeline and calls the other helpers
def annotate_pfam(sequences, threads):
    """
    Given pyHMMER format sequences and number of threads, returns a dataframe of PFAM annotations
    """
    pfam_dir = '/home/jwestrob/.config/Astra/PFAM'
    hmms = parse_hmms(pfam_dir)
    result_df = run_hmmsearch(hmms, sequences, threads)
    return result_df

def parse_pfam_output(pfam_output_file):
    domains = dict()
    with open(pfam_output_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                columns = line.split()
                if len(columns) < 7:  # ensure that there are at least 7 columns
                    continue
                seq_id = columns[0]
                domain = columns[6]
                if seq_id in domains:
                    domains[seq_id].append(domain)
                else:
                    domains[seq_id] = [domain]
    return domains

def parse_tsv_annotations(tsv_file):
    annotations = dict()
    with open(tsv_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            columns = line.strip().split('\t')
            seq_id = columns[0]
            hmm_name = columns[1]
            if seq_id in annotations:
                annotations[seq_id].append(hmm_name)
            else:
                annotations[seq_id] = [hmm_name]
    
    # Join multiple annotations with ';'
    for seq_id, hmm_names in annotations.items():
        annotations[seq_id] = ';'.join(hmm_names)

    return annotations

def process_protein_record(protein_rec, total_start, pfam_annotations, scaf_rec_new):
    start = int(protein_rec.description.decode().split(' # ')[1])
    new_start = max(0, start - total_start)
    startpos = SeqFeature.ExactPosition(new_start)

    end = int(protein_rec.description.decode().split(' # ')[2])
    new_end = end - total_start
    endpos = int(SeqFeature.ExactPosition(new_end))

    strand = int(protein_rec.description.decode().split(' # ')[3])

    rec_location = FeatureLocation(startpos, endpos)
    rec_feature = SeqFeature.SeqFeature(rec_location, type="CDS", strand=strand)

    protein_annotations = pfam_annotations[pfam_annotations['sequence_id'] == protein_rec.name.decode()]
    if not protein_annotations.empty:
        hmm_names = protein_annotations['hmm_name'].unique()
        rec_feature.qualifiers['locus_tag'] = ';'.join(hmm_names)
    else:
        rec_feature.qualifiers['locus_tag'] = protein_rec.name.decode()

    rec_feature.qualifiers['translation'] = protein_rec.textize().sequence
    rec_feature.qualifiers['codon_start'] = 1
    rec_feature.qualifiers['transl_table'] = 11
    scaf_rec_new.features.append(rec_feature)


def main():
    protein_recs = list(SeqIO.parse(protfile, 'fasta'))

    # Create a list to store the digital sequences
    digital_sequences = []
    alphabet = pyhmmer.easel.Alphabet.amino()

    for rec in protein_recs:
        name_bytes = rec.id.encode('utf-8')
        seq_bytes = alphabet.encode(str(rec.seq))
        description_bytes = rec.description.encode('utf-8')
        sequence = pyhmmer.easel.DigitalSequence(
            name=name_bytes,
            sequence=seq_bytes,
            alphabet=alphabet,
            description=description_bytes
        )
        digital_sequences.append(sequence)

    if args.all:
        # Get all unique scaffold IDs
        all_scaffold_ids = set('_'.join(rec.id.split('_')[:-1]) for rec in protein_recs)

        for scaffold_id in all_scaffold_ids:
            scaffold_proteins = grab_scaffold_proteins(digital_sequences, scaffold_id)
            process_scaffold_proteins(scaffold_proteins, scaffold_id, args.threads)

    else:
        if prot_id is not None:
            scaffold_id = '_'.join(prot_id.split('_')[:-1])
            scaffold_proteins = grab_scaffold_proteins(digital_sequences, scaffold_id)

            if num_orfs_up is None or num_orfs_down is None:
                print("Warning: num_orfs_up and/or num_orfs_down are not provided. Using default values.")
                num_orfs_up = num_orfs_up or 5
                num_orfs_down = num_orfs_down or 5

            neighborhood_recs = grab_neighborhood(scaffold_proteins, prot_id, num_orfs_up, num_orfs_down)
            process_scaffold_proteins(neighborhood_recs, scaffold_id, args.threads)
        else:
            print("Warning: prot_id is not provided. Processing all scaffolds.")
            all_scaffold_ids = set('_'.join(rec.id.split('_')[:-1]) for rec in protein_recs)

            for scaffold_id in all_scaffold_ids:
                scaffold_proteins = grab_scaffold_proteins(digital_sequences, scaffold_id)
                process_scaffold_proteins(scaffold_proteins, scaffold_id, args.threads, args.pfam, args.tsv)

def process_scaffold_proteins(protein_recs, scaffold_id, threads):
    if args.fasta:
        fastafile = args.fasta
        with open(fastafile, 'w') as f:
            for rec in protein_recs:
                SeqIO.write(rec, f, 'fasta')

    total_start = int(protein_recs[0].description.decode().split(' # ')[1])
    total_end = int(protein_recs[-1].description.decode().split(' # ')[2])
    dummy_sequence = Seq('A' * (total_end - total_start), alphabet=generic_dna)
    scaf_rec_new = SeqRecord(dummy_sequence, id=scaffold_id)

    start_fp = int(SeqFeature.ExactPosition(total_start))
    end_fp = int(SeqFeature.ExactPosition(total_end))
    source_feature_location = FeatureLocation(start_fp, end_fp)
    source_feature = SeqFeature.SeqFeature(source_feature_location,
                                            type='source', strand=1)
    scaf_rec_new.features.append(source_feature)

    # Create a list to store the digital sequences for the neighborhood
    neighborhood_digital_sequences = []

    for rec in protein_recs:
        name_bytes = rec.name
        sequence = pyhmmer.easel.DigitalSequence(
            name=name_bytes,
            sequence=rec.sequence,
            alphabet=rec.alphabet
        )
        neighborhood_digital_sequences.append(sequence)

    # Run PFAM annotations on the neighborhood sequences
    annotations = annotate_pfam(neighborhood_digital_sequences, threads)

    for protein_rec in protein_recs:
        process_protein_record(protein_rec, total_start, annotations, scaf_rec_new)

    scaf_rec_new.annotations['molecule_type'] = 'DNA'
    print("Writing neighborhood genbank to " + outfile)
    print('---')
    SeqIO.write(scaf_rec_new, outfile, 'genbank')

if __name__ == "__main__":
    main()

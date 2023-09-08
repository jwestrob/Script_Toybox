import warnings
import os, sys
from Bio import SeqIO, SeqFeature, BiopythonWarning
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import argparse

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

parser.add_argument('-pfam', metavar='[PFAM OUTPUT FILE]', help='pfam_scan.pl output file', default=None, required=True)


args = parser.parse_args()

protfile = os.path.abspath(args.protfile)
prot_id = args.prot_id  
num_orfs = args.num_orfs
num_orfs_up = args.up if args.up else num_orfs
num_orfs_down = args.down if args.down else num_orfs
outfile = args.outfile

def grab_neighborhood(scaffold_proteins, prot_id, num_orfs_up, num_orfs_down):
    # find the index of the protein with the provided ID
    for i, rec in enumerate(scaffold_proteins):
        if rec.id == prot_id:
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

def process_protein_record(protein_rec, total_start, domains, scaf_rec_new):
    start = int(protein_rec.description.split(' # ')[1])
    new_start = max(0, start - total_start)
    startpos = SeqFeature.ExactPosition(new_start)

    end = int(protein_rec.description.split(' # ')[2])
    new_end = end - total_start
    endpos = int(SeqFeature.ExactPosition(new_end))

    strand = int(protein_rec.description.split(' # ')[3])

    rec_location = FeatureLocation(startpos, endpos)
    rec_feature = SeqFeature.SeqFeature(rec_location, type="CDS", strand=strand)
    rec_feature.qualifiers['protein_id'] = protein_rec.id
    rec_feature.qualifiers['translation'] = protein_rec.seq

    if protein_rec.id in domains:
        rec_feature.qualifiers['protein_id'] = ",".join(domains[protein_rec.id])
    else:
        rec_feature.qualifiers['protein_id'] = protein_rec.id

    rec_feature.qualifiers['translation'] = str(protein_rec.seq)
    rec_feature.qualifiers['codon_start'] = 1
    rec_feature.qualifiers['transl_table'] = 11
    scaf_rec_new.features.append(rec_feature)

def grab_scaffold_proteins(protein_recs, scaffold_id):
    #Retrieves proteins from desired scaffold
    def scaffold_filter(rec):
        return '_'.join(rec.id.split('_')[:-1]) == scaffold_id

    return list(filter(scaffold_filter, protein_recs))

def main():
    protein_recs = list(SeqIO.parse(protfile, 'fasta'))

    if args.all:
        # Get all unique scaffold IDs
        all_scaffold_ids = set('_'.join(rec.id.split('_')[:-1]) for rec in protein_recs)

        for scaffold_id in all_scaffold_ids:
            scaffold_proteins = grab_scaffold_proteins(protein_recs, scaffold_id)
            process_scaffold_proteins(scaffold_proteins, scaffold_id)  # Pass scaffold_id here

    else:
        scaffold_id = '_'.join(prot_id.split('_')[:-1])
        scaffold_proteins = grab_scaffold_proteins(protein_recs, scaffold_id)
        neighborhood_recs = grab_neighborhood(scaffold_proteins, prot_id, num_orfs_up, num_orfs_down)
        process_scaffold_proteins(neighborhood_recs, scaffold_id)  # Pass scaffold_id here


def process_scaffold_proteins(protein_recs, scaffold_id):
    if args.fasta:
        fastafile = args.fasta
        with open(fastafile, 'w') as f:
            for rec in protein_recs:
                SeqIO.write(rec, f, 'fasta')

    if args.pfam is not None:
        pfam_output_file = args.pfam
        domains = parse_pfam_output(pfam_output_file)
    else:
        domains = {}

    total_start = int(protein_recs[0].description.split(' # ')[1])
    total_end = int(protein_recs[-1].description.split(' # ')[2])
    dummy_sequence = 'A' * (total_end - total_start)
    scaf_rec_new = SeqRecord(Seq(dummy_sequence), id=scaffold_id)

    start_fp = int(SeqFeature.ExactPosition(total_start))
    end_fp = int(SeqFeature.ExactPosition(total_end))
    source_feature_location = FeatureLocation(start_fp, end_fp)
    source_feature = SeqFeature.SeqFeature(source_feature_location,
                                            type='source', strand=1)
    scaf_rec_new.features.append(source_feature)

    for protein_rec in protein_recs:
        process_protein_record(protein_rec, total_start, domains, scaf_rec_new)

    scaf_rec_new.annotations['molecule_type'] = 'DNA'
    print("Writing neighborhood genbank to " + outfile)
    print('---')
    SeqIO.write(scaf_rec_new, outfile, 'genbank')



if __name__ == "__main__":
    main()

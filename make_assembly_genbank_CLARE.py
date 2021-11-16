from Bio import SeqIO, SeqFeature
#from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqFeature import FeatureLocation
import os, sys, pandas as pd
import time
from pathos.multiprocessing import ProcessingPool as Pool
import argparse


parser = argparse.ArgumentParser(
    description='Makes an enormous genbank file out of an assembly fasta and its proteins.')

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-nucfile', metavar='[NUCLEOTIDE FASTA FILE]',
                    help='(Metagenomic) fasta sequence file.', required=True)

requiredNamed.add_argument('-protfile', metavar='[PRODIGAL PROTEIN FASTA]',
                    help='FASTA format genes for the provided nucfile.', default=None, required=True)

requiredNamed.add_argument('-threads', metavar='[NUM THREADS]',
                    help='Number of threads to use. Give me lots please.', default=1, required=True)

requiredNamed.add_argument('-outfile', metavar='[OUTPUT FILENAME]', help='Output file')

args = parser.parse_args()

assembly_nuc = os.path.abspath(args.nucfile)
assembly_prot = os.path.abspath(args.protfile)
outfile = args.outfile
threads = int(args.threads)

#print(assembly_nuc, assembly_prot, outfile, threads)
#sys.exit()

t1 = time.time()
new_recs = []

protein_recs = list(SeqIO.parse(assembly_prot, 'fasta'))


def make_genbank_recs(rec):
    new_rec = rec
    #new_rec.seq.alphabet = generic_dna
    scaffold = new_rec.id
    
    scaffold_recs  = list(filter(lambda x: x.id.startswith(scaffold + '_'), protein_recs))
    
    
    for protein_rec in scaffold_recs:
        start = int(protein_rec.description.split(' # ')[1])
        startpos = SeqFeature.ExactPosition(start)
        end = int(protein_rec.description.split(' # ')[2])
        endpos = int(SeqFeature.ExactPosition(end))
        strand = int(protein_rec.description.split(' # ')[3])
        rec_location = FeatureLocation(startpos, endpos)
        rec_feature = SeqFeature.SeqFeature(rec_location, type="CDS", strand=strand)

        #Add ORF name without genome ID
        rec_feature.qualifiers['protein_id'] = protein_rec.id
        rec_feature.qualifiers['translation'] = protein_rec.seq
        rec_feature.qualifiers['locus_tag'] = protein_rec.description

        new_rec.features.append(rec_feature)
    return new_rec
    
p = Pool(20)

genbank_recs = list(p.map(make_genbank_recs, SeqIO.parse(assembly_nuc, 'fasta')))

if not outfile.endswith('.gbk'):
    outfile = outfile + '.gbk'

SeqIO.write(genbank_recs, outfile, 'genbank')
    
print("process took " + str(time.time() - t1) + " seconds.")   

import os, sys
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import argparse

parser = argparse.ArgumentParser(description='Makes a genbank file out of (one scaffold of) a fasta and its proteins for use with clinker. If you do not specify a start and end ORF, it will default to the entire scaffold. Only works with prodigal output XOXO')

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-protfile', metavar='[PRODIGAL PROTEIN FASTA]', help='FASTA protein sequences (PRODIGAL) for the provided nucfile.', default=None, required=True)

requiredNamed.add_argument('-prot_id', metavar='[PROTEIN ID]', help='ID of the protein of interest.', required=True)

requiredNamed.add_argument('-num_orfs', metavar='[NUMBER OF ORFS]', help='Number of ORFs upstream and downstream.', required=True, type=int)

requiredNamed.add_argument('-up', metavar='[NUMBER OF UPSTREAM ORFS]', help='Number of upstream ORFs.', required=False, type=int)

requiredNamed.add_argument('-down', metavar='[NUMBER OF DOWNSTREAM ORFS]', help='Number of downstream ORFs.', required=False, type=int)

requiredNamed.add_argument('-outfile', metavar='[OUTPUT FILENAME]', help='Output file')

parser.add_argument('-fasta', metavar='[FASTA FILENAME]', help='Optional output FASTA filename', default=None, required=False)


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
        raise ValueError("Protein ID not found.")

    # find the start and end indices for the range of proteins we are interested in
    start = max(0, target_index - num_orfs_up)
    end = min(len(scaffold_proteins), target_index + num_orfs_down + 1)  # add 1 because range is exclusive

    # select the proteins in this range
    subset_recs = scaffold_proteins[start:end]

    return subset_recs

def grab_scaffold_proteins(protein_recs, scaffold_id):
    #Retrieves proteins from desired scaffold
    def scaffold_filter(rec):
        return '_'.join(rec.id.split('_')[:-1]) == scaffold_id

    return list(filter(scaffold_filter, protein_recs))

def main():
    scaffold_id = '_'.join(prot_id.split('_')[:-1])

    protein_recs = list(SeqIO.parse(protfile, 'fasta'))

    scaffold_proteins = grab_scaffold_proteins(protein_recs, scaffold_id)

    neighborhood_recs = grab_neighborhood(scaffold_proteins, prot_id, num_orfs_up, num_orfs_down)

    if args.fasta:
        fastafile = args.fasta
        with open(fastafile, 'w') as f:
            for rec in neighborhood_recs:
                SeqIO.write(rec, f, 'fasta')

    #Retrieve coordinates of the start of the first and last ORFs to subset the scaffold sequence
    #Retrieves the start of the first ORF
    total_start = int(neighborhood_recs[0].description.split(' # ')[1])
    #Retrieves the end of the last ORF
    total_end = int(neighborhood_recs[-1].description.split(' # ')[2])

    #Create a dummy sequence of 'A's with the length equal to total_end - total_start
    dummy_sequence = 'A' * (total_end - total_start)
    scaf_rec_new = SeqRecord(Seq(dummy_sequence), id=scaffold_id)


    #Printouts in case something goes wrong and the coordinates are formatted improperly
    if len(scaf_rec_new.seq) == 0:
        print(len(scaf_rec))
        print(total_start)
        print(total_end)
        print("Coordinates formatted improperly! This is probably Jacob's fault, talk to him!")
        sys.exit(420)

    #Genbanks require a 'source' feature at the start. It does nothing.
    #But clinker doesn't work if you don't have it sooooooooo...
    start_fp = int(SeqFeature.ExactPosition(total_start))
    end_fp = int(SeqFeature.ExactPosition(total_end))

    source_feature_location = FeatureLocation(start_fp, end_fp)
    source_feature = SeqFeature.SeqFeature(source_feature_location,
                                            type='source', strand=1)

    scaf_rec_new.features.append(source_feature)

    #For initial debugging purposes; to be removed
    print("Total start: ", total_start)
    print("Total end: ", total_end)

    for protein_rec in neighborhood_recs:
        #Grab coordinate at the start of the ORF and format as a SeqFeature.ExactPosition
        start = int(protein_rec.description.split(' # ')[1])
        new_start = max(0, start - total_start)
        startpos = SeqFeature.ExactPosition(new_start)
        print("Start: ", start)
        print("New start: ", new_start)

        #Grab coordinate of the end of the ORF and format as a SeqFeature.ExactPosition
        end = int(protein_rec.description.split(' # ')[2])
        new_end = end - total_start
        endpos = int(SeqFeature.ExactPosition(new_end))
        print("End: ", end)
        print("New end: ", new_end)

        #Grab strand info for clinker display
        strand = int(protein_rec.description.split(' # ')[3])

        #Declare full FeatureLocation object
        rec_location = FeatureLocation(startpos, endpos)
        #Initialize feature
        rec_feature = SeqFeature.SeqFeature(rec_location, type="CDS",
                                            strand=strand)

        #I don't think this is necessary but if it doesn't work later I'll add it back in
        """
        location_feature = SeqFeature.SeqFeature(rec_location, type='gene',
                                                strand=strand)
        location_feature.qualifiers['locus_tag'] = protein_rec.id
        """

        #Add ORF name without genome ID
        rec_feature.qualifiers['protein_id'] = protein_rec.id
        rec_feature.qualifiers['translation'] = protein_rec.seq
        """
        If you have annotation information, modify this and add that information
        into the 'locus_tag' field, and you'll be able to see it in clinker.
        """
        rec_feature.qualifiers['locus_tag'] = protein_rec.id
        #Add the protein sequence; key for clinker alignments!
        rec_feature.qualifiers['translation'] = str(protein_rec.seq)
        rec_feature.qualifiers['codon_start'] = 1
        rec_feature.qualifiers['transl_table'] = 11
        scaf_rec_new.features.append(rec_feature)
        print('----------------------------')
    #Another useless bit that is required for no reason whatsoever!
    #If you don't have it SeqIO.write will fail
    scaf_rec_new.annotations['molecule_type'] = 'sleihgsldkjfskjfdngslkdjfaerkj'

    print("Writing neighborhood genbank to " + outfile)
    SeqIO.write(scaf_rec_new, outfile, 'genbank')


if __name__ == "__main__":
    main()

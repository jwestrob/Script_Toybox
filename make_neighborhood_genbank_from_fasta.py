import os, sys
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import argparse

parser = argparse.ArgumentParser(
    description='Makes a genbank file out of (one scaffold of) a fasta and its proteins for use with clinker.\
                 If you do not specify a start and end ORF, it will default to the entire scaffold. Only works with prodigal output XOXO')

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-nucfile', metavar='[NUCLEOTIDE FASTA FILE]',
                    help='Fasta nucleotide sequence file.', required=True)

requiredNamed.add_argument('-protfile', metavar='[PRODIGAL PROTEIN FASTA]',
                    help='FASTA protein sequences (PRODIGAL) for the provided nucfile.', default=None, required=True)

requiredNamed.add_argument('-scaffold_id', metavar='[SCAFFOLD ID]',
                    help='ID of the scaffold you want to make a genbank out of.')

requiredNamed.add_argument('-start_orfnum', metavar='[NEIGHBORHOOD_START_ORF]',
                    help='Number of the ORF that begins the neighborhood of interest. Default=1', default=1)

requiredNamed.add_argument('-end_orfnum', metavar='[NEIGHBORHOOD_END_ORF]',
                    help='Number of the ORF that ends the neighborhood of interest. Default=Last ORF.', default=-1)

requiredNamed.add_argument('-outfile', metavar='[OUTPUT FILENAME]', help='Output file')

args = parser.parse_args()

nucfile = os.path.abspath(args.nucfile)
protfile = os.path.abspath(args.protfile)
scaffold_id = args.scaffold_id
start_orfnum = int(args.start_orfnum)
end_orfnum = int(args.end_orfnum)

outfile = args.outfile

def grab_rec(nucfile, scaffold_id=scaffold_id):
    ### Parses nucleotide infile and retrieves scaffold of interest
    for rec in SeqIO.parse(nucfile, 'fasta'):
        if rec.id == scaffold_id:
            return rec
    #Were we unable to find the scaffold ID specified?
    print("Scaffold ID not found in nucfile. Please check formatting.")
    sys.exit(420)

def grab_neighborhood(scaffold_proteins, start=start_orfnum, end=end_orfnum):
    ### Subsets the protein recs to only include those from the neighborhood of interest.

    #Base case; you want the whole scaffold
    if start == 1 and end == -1:
        return protein_recs


    def neighborhood_filter(rec):
        #Retrieves proteins in desired neighborhood
        orfnum = int(rec.id.split('_')[-1])
        if orfnum >= start and orfnum <= end:
            return True
        return False

    subset_recs = list(filter(neighborhood_filter, scaffold_proteins))
    return subset_recs

def grab_scaffold_proteins(protein_recs):
    #Retrieves proteins from desired scaffold
    def scaffold_filter(rec):
        if '_'.join(rec.id.split('_')[:-1]).endswith(scaffold_id):
            return True
        return False

    return list(filter(scaffold_filter, protein_recs))


def main():

    scaffold_rec = grab_rec(nucfile)

    protein_recs = list(SeqIO.parse(protfile, 'fasta'))

    scaffold_proteins = grab_scaffold_proteins(protein_recs)

    neighborhood_recs = grab_neighborhood(scaffold_proteins)

    #Retrieve coordinates of the start of the first and last ORFs to subset the scaffold sequence
    #Retrieves the start of the first ORF
    total_start = int(neighborhood_recs[0].description.split(' # ')[1])
    #Retrieves the end of the last ORF
    total_end = int(neighborhood_recs[-1].description.split(' # ')[2])

    scaf_rec_new = SeqRecord(Seq(str(scaffold_rec.seq[total_start:total_end])), id=scaffold_rec.id)

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

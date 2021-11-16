import os, sys
from Bio.Alphabet import generic_dna, generic_protein
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

if '-h' in sys.argv:
    print("Example usage: get_neighborhood.py [center_orf_id] [neighborhood_size] [genbank_file] [output_genbank_filename]")
    print("This is made for use with ggkbase-formatted genbank files. Please don't use normal genbank files with it. XOXO")
    sys.exit()

#When you run the script, give it the name of the ORF you want to center around
desired_orfname = sys.argv[1]

#How many ORFs above and below the central ORF do you want to keep?
neighborhood_size = int(sys.argv[2])

#Load in the genbank file
recs = list(SeqIO.parse(sys.argv[3], 'genbank'))
if len(recs) > 1:
    print("There is more than one scaffold in your genbank file. Please fix that.")
    print("I am way too lazy to write a scaffold filter in this script")
    print("Why do you do this to me")
    sys.exit()

output_filename = sys.argv[4]

new_recs = []

#First feature is 'source' which contains no useful info. Delete
features = recs[0].features[1:]

prev_gene = None
new_features = []
for feature in features:
    if feature.type == 'tRNA':
        continue
    if feature.type == 'gene':
        #All this just to grab the ORF name... ugh
        prev_gene = feature.qualifiers['locus_tag'][0]
    if feature.type == 'CDS':
        feature.qualifiers['locus_tag'] = [prev_gene]
        new_features.append(feature)


new_rec = recs[0]
new_rec.features = new_features

#OK, now that the ggkbase formatting is fixed, let's move on...

desired_orfnum = int(desired_orfname.split('_')[-1])

def in_neighborhood_filter(feature, desired_orfnum, neighborhood_size):
    if feature.type == 'gene' or feature.type == 'tRNA':
        return False
    feature_orfname = feature.qualifiers['locus_tag'][0]

    try:
        feature_orfnum = int(''.join(feature_orfname.split('_')[-1].split()))
    except:
        print("ERROR")
        print(feature_orfname)
        sys.exit()

    #Is it in the neighborhood?
    if abs(feature_orfnum - desired_orfnum) <= neighborhood_size:
        return True
    else:
        return False


in_neighborhood_features = list(filter(lambda x: in_neighborhood_filter(x, desired_orfnum, neighborhood_size), features))

total_start = in_neighborhood_features[0].location.start
total_end = in_neighborhood_features[-1].location.end

#Declare new SeqRecord; subset the nucleotide sequence for visualization in clinker
neighborhood_rec = SeqRecord(
                    Seq(str(recs[0].seq[total_start:total_end]), alphabet=generic_dna)
                    )
neighborhood_rec.features = in_neighborhood_features



new_recs = []
for nrec in [neighborhood_rec]:
	new_rec = nrec
	new_features = []
	for feature in in_neighborhood_features:
		cur_strand = feature.location.strand
		orig_start = feature.location.start
		orig_end = feature.location.end

		new_start = max(0, orig_start - total_start)
		new_end = orig_end - total_start

		new_feature = feature
		new_feature.location = FeatureLocation(new_start, new_end, cur_strand)
		new_features.append(new_feature)

	new_rec.features = new_features
	new_recs.append(new_rec)

SeqIO.write(new_recs, output_filename, 'gb')

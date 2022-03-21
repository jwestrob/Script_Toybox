import os, sys
#from Bio.Alphabet import generic_dna, generic_protein
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

#Let's get the scaffold ID from the orfname.
if '|' in desired_orfname:
    #Check for extra weird ggkbase orf format
    if len(desired_orfname.split('|')) == 3:
        scaffold_id = '_'.join(desired_orfname.split('|')[0].split('_')[:-1])
    else:
        scaffold_id = '_'.join(desired_orfname.split('|')[1].split('_')[:-1])
else:
    scaffold_id = '_'.join(desired_orfname.split('_')[:-1])

print("Scaffold ID: ", scaffold_id)

#How many ORFs above and below the central ORF do you want to keep?
neighborhood_size = int(sys.argv[2])


#Let's run a preliminary check to fix the genbank file because
#ggkbase spits out poorly formatted dates and biopython
#was written by assholes who kill your script if the dates
#are formatted incorrectly. Pinche carajo

genbank_file = sys.argv[3]

with open(genbank_file, 'r') as infile:
    gb_lines = [x.rstrip() for x in infile.readlines()]

def fix_brians_dates(line):
    if line.startswith('LOCUS') and 'BCT' in line:

        split_line = line.split('BCT ')
        line_no_date = split_line[0]
        return line_no_date + 'BCT 01-JAN-1980'
    else:
        return line

fixed_lines = list(map(fix_brians_dates, gb_lines))

with open(genbank_file, 'w') as outfile:
    for element in fixed_lines:
        outfile.writelines(element + '\n')

#Load in the genbank file
recs = list(SeqIO.parse(sys.argv[3], 'genbank'))

if len(recs) > 1:

    try:
        scaffold_rec = list(filter(lambda x: x.id == scaffold_id, recs))[0]
    except:

        ids  = [rec.id for rec in recs]
        #check for NCBI format bullshit
        if ids[0].endswith('.1'):
            scaffold_rec = list(filter(lambda x: x.id == scaffold_id + '.1', recs))[0]
        else:
            print("Scaffold not found; please check the format of your inputs to ensure that they match.")
            sys.exit()

else:
    scaffold_rec = recs[0]

output_filename = sys.argv[4]

new_recs = []

#First feature is 'source' which contains no useful info. Delete
features = scaffold_rec.features[1:]

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
    try:
        feature_orfname = feature.qualifiers['locus_tag'][0]
    except:
        #Then it's a 16S
        #print(feature.qualifiers)
        return False

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


try:
	total_start = in_neighborhood_features[0].location.start
except:
	print("Num features:")
	print(len(in_neighborhood_features))

total_end = in_neighborhood_features[-1].location.end

#Declare new SeqRecord; subset the nucleotide sequence for visualization in clinker
neighborhood_rec = SeqRecord(
                    Seq(str(recs[0].seq[total_start:total_end])), annotations={"molecule_type": "DNA"}
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

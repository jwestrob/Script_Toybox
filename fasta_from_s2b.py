import os, sys, pandas as pd
from Bio import SeqIO

if '-h' in sys.argv:
	print('Example usage: python split_s2b.py [SCAF2BIN] [OUTDIR] [CONTIGFILE]')
#	print('Optional flag (put at the end): ggkbase')
	sys.exit()

if len(sys.argv) < 3:
	print('give me all the arguments man')
	print('run -h')
	sys.exit()


s2b = sys.argv[1]

outdir = sys.argv[2]

contigfile = sys.argv[3]


s2b_df = pd.read_csv(s2b, sep='\t', header=None)

#If pandas parses header stupidly
if s2b_df.iloc[0].values[0].lower().startswith('scaffold'):
	s2b_df = s2b_df.iloc[1:]

if len(s2b_df.columns) == 3:
	s2b_df.columns = ['scaffold', 'bin', 'taxonomy']

else:
	s2b_df.columns = ['scaffold', 'bin']

s2b_df_nounk = s2b_df[~s2b_df.bin.str.endswith('_UNK')]

bins = s2b_df_nounk.bin.unique().tolist()

bindict = dict(zip(s2b_df_nounk.scaffold.tolist(), s2b_df_nounk.bin.tolist()))

bin_list_dict = dict(zip(bins, [[] for bin in bins]))

for rec in SeqIO.parse(contigfile, 'fasta'):
	if rec.id in bindict:
		dest_bin = bindict[rec.id]
		bin_list_dict[dest_bin].append(rec)

if not os.path.exists(outdir):
	os.system('mkdir ' + outdir)

for bin in bin_list_dict.keys():
	SeqIO.write(bin_list_dict[bin], os.path.join(outdir, bin + '.fna'), 'fasta')

print('Tada')

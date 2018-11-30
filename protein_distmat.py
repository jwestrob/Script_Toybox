import os, sys, argparse, pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import pylab

def parse():
    parser = argparse.ArgumentParser(description='Take a set of protein sequences (fasta) and yield distance matrix.'+'\n' + \
                                    'More options to be considered soon; Feel free to suggest some.')

    parser.add_argument('-i', metavar='[infile]', help='File of protein input sequences.')
    parser.add_argument('-o', metavar='[outfile]', help='File to write distance matrix to (.npy format)')
    parser.add_argument('-outdir', metavar='[output directory]', default=None, \
                        help='Directory for blastp output/db. Default: Working directory')
    parser.add_argument('-tsne', action='store_true', default=False, help='Want to run a t-SNE on it?')
    parser.add_argument('-t', metavar='[threads]', default=1, help='Number of threads to use.')

    args = parser.parse_args()

    return args



def main():

    args = parse()

    infile = args.i
    outfile = str(args.o)
    outdir = str(args.outdir)
    t = int(args.t)
    tsne = str(args.tsne)

    if outdir is not None:
        if not os.path.exists(outdir):
            os.system('mkdir ' + outdir)
        wd = os.getcwd()
        infile_path = wd + '/' + infile

        #Move to directory for output
        os.chdir(outdir)
        #Make symbolic link to outfile
        os.system('ln -s ' + infile_path + ' .')
    else:
        outdir = os.getcwd()
    #Make BLASTP database
    os.system('makeblastdb -in ' + infile + ' -out ' + infile.split('.')[0] + ' -dbtype prot')

    print('blastp -db ' + infile.split('.')[0] + ' -query ' + infile + \
                ' -outfmt 6 -out ' + outdir + '/all-vs-all.tsv -num_threads ' + str(t))

    #Run BLASTP all-vs-all
    os.system('blastp -db ' + infile.split('.')[0] + ' -query ' + infile + \
                ' -outfmt 6 -out all-vs-all.tsv -num_threads ' + str(t))

    if outdir is not None:
        #Run along home
        os.chdir(wd)

    #Load in as Pandas DF
    allvall = pd.read_csv(outdir + '/all-vs-all.tsv', sep='\t', low_memory=False, header=None)

    #Leave out unused columns (we need Query, Hit, %ID)
    red_df = allvall[[0, 1, 2]]

    #Convert to numpy array
    red_arr = red_df.as_matrix()

    #Extract list of unique IDs from first column (Series)
    unique_entries = red_df[0].unique().tolist()

    #Get number of total unique IDs
    num_entries = len(unique_entries)

    #Extract X or Y coordinate on resulting distance matrix from hit ID
    name_dict = dict(zip(unique_entries, range(num_entries)))

    #Create empty dismat
    distmat = np.zeros((num_entries, num_entries))

    #Put in your distances
    for index, row in red_df.iterrows():
        row = red_df.iloc[index]
        x = name_dict[row[0]]
        y = name_dict[row[1]]
        dist = 1.0 - row[2]/100.0
        distmat[x][y] = dist
        distmat[y][x] = dist

    #Save distmat as .npy object
    np.save(outdir + '/' + outfile, distmat)

    np.savetxt(outdir + '/unique_ids.txt', np.array(unique_entries), delimiter='\t', fmt='%s')

    if tsne:
        two_embedding = TSNE(n_components=2, n_iter=2000, perplexity=50).fit_transform(distmat)
        pylab.scatter(two_embedding[:,0], two_embedding[:,1])
        pylab.title("2D t-SNE (protein)")
        pylab.show()

    os.system('Rscript /home/jacob/scripts/dale.R')
if __name__ == "__main__":
    main()

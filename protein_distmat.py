import os, sys, argparse, pandas as pd
import time
from pathos.multiprocessing import ProcessingPool as Pool
#from sklearn.manifold import TSNE
from MulticoreTSNE import MulticoreTSNE as TSNE
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
    parser.add_argument('-just_blast', action='store_true', default=False, help='Just run BLAST and make all-vs-all.tsv.')
    parser.add_argument('-already_blasted', action='store_true', default=False, help='Did you already run BLAST and then pandas yelled something stupid at you and dumped? use this flag')
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
    just_blast = args.just_blast
    blasted = args.already_blasted

    if outdir is not None:
        if not os.path.exists(outdir):
            os.system('mkdir ' + outdir)
        wd = os.getcwd()
        infile_path = wd + '/' + infile

        #Move to directory for output
        os.chdir(outdir)

        if not os.path.exists(os.path.join(outdir, infile)):
            #Make symbolic link to outfile
            os.system('ln -s ' + infile_path + ' .')
    else:
        outdir = os.getcwd()
    if not blasted:
        #Make BLASTP database
        os.system('makeblastdb -in ' + infile + ' -out ' + infile.split('.')[0] + ' -dbtype prot')

        print('blastp -db ' + infile.split('.')[0] + ' -query ' + infile + \
                ' -outfmt 6 -out ' + outdir + '/all-vs-all.tsv -num_threads ' + str(t))

        #Run BLASTP all-vs-all
        os.system('blastp -db ' + infile.split('.')[0] + ' -query ' + infile + \
                ' -outfmt 6 -out all-vs-all.tsv -num_threads ' + str(t))

    if just_blast:
        sys.exit()

    if outdir is not None:
        #Run along home
        os.chdir(wd)

    #Load in as Pandas DF
    allvall = pd.read_csv(outdir + '/all-vs-all.tsv', sep='\t', low_memory=False, header=None)

    #Leave out unused columns (we need Query, Hit, %ID)
    red_df = allvall[[0, 1, 2]]
    red_df.columns = ['seq1', 'seq2', 'pident']

    #Extract list of unique IDs from first column (Series)
    unique_entries = red_df['seq1'].unique().tolist()

    #Get number of total unique IDs
    num_entries = len(unique_entries)

    #Extract X or Y coordinate on resulting distance matrix from hit ID
    name_dict = dict(zip(unique_entries, range(num_entries)))

    print("Calculating distances...")

    t1 = time.time()
    red_df.seq1 = red_df.seq1.apply(lambda x: int(name_dict[x]))
    red_df.seq2 = red_df.seq2.apply(lambda x: int(name_dict[x]))
    red_df['dist'] = red_df.pident.apply(lambda x: 1.0 - (float(x) / 100.0))
    print("Distance calculation took " + str(time.time() - t1) + " seconds.")

    #Create empty dismat
    distmat = np.zeros((num_entries, num_entries))



    print("Distances calculated, writing to distmat...")

    def etch(distlist):
        idx1 = int(distlist[0])
        idx2 = int(distlist[1])
        dist = distlist[2]
        try:
            distmat[idx1, idx2] = dist
        except:
            print(idx1, idx2, dist)
            sys.exit()
        distmat[idx2, idx1] = dist
        return

    p = Pool(t)

    t1 = time.time()
    list(map(etch, red_df[['seq1', 'seq2', 'dist']].values  ))
    print("Etch took " + str(time.time() - t1) + " seconds.")

    """
    #Put in your distances
    for index, row in red_df.iterrows():
        row = red_df.iloc[index]
        x = name_dict[row[0]]
        y = name_dict[row[1]]
        dist = 1.0 - row[2]/100.0
        distmat[x][y] = dist
        distmat[y][x] = dist
    """

    #Save distmat as .npy object
    np.save(outdir + '/' + outfile, distmat)

    np.savetxt(outdir + '/unique_ids.txt', np.array(unique_entries), delimiter='\t', fmt='%s')

    if tsne:
        two_embedding = TSNE(n_jobs=t, n_components=2, n_iter=2000, perplexity=50).fit_transform(distmat)
        np.save(outdir + '/tsne_embedding_matrix.npy', two_embedding)
        pylab.scatter(two_embedding[:,0], two_embedding[:,1])
        pylab.title("2D t-SNE (protein)")
        #pylab.show()
        plt.savefig(outdir + '/tsne_embedding.png')

    # 305
    os.system('Rscript /home/jacob/scripts/dale.R')
if __name__ == "__main__":
    main()

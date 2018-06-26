from multiprocessing.dummy import Pool as ThreadPool
import numpy as np
import pandas as pd
import os, sys, csv, time
import argparse

parser = argparse.ArgumentParser(description="Take a list of fastq read IDs or contig IDs of length n, then a list containing a subset of those IDs, of length m. Returns a binary presence-absence vector of length n.")

parser.add_argument('-r', metavar='list of read IDs (one per line)', nargs='?', default=None, help="Read IDs. no > characters. One per line. Generate with: grep '^>' example.fasta > read_ids.fasta.")
parser.add_argument('-s', metavar='subset of read IDs (one per line)', help="Subset read IDs. Same conditions as above.")
parser.add_argument('-o', metavar='Name of .npy outfile', help='Name yo file')
parser.add_argument('-t', metavar='Number of threads to use', default=1, type=int, help="GIMME AN INT!")

args = parser.parse_args()
allreads = str(args.r)
subset = str(args.s)
outfile = str(args.o)
threads = int(args.t)

t1 = time.time()

all_read_IDs = np.loadtxt(str(args.r), dtype=str)

#Create list comprehension
all_read_IDs = [[i] for i in all_read_IDs]

subset_read_IDs = np.loadtxt(str(args.s), dtype=str)

#Add indices to all_read_IDs
for index, read in enumerate(all_read_IDs):
    all_read_IDs[index].append(index)

#Define label vector as all zeros
labels = np.zeros(len(all_read_IDs), dtype=np.int8)

#Define function (using global variables to hopefully cheat)

def check_for_membership(comp_list):
    #Takes in a list from a list comprehension of the form [Read_ID, index]
    if comp_list[0] in subset_read_IDs:
        labels[comp_list[1]] = 1
    return

pool = ThreadPool(threads)
results = pool.map(check_for_membership, all_read_IDs)

#Save to .npy file
np.savetxt(outfile, labels, fmt='%5.0f')

print("Time elapsed: " + str(time.time() - t1) + " using " + str(threads) + " threads.")

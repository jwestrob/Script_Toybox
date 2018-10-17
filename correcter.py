import os, sys, math
import numpy as np
import time
from Bio import SeqIO

def trimmer(filename):
    found = False
    new_fasta = []
    with open(filename) as inputfile:
        for line in inputfile:
            if line[0] == ">":
                found = True
            if found == True:
                new_fasta.append(line)
    return new_fasta

def writer(filename, new_fasta):
    with open(filename, 'wb') as ffile:
        for item in new_fasta:
            print>>ffile, item
    return

def quick_convert(fastaname):
    SeqIO.convert(fastaname, "fasta", 'corrected_' + fastaname, "fasta")
    os.system('mv corrected_' + fastaname + ' ' + fastaname)
    return


def main():
    t1 = time.time()
    for filename in os.listdir('.'):
        if filename.endswith(".fasta"):
            filename_arr = filename.split(".")
            new_fasta = trimmer(filename)
            quick_convert(filename)
            #writer(filename_arr[0], new_fasta)
            #os.system('rm ' + filename)
            #os.system('mv ' + filename_arr[0] + ' ' + filename)
    print "Process took " + str(time.time()-t1) + " seconds."

if __name__ == "__main__":
    main()

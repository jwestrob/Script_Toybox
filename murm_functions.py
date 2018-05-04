import pandas as pd
import numpy as np
import os, sys, csv
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Reduce a .csv of PubMLST metadata based on the files in a given directory.')

parser.add_argument('-p', metavar='--pubmlst', help='Path to PubMLST info .csv file')
parser.add_argument('-dir', metavar='--fastadir', help='Path to directory containing your favorite FASTA files.')

args = parser.parse_args()

pubmlst_path = args.pubmlst
fastadir = args.fastadir




if __name__ == "__main__":
	main()


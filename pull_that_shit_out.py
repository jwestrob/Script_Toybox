from Bio import SeqIO
import os, sys, csv
import argparse

parser = argparse.ArgumentParser(description='Take a set of names (or a single name) and create a subset fasta.')

parser.add_argument('-i', metavar='--infile', nargs=1, help='File to make subset of.')
parser.add_argument('-o', metavar='--outfile', nargs=1, help='File to add subset to.')
parser.add_argument('-term', metavar='--term', nargs='+', type=str, help='Terms to match in rec ID',\
                    default=None)
parser.add_argument('-ids', metavar='--ids', nargs='+', type=str, default=None,\
                    help='Give me a list of IDs to include.')
parser.add_argument('-idc', nargs='?', type=str, default=None,\
                    help='A .csv file containing all your IDs. preferably in one row.')

args = parser.parse_args()
infile = args.i[0]
outfile = args.o[0]

if args.term is not None:
	term = args.term[0]
else:
	term = None
ids = args.ids
idc = args.idc

in_recs = list(SeqIO.parse(infile, 'fasta'))
try:
    out_recs = list(SeqIO.parse(outfile, 'fasta'))
except:
    out_recs = []
orig_len = len(out_recs)

if term is not None:
    for rec in in_recs:
        if term in rec.id:
            out_recs.append(rec)

if ids is not None:
    for id in ids:
        for rec in in_recs:
            if id == rec.id:
                out_recs.append(rec)

if idc is not None:
    idc_list = []
    print(idc)
    with open(idc, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for row in csv_reader:
            idc_list.append(row)

    idc_list = [item for sublist in idc_list for item in sublist]

    for id in idc_list:
        for rec in in_recs:
            if id == rec.id:
                out_recs.append(rec)


SeqIO.write(out_recs, outfile, 'fasta')
if len(out_recs)-orig_len == 0:
	print("Big whoops, no seqs!")
	print(ids)
print("Wrote " + str(len(out_recs)-orig_len) + " to " + outfile + ". Congratulation.")

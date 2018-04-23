import os, csv, sys, collections, argparse

parser=argparse.ArgumentParser(description='Parse and filter Anvio summary files.')


parser.add_argument('-infile', metavar='blast log infile,', help="Path to blastn log file")
parser.add_argument('-outfile', metavar='outfile name', help="Name of output file")
parser.add_argument('-exclude_ids', metavar='', nargs='+', type=str, help='Names of the genomes to exclude, as found in the Anvio summary infile')
parser.add_argument('-print_ids', metavar='Print the genome IDs and exit.', nargs='?', help='Dont let your dreams be dreams')

args=parser.parse_args()

infile = str(args.infile)
outfile = str(args.outfile)
exclude_ids = args.exclude_ids
if(args.print_ids != None):
    print_ids = str(args.print_ids)
else:
    print_ids = 0

if print_ids != 0:
    with open(infile, 'r') as readfile:
        mlst_read = csv.reader(readfile, delimiter="\t")
        pan_list = list(mlst_read)
    id_list = []
    for index, row in enumerate(pan_list):
        if index != 0:
            if row[3] not in id_list:
                id_list.append(row[3])
    print(id_list)
    sys.exit()


with open(infile, 'r') as readfile:
    mlst_read = csv.reader(readfile, delimiter="\t")
    pan_list = list(mlst_read)

header = pan_list[0]

PC_list = []
bad_list = []

for index, row in enumerate(pan_list):
    if index != 0:
        if row[1] not in PC_list:
            PC_list.append(row[1])
        if row[3] in exclude_ids and row[1] not in bad_list:
            bad_list.append(row[1])

good_list = []

for PC in PC_list:
    if PC not in bad_list:
        good_list.append(PC)

final_list = []
final_list.append(header)
for index, row in enumerate(pan_list):
    if index != 0:
        if row[1] in good_list:
            final_list.append(row)

print(final_list[1:3])

with open(outfile, 'w') as writefile:
    writer = csv.writer(writefile, delimiter='\t')
    for element in final_list:
        writer.writerow(element)

print('boogie')

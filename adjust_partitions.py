import pandas as pd
import os, sys

import argparse

parser = argparse.ArgumentParser(description='Take a partitions.nex file and a column mapping file from trimal (generated with -colnumbering flag) and adjusts it so you still have your partitions')
parser.add_argument('-p', metavar='[PARTITION FILE]', help='Nexus file of partitions for your dataset')
parser.add_argument('-c', metavar='[COLUMN MAP]', help='Column mapping file (should be all on one line, comma delimited, beginning with #ColumnMap)')
parser.add_argument('-o', metavar='[OUTPUT NAME]', default='trimmed_partitions.nex', help='Name for your outfile. Default: trimmed_partitions.nex')

args = parser.parse_args()

partitions = str(args.p)
cmap = str(args.c)
outfile = str(args.o)

cmap = pd.read_csv(cmap, header=None).iloc[0].tolist()
cmap[0] = int(cmap[0].split('\t')[-1])
with open(partitions, 'r') as infile:
    lines = infile.readlines()
lines

columns = list(map(lambda x: x.split('= ')[-1].split(';')[0].split('-'), lines[2:-1]))
def closest(list, Number):
    closest = list[0]
    for element in list:
        if element > Number:
            return list.index(closest)
        else:
            closest = element

    return


out_lines = lines[0:1]

for index, column in enumerate(columns):
    new_positions = []
    #print(column)
    if index == 0:
        new_positions.append(str(1))
    else:
        new_positions.append(str(1+closest_col))
    closest_col = closest(cmap, int(column[1]))
    new_positions.append(str(closest_col))
    out_line = '\tcharset part' + str(index) + ' = ' + '-'.join(new_positions) + ';\n'
    out_lines.append(out_line)

out_lines.append('end;')

with open(outfile, 'w') as outfile:
    for line in out_lines:
        outfile.write(line)

import os
import csv
import sys
import time

time = time.time()

with open('golden_ids.csv', 'rb') as infile:
    hits_table = csv.reader(infile, delimiter=",")
    id_list = list(hits_table)

for filename in os.listdir('/home/luisa/Whole_Genomes_2017/individual_genomes'):
    if filename.endswith('.fasta') and filename.split(".")[0] in id_list:
        print "whoop"
        os.system('cp ' + filename + ' /home/luisa/Whole_Genomes_2017/Golden_Set')

print("boogie")

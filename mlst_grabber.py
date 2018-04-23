import os
import sys
import csv
import time

with open('golden_ids.csv', 'rb') as idfile:
    ids = csv.reader(idfile, delimiter=",")
    id_list = list(ids)


with open('pubmlst.csv', 'rb') as mlst:
    mlst_read = csv.reader(mlst, delimiter=",")
    mlst_list = list(mlst_read)

header = mlst_list[0]

flatten = lambda id_list: [item for sublist in id_list for item in sublist]
id_list_flat = flatten(id_list)


new_mlst_list = []

for row in mlst_list:
    if row[0] in id_list_flat:
        for id_row in id_list:
             if row[0] == id_row[0]:
                 for element in id_row:
                     if element != id_row[0]:
                         row[0] += "_"
                         row[0] += element
        new_mlst_list.append(row)




with open('golden_info.csv', 'wb') as writefile:
    writer = csv.writer(writefile, delimiter=',')
    writer.writerow(header)
    for element in new_mlst_list:
        writer.writerow(element)

print 'boogie'

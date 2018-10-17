import sys
import os
import csv

def parser():
    clean_vars = []
    with open('filtered_ids.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\n')
        for row in reader:
            clean_vars.append(row)

    return clean_vars

def query(data):
    count = 0
    for isolate in data:
	if(count%100==0):
         print("Download " + str(count) + " Complete.")
        url = "http://rest.pubmlst.org/db/pubmlst_spneumoniae_isolates/isolates/" + isolate[0] + "/contigs_fasta?header=original_designation"
        filename = isolate[0] + '.fasta'
        command = "curl -i " + url + " > " + filename
	os.system(command)
	count += 1
    return

def main():
    data = parser()
    print("Parsed!")
    query(data)

if __name__ == "__main__":
    main()

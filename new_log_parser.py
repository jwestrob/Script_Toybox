import csv
import argparse
import time
import sys

parser=argparse.ArgumentParser(description='Parse and filter blast+/blastn output (-outfmt "7 qseqid sseqid evalue bitscore length pident qseq sseq").\n Please keep in mind that if you want to use this on other BLAST results you will currently need to manually enter the names of the genes youd like to filter.')


parser.add_argument('-logfile', metavar='blast log infile,', help="Path to blastn log file")
parser.add_argument('-ident_cutoff', metavar='Percent identity cutoff', help='float; percentage between 0 and 100')
parser.add_argument('-length_cutoff', metavar='Alignment length cutoff', nargs='?', default=None, help='int; recommended to determine length of query sequences prior to use')
parser.add_argument('-alleles', metavar='Names of the alleles you want to find', nargs='+', type=str, help='Names of the alleles you want to find, as they appear in the BLAST log')


args=parser.parse_args()

logfile = str(args.logfile)

ident_cutoff = float(args.ident_cutoff)

if(args.length_cutoff != None):
    length_cutoff = int(args.length_cutoff)
else:
    length_cutoff = 0

alleles = args.alleles
num_alleles = len(alleles)

def appender(row, allele_id, all_alleles):
    for allele in alleles:
        if allele_id == allele:
            index = alleles.index(allele_id)
            all_alleles[index].append(row)
            return
    print "Unknown allele ID found"
    sys.exit()
    return
    """
    else:
        print("Found an ID that doesn't exist in the record; investigate\n")
        print("Offending ID: ", allele_id)
        sys.exit()
    """

def writer(all_alleles):
    for index, allele_list in enumerate(all_alleles):
        allele_id = alleles[index]
        csvname = allele_id + '_filtered.csv'
        with open(csvname, 'w') as writefile:
            csvwriter = csv.writer(writefile, delimiter="\t")
            for element in allele_list:
                csvwriter.writerow(element)
    return

def main():
    all_alleles = [[] for i in range(0, num_alleles)]
    with open(logfile, 'r') as csvfile:
        log = csv.reader(csvfile, delimiter="\t")
        for row in log:
            if '#' in row[0]:
                header = True
            else:
                allele_id = row[0].split("_")[-1]
                header = False
            if header != True:
                if float(row[5]) >= ident_cutoff and int(row[4]) >= length_cutoff:
                    row[0] = row[0].split("_")[-1]
                    row[1] = row[1].split("_")[-1]
                    appender(row, allele_id, all_alleles)
    writer(all_alleles)
    print "Complete!"

if __name__ == "__main__":
    main()

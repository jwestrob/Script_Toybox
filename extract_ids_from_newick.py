import sys, csv
from Bio import Phylo

try:
    newick = sys.argv[1]
    outfile = sys.argv[2]
except:
    print("Give me a newick file to parse and an outfile name!")
    sys.exit()

tree = Phylo.read(newick, "newick")

out_ids = [leaf.name for leaf in tree.get_terminals()]

with open(outfile, 'w') as writefile:
    csvwriter = csv.writer(writefile, delimiter=",")
    csvwriter.writerow(out_ids)

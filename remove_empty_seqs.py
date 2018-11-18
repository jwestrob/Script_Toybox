from Bio import SeqIO
from collections import defaultdict
import os, sys

#print(sys.argv)
#sys.exit()

dedup_records = defaultdict(list)
for record in SeqIO.parse(sys.argv[1], "fasta"):
    # Use the sequence as the key and then have a list of id's as the value
    dedup_records[str(record.seq)].append(record.id)
with open(sys.argv[1], 'w') as output:
    for seq, ids in dedup_records.items():
        # Join the ids and write them out as the fasta
        output.write(">{}\n".format('|'.join(ids)))
        output.write(seq + "\n")


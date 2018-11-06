#!/bin/bash

for fastafile in 3.6k*;
do
	mkdir scanning/$fastafile.dir
	python /home/jwestrob/scripts/gene_scanner.py -hmm giy-yig.hmm -p $fastafile -fo scanning/$fastafile.dir/$fastafile.giy-yig.faa \
		-ht scanning/$fastafile.dir/$fastafile.hitstable.tsv -luis -threads 6
done

echo "boogie"

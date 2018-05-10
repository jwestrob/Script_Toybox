import argparse
import os, sys

parser = argparse.ArgumentParser(description="Takes a fasta file of raw reads, a fasta file with contigs, and an index name. (Runs bowtie2-build in the current directory, so make sure you're ok with lots of files popping up.) Aligns raws to contigs, separates out mapped/unmapped reads in .bam and .fa format.")

parser.add_argument('-raws', metavar='fasta raws', nargs=1, help="Path to raw reads to align.")
parser.add_argument('-contigs', metavar='fasta contigs', nargs=1, help="Path to contigs.")
parser.add_argument('-index_name', metavar='Bowtie index name (optional)', default="ctg_idx", help="Name for bowtie2 index.")
parser.add_argument('-no_build', action='store_true', default=False, help="If you already ran bowtie2-build and don't want to do it again.")

args = parser.parse_args()

raws = args.raws[0]
contigs = args.contigs[0]
idx = args.index_name

if args.no_build == False:
	#Build bowtie2 index
	os.system('bowtie2-build ' + contigs + ' ' + idx)

#Align with bowtie2
os.system('bowtie2 -x ' + idx + ' -p 64 -r ' + raws + ' --un ctg_raw_aln_unmapped.sam --al ctg_raw_aln_mapped.sam --quiet')

#Convert mapped sam to bam
os.system('samtools view -S -b ctg_raw_aln_mapped.sam > ctg_raw_aln_mapped.bam')

#Convert unmapped sam to bam
os.system('samtools view -S -b ctg_raw_aln_unmapped.sam > ctg_raw_aln_unmapped.bam')

#Convert mapped reads to fasta
os.system('samtools bam2fq ctg_raw_aln_mapped.bam | seqtk seq -A mapped.fa')

#Convert unmapped reads to fasta
os.system('samtools bam2fq ctg_raw_aln_unmapped.bam | seqtk seq -A unmapped.fa')

print('Boogie')


import argparse
import os, sys

parser = argparse.ArgumentParser(description="Takes a fasta file of raw reads, a fasta file with contigs, and an index name. (Runs bowtie2-build in the current directory, so make sure you're ok with lots of files popping up.) Aligns raws to contigs, separates out mapped/unmapped reads in .bam and .fa format.")

#parser.add_argument('-raws', metavar='fasta raws', nargs=1, help="Path to raw reads to align.")
parser.add_argument('-forward', metavar='Forward reads', nargs=1, help="Path to forward reads.")
parser.add_argument('-reverse', metavar='Reverse reads', nargs=1, help="Path to reverse reads.")
parser.add_argument('-contigs', metavar='fasta contigs', nargs=1, help="Path to contigs.")
parser.add_argument('-index_name', metavar='Bowtie index name (optional)', default="ctg_idx", help="Name for bowtie2 index.")
parser.add_argument('-no_build', action='store_true', default=False, help="If you already ran bowtie2-build and don't want to do it again.")
parser.add_argument('-threads', default=1, help="Number of threads to use.")


args = parser.parse_args()

#raws = args.raws[0]
forward = args.forward[0]
reverse = args.reverse[0]
contigs = args.contigs[0]
idx = args.index_name
threads = str(args.threads)


if args.no_build == False:
	#Build bowtie2 index
	print('bowtie2-build --threads ' + threads + ' '  + contigs + ' ' + idx)
	os.system('bowtie2-build --threads ' + threads + ' ' + contigs + ' ' + idx)

#Align with bowtie2
print('bowtie2 -x ' + idx + ' -p ' + threads + ' -1 ' + forward + ' -2 ' + reverse + ' | shrinksam | sambam > contig_map.shrink.sort.bam')
os.system('bowtie2 -x ' + idx + ' -p ' + threads + ' -1 ' + forward + ' -2 ' + reverse + ' | shrinksam | sambam > contig_map.shrink.sort.bam')

#Separate unmapped read alignment
print('samtools view -@ ' + threads + ' -f 4 contig_map.shrink.sort.bam > unmapped.bam')
os.system('samtools view -@ ' + threads + ' -f 4 contig_map.shrink.sort.bam > unmapped.bam')

#Separate mapped read alignment
print('samtools view -@ ' + threads + ' -bSq 2 -F 4 contig_map.shrink.sort.bam > mapped.bam')
os.system('samtools view -@ ' + threads + ' -bSq 2 -F 4 contig_map.shrink.sort.bam > mapped.bam')




#Convert mapped reads to fastq
print('samtools fastq -@ ' + threads + ' mapped.bam -1 MAPPED.1.fastq.gz -2 MAPPED.2.fastq.gz -0 /dev/null -s /dev/null -n')
os.system('samtools fastq -@ ' + threads + ' mapped.bam -1 MAPPED.1.fastq.gz -2 MAPPED.2.fastq.gz -0 /dev/null -s /dev/null -n')

#Convert unmapped reads to fastq
print('samtools fastq -@ ' + threads + ' unmapped.bam -1 UNMAPPED.1.fastq.gz -2 UNMAPPED.2.fastq.gz -0 /dev/null -s /dev/null -n')
os.system('samtools fastq -@ ' + threads + ' unmapped.bam -1 UNMAPPED.1.fastq.gz -2 UNMAPPED.2.fastq.gz -0 /dev/null -s /dev/null -n')

print('Boogie')


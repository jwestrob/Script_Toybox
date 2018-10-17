import os, sys

for dirname in os.listdir('.'):
	if len(dirname.split('.')) == 4:
		sample = dirname.split('-')[0]
		os.chdir('/mnt/scratch/WGS_Harrison/dereplicated_genomes/' + dirname)
		fastaname = dirname.replace('-dir', '')
		os.system('bowtie2-build ' + fastaname + ' ' + sample)
		print('bowtie2 -x ' + sample + ' -p 64 -r ../../clean_raws/' + sample + '*  --no-unal -S ' + sample + '-aln.sam')
		os.system('bowtie2 -x ' + sample + ' -p 64 -r ../../clean_raws/' + sample + '*  --no-unal -S ' + sample + '-aln.sam')
		os.system('samtools view -S -b ' + sample + '-aln.sam -o ' + sample + '.sorted.mapped.bam')
		os.system('samtools bam2fq ' + sample + '.sorted.mapped.bam | seqtk seq -A > recruited_reads.fa')
		os.system('/home/jaw293/idba/bin/idba_ud --mink 40 --maxk 100 --min_count 1 --step 15 --min_contig 500 -r recruited_reads.fa --num_threads 64')
		os.system('cp out/contig.fa ' + sample + '-reassembly.fa')
		os.chdir('/mnt/scratch/WGS_Harrison/dereplicated_genomes')
	else:
		continue


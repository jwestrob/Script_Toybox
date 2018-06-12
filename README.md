# Script_Toybox
A bunch of stuff I've made over time working in the Hiller lab. May or may not be useful to me or anyone in the future; some of these were relatively early creations of mine.

## Align_and_separate.py
For when you have raw reads (fasta format) and a contigs file (fasta format) and you'd like to separate out the reads mapping to each contig.
Very useful for metagenomic re-assembly (post-binning,  pre-scaffolding).

## blastlog_to_MFA.py
Takes a blastlog from command-line blastn (using -outfmt 7 with query/hit nucleotide sequences returned) and gives you a multifasta ready to be aligned.
TODO: Include blastn command line options for input file. If you want to run this, see "parser_and_etc.py".

## contig_finder.py
If you run this in a directory with FASTA files, it'll print a list of filenames followed by the number of contigs in that file.

## correcter.py
Don't worry about this one. It was to fix some genomes I downloaded from PubMLST.

## file_sorter.py
Don't worry about this one. It was for operations on a database of pneumo genomes.

## GENBANK_PARSER.py
Made to operate on blastn results from pneumococcal insertion sequences. Input file format is... ill-defined. 
TODO: Include specification of input file format.

## get_concatenated_proteins.sh
Part of my anvio-based phylgenomics workflow. Anvi'o is required. 
Pulls out a concatenated alignment of large ribosomal proteins L1-L6 in nucleotide and amino acid FASTA format.

## get_extended_proteins.sh
Part of my anvio-based phylogenomics workflow. Anvi'o is required.
Same as get_concatenated_proteins.sh, but pulls out every ribosomal protein instead of LSU L1-L6.

## mlst_grabber.py
Don't worry about this one. It's for pulling out a subset of metadata from a larger .csv file. Not necessary if you know pandas lmao

## murm_functions.py
Don't worry about this one. I made it thinking we'd need it for analysis but we didn't need it. I'm keeping it around just in case.

## NCBI_fetcher.py
Given a comma-deliimited list of NCBI accession IDs, downloads the requested FASTAS using Biopython's Entrez module.

## parser_and_etc.py
BLASTn log parser. Run "python parser_and_etc.py -h" to see command-line parameters for BLASTn to get the proper format for your log.
Allows for lots of cutoffs and specifications so you can find whatever kind of protein you want. Enjoy.

## plot_checkm.R
Plots contamination/completeness measures generated by checkM. (dRep outputs these as part of its pipeline; it's worth including if you run the program yourself also.)
If by chance you're running checkM and want to use this but are unsure of the parameters you need to generate the proper format table, post an issue.

## pubmlst_parser.py
This is not for you unless you're in the Hiller lab, in which case just email me.

## query.py
This pulls something off pubmlst but I didn't write down exactly what or why. I'll get back to this.
TODO: Find out what this is

## read_length_histogram.py
Someone else's code I pulled off stack overflow. I DID NOT WRITE THIS SCRIPT! It will give you a histogram of contig lengths, but Anvi'o will do that too, and much better.

## run_anvio_dbandstats.py
Part of my anvio-based phylogenomics workflow. Anvi'o is required.
Runs several initial analyses for Anvi'o contigs databases, including ORF calling, tetramer frequency calculation, HMM scans, and NCBI COG annotations.

## run_bowtie2.py
Runs bowtie2-build and bowtie2 to make alignments. This version is highly incomplete; specific paths are included for only one analysis.

## Run_DADA2.R
Runs the DADA2 pipeline on a given set of FASTQ files (QC not required; part of the pipeline). 
Outputs a .pdf of plots, an OTU matrix, a taxonomy matrix, an .rds file containing a phyloseq object, a .fasta of 16S gene sequences, and an alignment of those seqeuences ready for phylogenetic analysis.

## run_fastx_qual.py
A very antagonistic, unhelpful program for running fastq_quality_filter without writing out the entire command. I recommend you ignore this script and do it yourself.

## scaffolds2bin.py
Takes scaffolds2bin.txt files from DASTool and converts them to FASTA files.

## weirdo_finder.py
Takes anvi'o pangenomics summary files and finds proteins unique to a particular strain. I used it in my Streptococcus thermophilus research.
If you want to use this for something, post an issue and I'll help out.
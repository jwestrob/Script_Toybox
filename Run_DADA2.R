if (!require("pacman")) install.packages("pacman")
pacman::p_load(argparse, dplyr, hash, ggplot2, dada2, phyloseq, stringr, seqinr)
require(Biostrings)

#load argument parser
parser <- ArgumentParser(description="Run the DADA2 pipeline on a given set of FASTQ sequences (unfiltered). Gives you a PDF with composition barplots and alpha diversity (-p) and a summary table containing taxon abundances (-s) as well as some other optional things if you feelin frisky. Please ensure that your fastq files end with _R1_001.fastq/_R2_001.fastq; otherwise go in and change the 'pattern' argument on lines 43/44.")

#Specify options
parser$add_argument("-s", "--sequences", type="character", nargs=1,
    help="Path to fastq files (Ideally the full path to the folder containing all your files).")
parser$add_argument("-st", "--summarytable", type="character", nargs=1,
    help="Name for your summary table (.tsv will be appended).")
parser$add_argument("-p", "--pdfout", type="character", nargs=1,
    help="PDF outfile name for plots.")
parser$add_argument("-ps", "--phyloseq", type="character", nargs=1,
    help="Name of file to save phyloseq object to. include .rds extenison.")
parser$add_argument("-f", "--fastaprefix", type="character", nargs=1,
    help="Name for fasta output files (used as prefix for aligned and unaligned FASTAs). REQUIRED")
parser$add_argument("-RT", "--TRAIN", type="character", nargs=1,
    help="Path to SILVA training data. REQUIRED")
parser$add_argument("-RS", "--SPECIES", type="character", nargs=1,
    help="Path to SILVA species-level taxonomy data. REQUIRED")
parser$add_argument("-t", "--taxa", type="character", default=NULL,
    help="Optional: Print original taxa table (ASV sequence format). Include .tsv extension.")
parser$add_argument("-o", "--otu", type="character", default=NULL,
    help="OPTIONAL: Print original OTU/ASV table (ASV sequence format). Include .tsv extension.")


#load arguments as variables (in array "args")
args <- parser$parse_args()

fastqpath <- args$sequences
fastaprefix <- args$fastaprefix
summarytable <- args$summarytable
pdfout <- args$pdfout
silva_train <- args$TRAIN
silva_species <- args$SPECIES
phylo_file <- args$phyloseq

if(length(args$taxa) == 0){
  tax = FALSE
}else{
  tax = TRUE
}
if(length(args$otu) == 0){
  otu = FALSE
}else{
  otu = TRUE
}

if(tax){
  tax_table = args$taxa
}else{
  tax_table = NULL
}
if(otu){
  otu_table = args$otu
}else{
  otu_table = NULL
}

########################################
#              READ DATA               #
########################################

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(fastqpath, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastqpath, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
pdf(pdfout)
print(plotQualityProfile(fnFs) + ggtitle('Forward read quality profiles.'))
print(plotQualityProfile(fnRs) + ggtitle('Reverse read quality profiles.'))

filt_path <- file.path(fastqpath, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

#Perform error rate estimation

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

print(plotErrors(errF, nominalQ=TRUE) + ggtitle('Error rates - Forward reads'))
print(plotErrors(errR, nominalQ=TRUE) + ggtitle('Error rates - Reverse reads'))

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

print("Sequence length distribution: ", table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithread=TRUE)
taxa <- addSpecies(taxa, silva_species)

samples.out <- rownames(seqtab.nochim)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(taxa))

print(plot_richness(ps, measures=c("Shannon", "Simpson")) + theme_bw() + ggtitle("Alpha diversity; Shannon and Simpson"))

#ord.PCA.bray <- ordinate(ps, method="NMDS", distance="manhattan")

#print(plot_ordination(ps, ord.PCA.bray, title="Bray PCA") + theme_bw() + ggtitle("NMDS - Manhattan Distance"))

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
top80 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:80]

ps.alltaxa <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

#Plot non-normalized top 20 Genera
ps.top20 <- prune_taxa(top20, ps.alltaxa)
print(plot_bar(ps.top20, fill="Genus") + theme_bw() + ggtitle("Non-normalized composition of top 20 genera") + theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90, hjust = 1)))

#Plot non-normalized top 20 Families
print(plot_bar(ps.top20, fill="Family") + theme_bw() + ggtitle("Non-normalized composition of top 20 families") + theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90, hjust = 1)))

#Normalize top 20 and plot
ps.top20 <- transform_sample_counts(ps.top20, function(OTU) OTU/sum(OTU))
plot_bar(ps.top20, fill="Genus") + ggtitle("Normalized composition of top 20 genera") +
    guides(fill=guide_legend(nrow=40,byrow=TRUE)) + theme_bw() +
    theme(legend.key.size = unit(2, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))
plot_bar(ps.top20, fill="Family") + ggtitle("Normalized composition of top 20 genera") +
    guides(fill=guide_legend(nrow=40,byrow=TRUE)) + theme_bw() +
    theme(legend.key.size = unit(2, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))


#Plot non-normalized top 80 Genera
ps.top80 <- prune_taxa(top80, ps.alltaxa)
print(plot_bar(ps.top80, fill="Genus") + theme_bw() + ggtitle("Non-normalized composition of top 80 genera") +
    theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=guide_legend(nrow=40, byrow=TRUE)))
print(plot_bar(ps.top80, fill="Family") + theme_bw() + ggtitle("Non-normalized composition of top 80 families") +
    theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=guide_legend(nrow=40, byrow=TRUE)))

#Normalize top 80 and plot
ps.top80 <- transform_sample_counts(ps.top80, function(OTU) OTU/sum(OTU))

plot_bar(ps.top80, fill="Genus") + ggtitle("Normalized composition of top 80 genera") +
    guides(fill=guide_legend(nrow=40,byrow=TRUE)) + theme_bw() +
    theme(legend.key.size = unit(2, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))
plot_bar(ps.top80, fill="Family") + ggtitle("Normalized composition of top 80 genera") +
    guides(fill=guide_legend(nrow=40,byrow=TRUE)) + theme_bw() +
    theme(legend.key.size = unit(2, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))

Species <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.Species <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

print(plot_bar(ps.Species, fill="Species") + theme_bw() + ggtitle("All species (normalized)") +
    theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=guide_legend(nrow=40, byrow=TRUE)))

if(tax){
  #Convert tax_table to matrix
  tax_mat = as(tax_table(ps), 'matrix')

  #Write that nonsense down
  write(tax_mat, file=tax_table, sep='\t')
}
if(otu){
  #Same with OTU table
  otu_mat = as(otu_table(ps), 'matrix')

  #Write that nonsense down again
  write(tax_mat, file=otu_table, sep='\t')
}

saveRDS(ps, phylo_file)

h <- hash()
for(rowname in row.names(tax_mat)){
    taxonomy_nonsense = sapply(colnames(tax_mat), function(x) paste(tax_mat[rowname, x], sep=';'))
    tax_info = paste(taxonomy_nonsense, collapse=';')
    h[rowname] = tax_info
}

new_header <- lapply(colnames(otu_mat), function(x) h[[x]])
otu_mat_newheader = otu_mat
colnames(otu_mat_newheader) = new_header

Total = c()
for(i in 1:nrow(otu_mat_newheader)){
    sum = sum(otu_mat_newheader[i,])
    Total[[i]] = sum
}

otu_mat_wsums = cbind(Total, otu_mat_newheader)
#print(head(otu_mat_wsums))
otu_mat_wsums = cbind(Sample_ID = rownames(otu_mat_wsums), otu_mat_wsums)

for(i in 1:length(colnames(otu_mat_wsums))){
    colnames(otu_mat_wsums)[i] = paste(strsplit(colnames(otu_mat_wsums)[[i]],'.',fixed=TRUE)[[1]],collapse=';')
}

rownames(otu_mat_wsums) <- NULL

to_save <- data.frame(otu_mat_wsums)
write.table(to_save, file=summarytable, sep='\t', quote=FALSE)
dev.off()

aln_df = data.frame()
for(i in keys(h)){
    aln_df[i,"16S"] = i
    aln_df[i,"Tax"] = h[[i]]
}
rownames(aln_df) <- NULL
head(aln_df)

seqs <- list()
tax <- list()
#Uniquify taxonomy IDs for later processing with iq-tree
tax_char <- lapply(aln_df["Tax"], function(x) as.character(x))

wombo <- make.unique(as.character(unlist(tax_char)))
print(wombo)
aln_df["New_Tax"] <- wombo

for(i in 1:nrow(aln_df["16S"])){
    seqs[i] <- aln_df[i, "16S"]
    tax[i] <- aln_df[i, "New_Tax"]
}

#Remove list comprehension
seqs <- unlist(seqs)
#Turn into biostrings-friendly format
seqs <- sapply(seqs, function(x) x <- DNAString(x))

#Save seqeuences as fasta
unaligned <- paste0(fastaprefix, '_unaligned.fa', sep='')
write.fasta(seqs, tax, unaligned, open='w')
#Align with MUSCLE
aligned <- paste0(fastaprefix, '_ALIGNED.fa', sep='')
command <- paste('muscle -in ', unaligned, ' -out ', aligned, sep='')
system(command)

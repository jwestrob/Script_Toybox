if (!require("pacman")) install.packages("pacman")
pacman::p_load(argparse, tidyverse, ggplot2, reshape2)


parser <- ArgumentParser(description="Takes results from checkm (results.tsv) and gives you contamination/completeness as a bar plot.")

parser$add_argument("-r", "--results", type="character", nargs=1,
    help="Path to results.tsv file.")

parser$add_argument("-p", "--pdf", type="character", nargs=1,
    help="Name of output pdf (include .pdf extension please)")

parser$add_argument("-t", "--title", type="character", nargs=1, 
    help="Title for your plot.")

args <- parser$parse_args()

checkm <- args$results
pdfout <- args$pdf
title = args$title

results <- read.csv(checkm, sep='\t')
head(results)

dfm <- melt(results[,c('Bin.Id','Completeness','Contamination')],id.vars = 1)

p <- ggplot(dfm,aes(x = Bin.Id,y = value)) +
    geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)

pdf(pdfout)
print(p)

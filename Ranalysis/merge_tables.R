# Use this script to make the gene and transcript abundance table for EQA
library(ggplot2)
require(reshape2)
require(scales)
require(plyr)
library(stringr)
#library(tidyverse)

# swt working directory
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/stringtie/merge")
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/stringtie/merge")

gene_gtf <- read.csv("EQA_gene.txt", sep = " ", header = FALSE)
gene_count <- read.csv("gene_count_matrix.csv", sep = ",", header = TRUE, check.names = FALSE)
colnames(gene_count) <- c("gene_id","202201",  "202210", "202211", "202212", "202213", "202214", "202215", "202216", "202217", "202218", 
                          "202219", "202202",  "202220", "202221", "202222", "202223", "202224", "202203",  "202204",  "202205",  "202206",
                          "202207",  "202208",  "202209")
transcript_gtf <- read.csv("EQA_transcript.txt", sep = " ", header = FALSE)
transcript_count <- read.csv("transcript_count_matrix.csv", sep = ",", header = TRUE, check.names = FALSE)
colnames(transcript_count) <- c("transcript_id","202201",  "202210", "202211", "202212", "202213", "202214", "202215", "202216", "202217", "202218", 
                                "202219", "202202",  "202220", "202221", "202222", "202223", "202224", "202203",  "202204",  "202205",  "202206",
                                "202207",  "202208",  "202209")


## make the gene.counts.abundance file
# split the gene_id to only save the ensemble id
gene_count$gene <- str_split_fixed(gene_count$gene_id, "[|]", 2)[,1]
#remove the original gene_id column
gene_count <- gene_count[,-which(colnames(gene_count) %in% c("gene_id"))]
#reshape the data frame 
gene_count.m <- reshape2::melt(gene_count)
#rename the two data frames' headers for merging
colnames(gene_gtf) <- c("sample", "gene_id", "FPKM", "TPM")
colnames(gene_count.m) <- c("gene_id", "sample", "counts")
#merge and reorder the columns
gene_merge <- merge(gene_gtf, gene_count.m, by.x = c("sample", "gene_id"), by.y = c("sample", "gene_id"), all=FALSE, all.y=TRUE)
#gene_merge <- merge(gene_gtf, gene_count.m, all=FALSE)
gene_merge <- gene_merge[,c("sample", "gene_id", "counts", "FPKM", "TPM")]
#export the final text
write.table(gene_merge, "gene.counts.abundance.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## make the transcript.counts.abundance file
transcript_count.m <- reshape2::melt(transcript_count)
colnames(transcript_gtf) <- c("sample", "transcript_id", "FPKM", "TPM")
colnames(transcript_count.m) <- c("transcript_id", "sample", "counts")
transcript_merge <- merge(transcript_gtf, transcript_count.m, by = c("sample", "transcript_id"))
transcript_merge <- transcript_merge[,c("sample", "transcript_id", "counts", "FPKM", "TPM")]
write.table(transcript_merge, "transcript.counts.abundance.txt", sep = "\t", row.names = FALSE)

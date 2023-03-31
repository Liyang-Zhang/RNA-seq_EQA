# Use this script to do PCA analysis
# Using DeSeq2 package
library(DESeq2)
library(ggplot2)
library(ggrepel)

# set the working directory
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/counts/")
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/counts")

#Use Mov10 as an example
count=read.table("/dssg/home/acct-medkwf/medkwf4/ngs/hbctraining/rnaseq/results/counts/Mov10.clean.matrix", header = TRUE, quote = '\t')
head(count)
countMatrix=as.matrix(count[2:13])
rownames(countMatrix)=count[,1]
sampleNames=c(colnames(count)[2:13])
table2=data.frame(name=c("SRR960459", "SRR960460", "SRR960461", "SRR960462", "SRR960463", "SRR960464", "SRR960465", "SRR960466", "SRR960467", "SRR960468", "SRR960469", "SRR960470"),condition=c(rep("high", 6), rep("normal", 6)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)
pdf("/dssg/home/acct-medkwf/medkwf4/ngs/hbctraining/rnaseq/results/counts/allGene.PCA.pdf")
plotPCA(rld, intgroup=c("name","condition"))
dev.off()

#Use STAR and featureCounts matrix
count=read.table("all_featurecounts_clean.Rmatrix", header = TRUE, quote = '\t', check.names = FALSE)
head(count)
countMatrix=as.matrix(count[2:25])
rownames(countMatrix)=count[,1]
sampleNames=c(colnames(count)[2:25])
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

nudge <- position_nudge(y = 1)
options(ggrepel.max.overlaps = Inf)
plotPCA(rld, intgroup=c("name","condition")) +
  geom_text_repel(aes(label=name))

# pdf("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/counts/RNAseq_EQA_PCA.pdf")
# plotPCA(rld, intgroup=c("name","condition")) +
#   geom_text_repel(aes(label=name))
# dev.off()


##############################################################################################
# PCA analysis for stringtie transcriptome and gene level
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/stringtie/merge")
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/stringtie/merge")
count_stringtie_trans <- read.csv("transcript_count_matrix.csv", header = TRUE, check.names = FALSE)
count_stringtie_genes <- read.csv("gene_count_matrix.csv", header = TRUE, check.names = FALSE)
head(count_stringtie_genes)
countMatrix=as.matrix(count_stringtie_genes[2:25])
rownames(countMatrix)=count_stringtie_genes[,1]
sampleNames=c(colnames(count_stringtie_genes)[2:25])
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

nudge <- position_nudge(y = 1)
options(ggrepel.max.overlaps = Inf)
p <- plotPCA(rld, intgroup=c("name","condition")) +
    geom_text_repel(aes(label=name))

ggsave(p, filename = "PCA_geneLevel_stringtie_cn500.pdf",width = 11,height =8)

##############################################################################################
#PCA analysis for transcriptome level results (salmon)
library(tximport)
library(GenomicFeatures)
library(data.table)

setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/salmon")
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/salmon")

##载入transcript_id和symbol的对应关系文件
t2s <- fread("/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/t2g.all.txt", data.table = F, header = F); head(t2s)

##找到所有quant.sf文件所在路径  导入salmon文件处理汇总
files <- list.files(pattern="*quant.sf",recursive=T, full.names = T)
txi <- tximport(files, type = "salmon", tx2gene = t2s)

##提取文件夹中的样品名作为counts行名
cn <- sapply(strsplit(files,'\\/'), function(x) x[length(x)-1])
colnames(txi$counts) <- gsub('.salmon','',cn); colnames(txi$counts)

##提取出counts/tpm表达矩阵
counts <- as.data.frame(apply(txi$counts,2,as.integer)) #将counts数取整
rownames(counts) <- rownames(txi$counts) 
tpm <- as.data.frame(txi$abundance)  ###abundance为基因的Tpm值
colnames(tpm) <- colnames(txi$counts)

##进行PCA分析
countMatrix=as.matrix(counts)
sampleNames=c(colnames(counts))
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

options(ggrepel.max.overlaps = Inf)
plotPCA(rld, intgroup=c("name","condition")) +
  geom_text_repel(aes(label=name))

##进行PCC分析
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results")
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results")
dir.create("Ranalysis")
setwd("Ranalysis")
dat <- log2(tpm+1)
#取高表达量前500基因
dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),]#取高表达量前500基因
M <- cor(dat_500)
p2 <-pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        fontsize=7)
ggsave(p2,filename = 'check_cor_top500.pdf',width = 7.5,height =6)
#取所有基因
dat_all <- dat[names(sort(apply(dat,1,mad),decreasing = T)),]
M <- cor(dat_all)
p3 <-pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        fontsize=7)
ggsave(p3,filename = 'check_cor_all.pdf',width = 7.5,height =6)

####################################################################
#只对spikein做分析
counts_spikein <- counts[1:92,]
countMatrix_spikein=as.matrix(counts_spikein)
sampleNames_spikein=c(colnames(counts_spikein))
table2=data.frame(name=sampleNames_spikein,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames_spikein
dds <- DESeqDataSetFromMatrix(countMatrix_spikein, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

options(ggrepel.max.overlaps = Inf)
plotPCA(rld, intgroup=c("name","condition")) +
  geom_text_repel(aes(label=name))

tpm_spikein <- tpm[1:92,]
dat_spikein <- log2(tpm_spikein+1)
dat_spikein_all <- dat[names(sort(apply(dat_spikein,1,mad),decreasing = T)),]
M <- cor(dat_spikein_all)
p3 <-pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        fontsize=7)
ggsave(p3,filename = 'check_cor_spikein.pdf',width = 7.5,height =6)

####################################################################
#对瑞普的结果做PCA分析
library(tidyr)
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/anotherlab/results")
count_repu_genes=read.csv2("genes.counts.abundance.txt", header = TRUE, sep = '\t', check.names = FALSE)
count_repu_genes <- count_repu_genes[,c(1,2,3)]
count_repu_genes <- pivot_wider(count_repu_genes, names_from = "#sample", values_from = "counts")

count_repu_genes.matrix=as.matrix(count_repu_genes[2:25])
sampleNames=c(colnames(count_repu_genes)[2:25])
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(count_repu_genes.matrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

nudge <- position_nudge(y = 1)
options(ggrepel.max.overlaps = Inf)
p <- plotPCA(rld, intgroup=c("name","condition")) +
    geom_text_repel(aes(label=name))
ggsave(p,filename = 'repu_gene_PCA.pdf',width = 11,height =8)


count_repu_transcript <- read.csv2("transcripts.counts.abundance.txt", header = TRUE, sep = '\t', check.names = FALSE)
count_repu_transcript <- count_repu_transcript[,c(1,2,3)]
count_repu_transcript <- pivot_wider(count_repu_transcript, names_from = "#sample", values_from = "counts")
count_repu_transcript.matrix=as.matrix(count_repu_transcript[2:25])
sampleNames=c(colnames(count_repu_transcript)[2:25])
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(count_repu_transcript.matrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

nudge <- position_nudge(y = 1)
options(ggrepel.max.overlaps = Inf)
p2 <- plotPCA(rld, intgroup=c("name","condition")) +
  geom_text_repel(aes(label=name))
ggsave(p2,filename = 'repu_transcript_PCA.pdf',width = 11,height =8)

#对自己的回传结果做PCA分析
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/salmon")
count_transcript <- read.csv2("transcript.counts.abundance.txt", header = TRUE, sep = '\t', check.names = FALSE)
count_transcript <- count_transcript[,c(1,2,3)]
count_transcript$counts <- as.integer(count_transcript$counts)
count_transcript <- pivot_wider(count_transcript, names_from = "#sample", values_from = "counts")
count_transcript.matrix=as.matrix(count_transcript[2:25])
sampleNames=c(colnames(count_transcript)[2:25])
table2=data.frame(name=sampleNames,condition=as.character(seq(1, 24, by=1)))
rownames(table2)=sampleNames
dds <- DESeqDataSetFromMatrix(count_transcript.matrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds)

nudge <- position_nudge(y = 1)
options(ggrepel.max.overlaps = Inf)
p2 <- plotPCA(rld, intgroup=c("name","condition")) +
  geom_text_repel(aes(label=name))
ggsave(p2,filename = 'salmon_transcript_PCA.pdf',width = 11,height =8)

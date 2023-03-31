###########################################################################################
##使用DEseq2对STAR比对，stringtie计数的结果做差异表达分析
library(DESeq2)
library(ggplot2)
library(ggrepel)
library("BiocParallel")
library(dplyr)


##去除所有环境变量
#rm(list = ls()) 
#options(stringsAsFactors = F)

setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/stringtie/merge")
#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/stringtie/merge")
##读取gene count matrix 和 tpm count matrix
gene_count <- read.csv("gene_count_matrix.csv", sep = ",", header = TRUE, check.names = FALSE)
tpm_count <- read.csv("gene_tpm_matrix.csv", sep = ",", header = TRUE, check.names = FALSE)

##提取出counts/tpm表达矩阵
rownames(gene_count) <- gene_count$gene_id
counts <- gene_count[,-1]
rownames(tpm_count) <- tpm_count$gene_id
tpm <- tpm_count[,-1]

##针对reads和TPM对原始matrix进行过滤，一共61644
keep_tpm <- rowSums(tpm) >= 0.5*ncol(tpm) ##TPM平均大于0.5
keep_counts <- rowSums(counts) >= 1*ncol(counts) ##reads平均大于1.5

keep <- keep_counts[order(factor(names(keep_counts)))] & keep_tpm[order(factor(names(keep_tpm)))] ##26370
counts <- counts[names(keep)[keep],]
#keep <- keep_tpm & keep_counts ##问题出在这里没有sort，两个逻辑向量的顺序不一致，取并集是出错
#counts_test <- counts[keep_tpm,]  ## 这种保留方式会出bug，不把True保留，还不知道问题出在哪

##制作matrix，进行分组，注意stringtie的列名
countMatrix=as.matrix(counts)
sampleNames=c(colnames(counts))
table2=data.frame(name=sampleNames,
                  condition=c("T1","T3","T3","T3","T4","T4","T4","T5","T5","T5","T6","T1","T6","T6","C2","C2","C2","T1",
                              "C1","C1","C1","T2","T2","T2"))
rownames(table2)=sampleNames

##构建dds DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = table2,
                              design = ~ condition)

dds <- DESeq(dds,quiet = F)

##分组输出结果
res <- results(dds, contrast = c("condition", "T6", "C2"))  #选择需要进行DEG分析的组
resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)

##筛选条件设置
log2FC_cutoff = log2(2)
padj_cutoff = 0.001
##选取差异分析结果，修改输出的表格
need_DEG <- DEG_DEseq2[,c(2,5,6)] #选取log2FoldChange, p, padj信息
colnames(need_DEG) <- c('log2FC','p_value', 'Q_value') 
need_DEG$gene_id <- row.names(need_DEG)
need_DEG$sample_pair <- 'T6/C2'
need_DEG$type  <- as.factor(ifelse(need_DEG$Q_value < padj_cutoff & abs(need_DEG$log2FC) > log2FC_cutoff,
                                   ifelse(need_DEG$log2FC > log2FC_cutoff ,'Upregulate','Downregulate'),'NonDE'))

#拆分原始gene_id 到 HGNC symbol 和 ensembl id
need_DEG$gene_symbol <- str_split_fixed(need_DEG$gene_id, "[|]", 2)[,2]
need_DEG$gene_id <- str_split_fixed(need_DEG$gene_id, "[|]", 2)[,1]

#改变列序
#need_DEG_T1C1 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T2C2 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T3C2 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T4C2 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T5C2 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]
need_DEG_T6C2 <- need_DEG[,c("sample_pair", "gene_id", "gene_symbol", "log2FC", "p_value", "Q_value", "type")]

# 整合总表
need_DEG_merge_all <- do.call("rbind", list(need_DEG_T1C1, need_DEG_T2C2, need_DEG_T3C2, need_DEG_T4C2, need_DEG_T5C2, need_DEG_T6C2))

write.table(need_DEG_merge_all, "Differential.expressed.genes.9.23.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#结果做火山图
title <- paste0('T6/C2',
                '\n Up :  ',nrow(need_DEG_T6C2[need_DEG_T6C2$type =='Upregulate',]) ,
                '\n Down : ',nrow(need_DEG_T6C2[need_DEG_T6C2$type =='Downregulate',]),
                '\n FoldChange >= ',round(2^log2FC_cutoff,3))

g <- ggplot(data=need_DEG_T6C2, 
            aes(x=log2FC, y=-log10(Q_value), 
                color=type)) +
  #点和背景
  geom_point(alpha=0.4, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( P.adjust )") +
  #标题文本
  ggtitle( title ) +
  #分区颜色                  
  scale_colour_manual(values = c('blue','grey','red'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(), #不显示图例标题
        legend.position="right")  #图例位置
g

ggsave(g,filename = 'volcano_padj_T6C2_new.pdf',width =8,height =7.5)

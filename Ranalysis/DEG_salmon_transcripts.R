###########################################################################################
##使用DEseq2对Salmon结果做差异表达分析
library(DESeq2)
library(ggplot2)
library(ggrepel)
library("BiocParallel")
library(dplyr)


#去除所有环境变量
rm(list = ls()) 
options(stringsAsFactors = F)

#setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/CG0701-R22024472-202208251343-1/results/salmon")
setwd("/dssg/home/acct-medkwf/medkwf4/data/cgdata/bcl/cn500/results/salmon")
#不进行转录本水平到基因水平的转换
t2t <- fread("/dssg/home/acct-medkwf/medkwf4/reference/rna_EQA/t2t.all.txt", data.table = F, header = F); head(t2t)

##找到所有quant.sf文件所在路径  导入salmon文件处理汇总
files <- list.files(pattern="*quant.sf",recursive=T, full.names = T)
txi <- tximport(files, type = "salmon", tx2gene = t2t)
##提取文件夹中的样品名作为counts列名
cn <- sapply(strsplit(files,'\\/'), function(x) x[length(x)-1])
colnames(txi$counts) <- gsub('.salmon','',cn); colnames(txi$counts)

##提取出counts/tpm表达矩阵
counts <- as.data.frame(apply(txi$counts,2,as.integer)) #将counts数取整
rownames(counts) <- rownames(txi$counts) 
tpm <- as.data.frame(txi$abundance)  ###abundance为基因的Tpm值
colnames(tpm) <- colnames(txi$counts)

##针对reads和TPM对原始matrix进行过滤
keep_tpm <- rowSums(tpm) >= 0.5*ncol(tpm) ##TPM平均大于0.5, 留下87616
keep_counts <- rowSums(counts) >= 1*ncol(counts) ##reads平均大于1.5, 留下108950
# keep <- keep_tpm & keep_counts ##留下86258
# counts <- counts[keep,]
keep <- keep_counts[order(factor(names(keep_counts)))] & keep_tpm[order(factor(names(keep_tpm)))] ##86258
counts <- counts[names(keep)[keep],]

##制作matrix，进行分组
countMatrix=as.matrix(counts)
sampleNames=c(colnames(counts))
table2=data.frame(name=sampleNames,
                  condition=c(rep("T1",3), rep("C1",3), rep("T2",3), rep("T3",3), rep("T4",3), rep("T5",3),rep("T6",3),rep("C2",3)))
rownames(table2)=sampleNames

##构建dds DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = table2,
                              design = ~ condition)

dds <- DESeq(dds,quiet = F)

res <- results(dds, contrast = c("condition", "T1", "C1"))  #选择需要进行DEG分析的组
resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)

##筛选条件设置
log2FC_cutoff = log2(2)
padj_cutoff = 0.001
##选取差异分析结果，修改输出的表格
need_DEG <- DEG_DEseq2[,c(2,5,6)] #选取log2FoldChange, p, padj信息
colnames(need_DEG) <- c('log2FC','p_value', 'Q_value') 
need_DEG$transcript_id <- row.names(need_DEG)
need_DEG$sample_pair <- 'T1/C1'
need_DEG$type  <- as.factor(ifelse(need_DEG$Q_value < padj_cutoff & abs(need_DEG$log2FC) > log2FC_cutoff,
                                   ifelse(need_DEG$log2FC > log2FC_cutoff ,'Upregulate','Downregulate'),'NonDE'))

#得到HGNC基因信息
# library(biomaRt)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# gene_symbol<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), 
#                     filters= 'ensembl_gene_id', values = need_DEG$gene_id, mart = ensembl)

#need_DEG_merge_T1C1 <- merge(need_DEG, gene_symbol, by.x = c("gene_id"), by.y = c("ensembl_gene_id"), all=FALSE, all.x=TRUE, all.y=FALSE)

#改变列序
need_DEG_T1C1 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T2C2 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T3C2 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T4C2 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T5C2 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]
#need_DEG_T6C2 <- need_DEG[,c("sample_pair", "transcript_id", "log2FC", "p_value", "Q_value", "type")]

# 整合总表
need_DEG_merge_all <- do.call("rbind", list(need_DEG_T1C1, need_DEG_T2C2, need_DEG_T3C2, need_DEG_T4C2, need_DEG_T5C2, need_DEG_T6C2))
write.table(need_DEG_merge_all, "Different.expressed.transcripts.9.23.txt", sep = "\t", row.names = FALSE, quote = FALSE)

############################################################################################
#结果做火山图
title <- paste0('T1/C1', 
                ' Up :  ',nrow(need_DEG_T1C1[need_DEG_T1C1$type =='Upregulate',]) ,
                '\n Down : ',nrow(need_DEG_T1C1[need_DEG_T1C1$type =='Downregulate',]),
                '\n FoldChange >= ',round(2^log2FC_cutoff,3))

g <- ggplot(data=need_DEG_T1C1, 
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

ggsave(g,filename = 'volcano_padj_T6C2.pdf',width =8,height =7.5)
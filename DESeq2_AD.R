setwd()
getwd()

#Load packages
library(DESeq2)
library(dplyr)
library(reshape)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

############# Healthy vs MCI #############

coldata.1=read.table("coldata.HvsM.txt", header=T, row.names=1)
cts.1=read.table("rawcounts.HvsM.txt", header=T, row.names = 1)
head(cts.1)

all(rownames(coldata.1) %in% colnames(cts.1))

dds.1 <- DESeqDataSetFromMatrix(countData = cts.1,
colData = coldata.1,
design = ~ condition)
dds.1

keep.1 <- rowSums(counts(dds.1)) >= 5
dds.1 <- dds.1[keep.1,]


dds.1$condition <- factor(dds.1$condition, levels = c("Healthy","Mild_cog"))
dds.1 <- DESeq(dds.1)
res.1 <- results(dds.1)

res_HvsM <- results(dds.1, contrast=c("condition","Healthy","Mild_cog"))
resultsNames(dds.1)


resOrdered <- res_HvsM[order(res_HvsM$pvalue),]
summary(res.1)

res.1 <- results(dds.1, alpha=0.05)
summary(res.1)


degs.HvsM <- as.data.frame(subset(res.1,(abs(res.1$log2FoldChange) >= 1) & (res.1$padj <= 0.05)))
nrow(degs.HvsM) # number of degs (differentially expressed genes)
write.csv(degs.HvsM, file="DEGs_HvsM_results.csv")

################### Healthy vs AD ###################


coldata.2=read.table("coldata.HvsAD.txt", header=T, row.names=1)
cts.2=read.table("rawcounts.HvsAD.txt", header=T, row.names = 1)
head(cts.1)


all(rownames(coldata.2) %in% colnames(cts.2))


dds.2 <- DESeqDataSetFromMatrix(countData = cts.2,
                                colData = coldata.2,
                                design = ~ condition)
dds.2

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]


dds.2$condition <- factor(dds.2$condition, levels = c("Healthy","Alzheimers"))
dds.2 <- DESeq(dds.2)
res.2 <- results(dds.2)

res_HvsAD <- results(dds.2, contrast=c("condition","Healthy","Alzheimers"))
resultsNames(dds.2)


resOrdered <- res_HvsAD[order(res_HvsAD$pvalue),]
summary(res)

res.2 <- results(dds.2, alpha=0.05)
summary(res.2)


degs.HvsAD <- as.data.frame(subset(res.2,(abs(res.2$log2FoldChange) >= 2) & (res.2$padj <= 0.05)))
nrow(degs.HvsAD) # number of degs (differentially expressed genes)
write.csv(degs.HvsAD, file="DEGs_HvsAD_results.csv")


############## Combine res of H vs MCI and H vs AD #############

#Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 2
threshold.1 <- res.1$padj < padj.cutoff & abs(res.1$log2FoldChange) > lfc.cutoff
threshold.2 <- res.2$padj < padj.cutoff & abs(res.2$log2FoldChange) > lfc.cutoff

length(which(threshold.1))
length(which(threshold.2))

res.1$threshold_HvsM <- threshold.1
res.2$threshold_HvsAD <- threshold.2

sig.HvsM <- data.frame(subset(res.1, threshold_HvsM==TRUE))
sig.HvsAD <- data.frame(subset(res.2, threshold_HvsAD==TRUE))

sig.HvsM_ordered <- sig.HvsM[order(sig.HvsM$padj), ]
sig.HvsAD_ordered <- sig.HvsAD[order(sig.HvsAD$padj), ]

top20_sig.HvsM_genes <- rownames(sig.HvsM_ordered[1:20, ])
top20_sig.HvsAd_genes <- rownames(sig.HvsAD_ordered[1:20, ])


normalized_counts.1 <- counts(dds.1, normalized=T)
normalized_counts.2 <- counts(dds.2, normalized=T)

top20_sig.HvsM_norm <- as.data.frame(normalized_counts.1[top20_sig.HvsM_genes, ])
top20_sig.HvsAD_norm <- as.data.frame(normalized_counts.2[top20_sig.HvsAd_genes, ])


top20_sig.HvsM_norm$gene=rownames(top20_sig.HvsM_norm)
top20_sig.HvsAD_norm$gene=rownames(top20_sig.HvsAD_norm)

top20.sig.HvsM_AD=merge(top20_sig.HvsAD_norm, top20_sig.HvsM_norm, by=c("gene", "H1", "H2"), all=T)

#Missing values of genes from HvsM
A=top20.sig.HvsM_AD[is.na(top20.sig.HvsM_AD$AD1),]
A2=as.data.frame(normalized_counts.2[A$gene, ][,c(3:4)])
B=merge(top20.sig.HvsM_AD, A2, by=c("gene","AD1", "AD2"), all=T)
head(B)

library(dplyr)
coalesce_by_column <- function(df) {
return(coalesce(df[1], df[2]))
}


C=as.data.frame(B %>%
group_by(gene) %>%
summarise_all(coalesce_by_column))
head(C)

#Missing values of genes from HvsAD
a=top20.sig.HvsM_AD[is.na(top20.sig.HvsM_AD$MC1),]
nrow(a)
a2=as.data.frame(normalized_counts.1[a$gene, ][,c(3:4)])
a2$gene=rownames(a2)
b=merge(C, a2, by=c("gene","MC1", "MC2"), all=T)

c=as.data.frame(b %>%
group_by(gene) %>%
summarise_all(coalesce_by_column))
head(c)
c

top20.sig.HvsM_AD.final=c
top20.sig.HvsM_AD.final$gene->rownames(top20.sig.HvsM_AD.final)
top20.sig.HvsM_AD.final$gene=NULL
head(top20.sig.HvsM_AD.final)

#Heatmap
anno.col=read.table("coldata.txt", header=T, row.names=1)
anno.col


ppi=330

png("H_MCI_AD_logFC.png", width=5*ppi, height=8*ppi, res=ppi)
pheatmap(top20.sig.HvsM_AD.final,  cluster_rows = T, show_rownames=T,
border_color=NA, fontsize = 10, scale="row",
fontsize_row = 10, height=20, annotation_col = anno.col)
dev.off()

#Volcanoplot
res.HvsM_df <- data.frame(res.1)
head(res.HvsM_df)
res.HvsAD_df <- data.frame(res.2)
head(res.HvsAD_df)



res.HvsD_df_ordered <- res.HvsAD_df[order(res.HvsAD_df$padj), ]
res.HvsAD_df_ordered$genelabels <- rownames(res.HvsAD_df_ordered) %in% rownames(res.HvsAD_df_ordered[1:10,])
res.HvsAD_df_ordered <- res.HvsAD_df[order(res.HvsAD_df$padj), ]
res.HvsAD_df_ordered$genelabels <- rownames(res.HvsAD_df_ordered) %in% rownames(res.HvsAD_df_ordered[1:10,])


png("H_MCI_volcano.png", width=8*ppi, height=6*ppi, res=ppi)

ggplot(res.HvsM_df_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res.HvsM_df_ordered),""))) +
  ggtitle("Healthy vs Mild cognitive impairment") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +theme_classic()

dev.off()



png("H_AD_volcano.png", width=8*ppi, height=6*ppi, res=ppi)

ggplot(res.HvsAD_df_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res.HvsAD_df_ordered),""))) +
  ggtitle("Healthy vs Alzheimer's Disease") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + theme_classic()

dev.off()





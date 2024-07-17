
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(DESeq2)
BiocManager::install("apeglm")

setwd("C:/Users/aurin/Desktop/tsg2")
getwd()

## Preparing the data

counts_nb2<-read.table(file = "Galaxy40_nb2.tabular", header = TRUE, sep = "\t")
colnames(counts_nb2)[2]='nb2'

counts_nb3<-read.table(file = "Galaxy38_nb3.tabular", header = TRUE, sep = "\t")
colnames(counts_nb3)[2]='nb3'

counts_03<-read.table(file = "Galaxy34_pa03.tabular", header = TRUE, sep = "\t")
colnames(counts_03)[2]='03'

counts_05<-read.table(file = "Galaxy36_pa05.tabular", header = TRUE, sep = "\t")
colnames(counts_05)[2]='05'

counts_matrix <- counts_nb2 %>% 
  left_join(counts_nb3, by=c('Geneid')) %>%
  left_join(counts_03, by=c('Geneid')) %>%
  left_join(counts_05, by=c('Geneid'))


counts_matrix <- counts_matrix %>% column_to_rownames(., var = 'Geneid')

coldata <- data.frame(cells=c('nb2', 'nb3', 'pa03', 'pa05'),
                      state=c('normal', 'normal', 'cancer', 'cancer')
)


## construct a deseq dataset object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, 
                              colData = coldata, 
                              design = ~ state)
dds


## prefiltering - keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

## set the factor level
dds$state <- relevel(dds$state, ref = "normal")

## run DESeq 
dds <- DESeq(dds)
res <- results(dds)
res

summary(res)

## MA plot of the results
plotMA(res)

# shrinkage of low counts
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")	
plotMA(resLFC, ylim=c(-2,2))

## converting res object to a df
res_df <- as.data.frame(res)
write.csv(res_df, 'res_df.csv')
sig <- res_df[which(res_df$padj <= 0.05),] 
up <- sig[which(sig$log2FoldChange > 0),]
down <- sig[which(sig$log2FoldChange < 0),]

write.csv(sig, 'sig_diffrential_exp.csv')
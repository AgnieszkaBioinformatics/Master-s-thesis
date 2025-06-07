library(DESeq2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)

## differentially bounded sites: chipseeker annotations
ann <- "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\chipseeker\\"
# lhp1_3h1 vs lhp1
# lhp1 vs WT
# lhp1 vs 3h1


## RNA_seq results -> profiling on genes from chipseq
count <- counts(dds)
count <- as.data.frame(count)

size_sum(count)   # "[26,137 × 12]"
## read in lhp1

########### lhp1_3h1 vs lhp1
#select counts for genes differentially expressed in lhp1_3h1 vs lhp1
select_lhp1_3h1_vs_lhp1 <- subset(count, rownames(count) %in% lhp1_3h1_vs_lhp1$geneId)    # "[140 × 12]"
selected_genes1 <- rownames(select_lhp1_3h1_vs_lhp1)

df <- as.data.frame(colData(dds)[,c("samples", "genotype")])


pheatmap(assay(rlog(dds))[selected_genes1,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, show_colnames = FALSE, annotation_col=df, scale="row")



########## lhp1 vs WT
select_lhp1_vs_wt<- subset(count, rownames(count) %in% lhp1_vs_WT$geneId)    # "[140 × 12]"
selected_genes2 <- rownames(select_lhp1_vs_wt)

pheatmap(assay(rlog(dds))[selected_genes2,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, show_colnames = FALSE, annotation_col=df, scale="row")



########## lhp1 vs 3h1
select_lhp1_vs_3h1<- subset(count, rownames(count) %in% lhp1_vs_3h1$geneId)    # "[140 × 12]"
selected_genes3 <- rownames(select_lhp1_vs_3h1)

pheatmap(assay(rlog(dds))[selected_genes3,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, show_colnames = FALSE, annotation_col=df, scale="row")


########## X3h1_vs_WT
select_X3h1_vs_WT<- subset(count, rownames(count) %in% X3h1_vs_WT$geneId)    # "[140 × 12]"
selected_genes4 <- rownames(select_X3h1_vs_WT)

pheatmap(assay(rlog(dds))[selected_genes4,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, show_colnames = FALSE, annotation_col=df, scale="row")


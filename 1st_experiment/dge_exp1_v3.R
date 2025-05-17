library(DESeq2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)

alpha <- 0.01
lfcThreshold <- 1

                      ### reading in the data and preparing for DESEQ ###
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\counts_v2\\counts_gene"
setwd(outdir)
samples <- c('3h1_2_sorted_s', '3h1_3_sorted_s', '3h1_4_sorted_s', 
             'Ihp1_1_sorted_s', 'Ihp1_3_sorted_s', 'Ihp1_4_sorted_s',
             'Ihp1_3h1_2_sorted_s', 'Ihp1_3h1_3_sorted_s', 'Ihp1_3h1_4_sorted_s',
             'wt_1_sorted_s', 'wt_3_sorted_s', 'wt_4_sorted_s')


for (sample in samples) {
  file_dir <- file.path(paste0(sample, ".txt"))
  df <- read.table(file_dir, header = TRUE, sep = "\t")
  
  df_long <- df %>%
    separate_rows(Chr, Start, End, Strand, sep = ";") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  name <- paste0("df_", sample)
  assign(name, df_long)
}

df_list <- c('df_3h1_2_sorted_s', 'df_3h1_3_sorted_s', 'df_3h1_4_sorted_s', 
             'df_Ihp1_1_sorted_s', 'df_Ihp1_3_sorted_s', 'df_Ihp1_4_sorted_s',
             'df_Ihp1_3h1_2_sorted_s', 'df_Ihp1_3h1_3_sorted_s', 'df_Ihp1_3h1_4_sorted_s',
             'df_wt_1_sorted_s', 'df_wt_3_sorted_s', 'df_wt_4_sorted_s')

for (df_name in df_list) {
  df <- get(df_name)
  
  # Remove columns 2 to 6
  df <- df[, -c(2:6)]
  
  # Remove duplicates based on the 'Geneid' column
  df <- df[!duplicated(df$Geneid), ]
  
  # Assign the modified data frame back to the original name
  assign(df_name, df)
}


combined <- df_3h1_2_sorted_s

for (df in df_list) {
  df <- get(df)
  combined <- merge(x = combined, y = df)
}

rownames(combined) <- combined$Geneid
combined$Geneid <- NULL

# saving
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\3h1_lhp1\\dge\\v3"
setwd(outdir)
write.csv(combined, 'combined_3h1_lhp1.csv')

combined <- read.csv('combined_3h1_lhp1.csv')
rownames(combined) <- combined$X
combined$X <- NULL



                          ##### PCA on all samples ######
colors_sample <- data.frame(samples = colnames(combined),
                            names = c("3h1", "3h1", "3h1",
                                      "lhp1", "lhp1", "lhp1",
                                      "lhp1_3h1", "lhp1_3h1", "lhp1_3h1",
                                      "wt", "wt", "wt"))

pca_result <- prcomp(t(combined))  # Transpose to have samples as rows

# Calculate the percentage of variance explained by each principal component
percent_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create a data frame of PCA results
pca_df <- as.data.frame(pca_result$x[, 1:2])  # Select the first two principal components
pca_df$samples <- rownames(pca_df)

# Merge with colors_sample data frame to get color groups
pca_df <- merge(pca_df, colors_sample, by = "samples")

# Plot PCA with ggplot2   -----------------> coś tu nie działa
ggplot(pca_df, aes(x = PC1, y = PC2, color = names)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of 3h1, lhp1, lhp1_3h1, and wt",
       x = sprintf("Principal Component 1 (%.2f%%)", percent_variance[1]),
       y = sprintf("Principal Component 2 (%.2f%%)", percent_variance[2])) +
  theme_minimal() +
  scale_color_manual(values = c("3h1" = "blue", "lhp1" = "red",
                                "lhp1_3h1" = "green", "wt" = "purple"))



###########


                                  ####### DEG #######


coldata_3h1 <- data.frame(cells = colnames(combined[,c(1:3, 10:12)]),
                          state = c("3h1", "3h1", "3h1",
                                    "WT", "WT", "WT"))

coldata_Ihp1 <- data.frame(cells = colnames(combined[,c(4:6, 10:12)]),
                           state = c("Ihp1", "Ihp1", "Ihp1",
                                     "WT", "WT", "WT"))

coldata_Ihp1_3h1 <- data.frame(cells = colnames(combined[,c(7:9, 10:12)]),
                               state = c("Ihp1_3h1", "Ihp1_3h1", "Ihp1_3h1",
                                         "WT", "WT", "WT"))

rownames(coldata_3h1) <- coldata_3h1$cells

rownames(coldata_Ihp1) <- coldata_Ihp1$cells

rownames(coldata_Ihp1_3h1) <- coldata_Ihp1_3h1$cells

coldata_all <- data.frame(cells = colnames(combined),
                          samples = c("3h1_2", "3h1_3", "3h1_4",
                                      "Ihp1_1", "Ihp1_3", "Ihp1_4",
                                      "Ihp1_3h1_2", "Ihp1_3h1_3", "Ihp1_3h1_4",
                                      "WT_1", "WT_3", "WT_4"),
                          
                          genotype = c("3h1", "3h1", "3h1",
                                    "Ihp1", "Ihp1", "Ihp1",
                                    "Ihp1_3h1", "Ihp1_3h1", "Ihp1_3h1",
                                    "WT", "WT", "WT"),
                          
                          replicate = c("2", "3", "4", 
                                        "1", "3", "4",
                                        "2", "3", "4",
                                        "1", "3", "4"))


rownames(coldata_all) <- coldata_all$cells

## DESeq object 

dds_3h1 <- DESeqDataSetFromMatrix(countData = combined[,c(1:3, 10:12)],
                                  colData = coldata_3h1,
                                  design = ~ state)

dds_Ihp1 <- DESeqDataSetFromMatrix(countData = combined[,c(4:6, 10:12)],
                                   colData = coldata_Ihp1,
                                   design = ~ state)

dds_Ihp1_3h1 <- DESeqDataSetFromMatrix(countData = combined[,c(7:9, 10:12)],
                                       colData = coldata_Ihp1_3h1,
                                       design = ~ state)

# on all data
dds <- DESeqDataSetFromMatrix(countData = combined,
                              colData = coldata_all,
                              design = ~ genotype)


### filtering out low counts ###
keep_3h1 <- rowSums(counts(dds_3h1)) >= 10
dds_3h1 <- dds_3h1[keep_3h1,]

keep_Ihp1 <- rowSums(counts(dds_Ihp1)) >= 10
dds_Ihp1 <- dds_Ihp1[keep_Ihp1,]

keep_Ihp1_3h1 <- rowSums(counts(dds_Ihp1_3h1)) >= 10
dds_Ihp1_3h1 <- dds_Ihp1_3h1[keep_Ihp1_3h1,]

keep_dds <- rowSums(counts(dds)) >= 10
dds <- dds[keep_dds,]

### setting the reference ###
dds_3h1$state <- relevel(dds_3h1$state, ref = "WT")
dds_Ihp1$state <- relevel(dds_Ihp1$state, ref = "WT")
dds_Ihp1_3h1$state <- relevel(dds_Ihp1_3h1$state, ref = "WT")

dds$genotype <- relevel(dds$genotype, ref = "WT")

### running the analysis
dds_3h1 <- DESeq(dds_3h1)
res_3h1 <- results(dds_3h1, alpha = alpha, lfcThreshold = lfcThreshold)

dds_Ihp1 <- DESeq(dds_Ihp1)
res_Ihp1 <- results(dds_Ihp1, alpha = alpha, lfcThreshold = lfcThreshold)

dds_Ihp1_3h1 <- DESeq(dds_Ihp1_3h1)
res_Ihp1_3h1 <- results(dds_Ihp1_3h1, alpha = alpha, lfcThreshold = lfcThreshold)

dds <- DESeq(dds) 
res_dds <- results(dds, alpha = alpha, lfcThreshold = lfcThreshold)

# plotting (only single mutants)
plotMA(res_3h1)
rld_3h1 <- rlog(dds_3h1, blind=FALSE)
plotPCA(rld_3h1, intgroup=c("state"), ntop = 1000)

plotMA(res_Ihp1)
rld_Ihp1 <- rlog(dds_Ihp1, blind=FALSE)
plotPCA(rld_Ihp1, intgroup=c("state"), ntop = 1000)

plotMA(res_Ihp1_3h1)
rld_Ihp1_3h1 <- rlog(dds_Ihp1_3h1, blind = FALSE)
plotPCA(rld_Ihp1_3h1, intgroup=c("state"), ntop = 1000)


# plot all data
rld_all <- rlog(dds, blind = FALSE)
plotPCA(rld_all, intgroup=c("genotype"), ntop = 1000)

#plotCounts(dds, gene="AT1G06760", intgroup="genotype") -> nic odkrywczego, lepiej heatmape


# saving results

res_df_3h1 <- as.data.frame(res_3h1)
write.csv(res_df_3h1, 'res_df_3h1.csv')
sig_3h1 <- res_df_3h1[which(res_df_3h1$padj <= 0.05),] 
write.csv(sig_3h1, 'sig_de_3h1.csv')

res_df_Ihp1 <- as.data.frame(res_Ihp1)
write.csv(res_df_Ihp1, 'res_df_Ihp1.csv')
sig_Ihp1 <- res_df_Ihp1[which(res_df_Ihp1$padj <= 0.05),] 
write.csv(sig_Ihp1, 'sig_de_Ihp1.csv')

res_df_Ihp1_3h1 <- as.data.frame(res_Ihp1_3h1)
write.csv(res_df_Ihp1_3h1, 'res_df_Ihp1_3h1.csv')
sig_Ihp1_3h1 <- res_df_Ihp1_3h1[which(res_df_Ihp1_3h1$padj <= 0.05),] 
write.csv(sig_Ihp1_3h1, 'sig_de_Ihp1_3h1.csv')




                        ######## lhp1 profile clustering #########

count <- counts(dds)
count <- as.data.frame(count)

size_sum(count)   # "[26,137 × 12]"
## read in lhp1
sig_lhp1 <- read.csv("sig_de_Ihp1.csv")

### select counts for genes differentially expressed in lhp1
select <- subset(count, rownames(count) %in% sig_lhp1$X)    # "[140 × 12]"
selected_genes <- rownames(select)

df <- as.data.frame(colData(dds)[,c("samples", "genotype")])

pheatmap(assay(rlog(dds))[selected_genes,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, show_colnames = FALSE, annotation_col=df, scale="row")







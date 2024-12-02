library(DESeq2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(regionReport)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

library(FactoMineR)
library(factoextra)


alpha <- 0.05

### reading in the data ###
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\exp2\\counts"
setwd(outdir)
samples <- c('BG9_sorted', 'BG18_sorted_s', 'BG32_sorted_s', 
             'brm1_31_sorted_s', 'brm1_37_sorted_s', 'brm1_38_sorted_s',
             'brm5_41_sorted_s', 'brm5_43_sorted_s', 'brm5_44_sorted_s',
             'Col4_sorted_s', 'Col19_sorted_s', 'Col20_sorted_s',
             'KR_5_sorted_s', 'KR_8_sorted_s', 'KR_21_sorted_s',
             'LD_7_sorted_s', 'LD_22_sorted_s', 'LD_30_sorted_s')

markers <- c("AT2G46020", "AT4G32280", "AT4G19170", "AT3G01970", "AT1G69490", 
             "AT2G40080", "AT3G07650", "AT2G47260")


for (sample in samples) {
  file_dir <- file.path(paste0(sample, ".txt"))
  df <- read.table(file_dir, header = TRUE, sep = "\t")
  
  df_long <- df %>%
    separate_rows(Chr, Start, End, Strand, sep = ";") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  name <- paste0("df_", sample)
  assign(name, df_long)
}

df_list <- c('df_BG9_sorted', 'df_BG18_sorted_s', 'df_BG32_sorted_s', 
             'df_brm1_31_sorted_s', 'df_brm1_37_sorted_s', 'df_brm1_38_sorted_s',
             'df_brm5_41_sorted_s', 'df_brm5_43_sorted_s', 'df_brm5_44_sorted_s',
             
             'df_Col4_sorted_s', 'df_Col19_sorted_s', 'df_Col20_sorted_s',
             
             'df_KR_5_sorted_s', 'df_KR_8_sorted_s', 'df_KR_21_sorted_s',
             'df_LD_7_sorted_s', 'df_LD_22_sorted_s', 'df_LD_30_sorted_s')



for (df_name in df_list) {
  df <- get(df_name)
  
  # Remove columns 2 to 6
  df <- df[, -c(2:6)]
  
  # Remove duplicates based on the 'Geneid' column
  df <- df[!duplicated(df$Geneid), ]
  
  # Assign the modified data frame back to the original name
  assign(df_name, df)
}

combined <- df_BG9_sorted

for (df in df_list) {
  df <- get(df)
  combined <- merge(x = combined, y = df)
}

rownames(combined) <- combined$Geneid
combined$Geneid <- NULL

## PCA of col bg ld kr

selected <- combined[,c(1:3, 10:12, 13:15, 16:18)]

combined_data <- scale(selected)

colors_sample <- data.frame(samples = colnames(combined_data),
                            names = c("BG9", "BG18", "BG32",
                                      "Col4", "Col19", "Col20",
                                      "KR5", "KR8", "KR21",
                                      "LD7", "LD22", "LD30"),
                            
                            sample_number = c("BG", "BG", "BG",
                                     "Col", "Col", "Col",
                                     "KR", "KR", "KR",
                                     "LD", "LD", "LD"))

pca_result <- prcomp(t(combined_data))  # Transpose to have samples as rows

# Calculate the percentage of variance explained by each principal component
percent_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create a data frame of PCA results
pca_df <- as.data.frame(pca_result$x[, 1:2])  # Select the first two principal components
pca_df$samples <- rownames(pca_df)

# Merge with colors_sample data frame to get color groups
pca_df <- merge(pca_df, colors_sample, by = "samples")

# Plot PCA with ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = sample_number)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of BRM1, BRM5, BG, and Col",
       x = sprintf("Principal Component 1 (%.2f%%)", percent_variance[1]),
       y = sprintf("Principal Component 2 (%.2f%%)", percent_variance[2])) +
  theme_minimal() +
  scale_color_manual(values = c("BG" = "blue", "Col" = "red", "KR" = "green", "LD" = "purple"))



### preparing the input ###
rownames(merged_bg) <- merged_bg$Geneid
merged_bg$Geneid <- NULL

rownames(merged_brm1) <- merged_brm1$Geneid
merged_brm1$Geneid <- NULL

rownames(merged_brm5) <- merged_brm5$Geneid
merged_brm5$Geneid <- NULL

#rownames(merged_crt) <- merged_crt$Geneid
#merged_crt$Geneid <- NULL

rownames(merged_kr) <- merged_kr$Geneid
merged_kr$Geneid <- NULL

rownames(merged_ld) <- merged_ld$Geneid
merged_ld$Geneid <- NULL


## PCA of  brm1 brm5 bg col


combined_data <- scale(combined_data)

colors_sample <- data.frame(samples = colnames(combined_data),
                            names = c("BG", "BG", "BG",
                                      "BRM1", "BRM1", "BRM1",
                                      "BRM5", "BRM5", "BRM5",
                                      "Col", "Col", "Col"))

pca_result <- prcomp(t(combined_data))  # Transpose to have samples as rows

# Calculate the percentage of variance explained by each principal component
percent_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create a data frame of PCA results
pca_df <- as.data.frame(pca_result$x[, 1:2])  # Select the first two principal components
pca_df$samples <- rownames(pca_df)

# Merge with colors_sample data frame to get color groups
pca_df <- merge(pca_df, colors_sample, by = "samples")

# Plot PCA with ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = names)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of BRM1, BRM5, BG, and Col",
       x = sprintf("Principal Component 1 (%.2f%%)", percent_variance[1]),
       y = sprintf("Principal Component 2 (%.2f%%)", percent_variance[2])) +
  theme_minimal() +
  scale_color_manual(values = c("BRM1" = "blue", "BRM5" = "red", "BG" = "green", "Col" = "purple"))



# PCA of col bg ld kr
merged_bg$dataset <- NULL
merged_col$dataset <- NULL
merged_kr$BG9_sorted.bam <- NULL
merged_kr$BG18_sorted_s.bam <- NULL
merged_kr$BG32_sorted_s.bam <- NULL

combined_data2 <- cbind(merged_bg, merged_ld, merged_kr, merged_col) 
combined_data2$Geneid <- NULL
rownames(combined_data2) <- NULL

combined_data2 <- scale(combined_data2)

colors_sample2 <- data.frame(samples = colnames(combined_data2),
                            names = c("BG", "BG", "BG",
                                      "LD", "LD", "LD",
                                      "KR", "KR", "KR",
                                      "Col", "Col", "Col"))

pca_result2 <- prcomp(t(combined_data2))  # Transpose to have samples as rows

# Calculate the percentage of variance explained by each principal component
percent_variance2 <- pca_result2$sdev^2 / sum(pca_result2$sdev^2) * 100

# Create a data frame of PCA results
pca_df2 <- as.data.frame(pca_result2$x[, 1:2])  # Select the first two principal components
pca_df2$samples <- rownames(pca_df2)

# Merge with colors_sample data frame to get color groups
pca_df2 <- merge(pca_df2, colors_sample2, by = "samples")

# Plot PCA with ggplot2
ggplot(pca_df2, aes(x = PC1, y = PC2, color = names)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of LD, KR, BG, and Col",
       x = sprintf("Principal Component 1 (%.2f%%)", percent_variance2[1]),
       y = sprintf("Principal Component 2 (%.2f%%)", percent_variance2[2])) +
  theme_minimal() +
  scale_color_manual(values = c("LD" = "blue", "KR" = "red", "BG" = "green", "Col" = "purple"))




  

coldata_bg <- data.frame(cells = colnames(merged_bg),
                          state = c("BG", "BG", "BG",
                                    "Col", "Col", "Col"))

coldata_brm1 <- data.frame(cells = colnames(merged_brm1),
                           state = c("brm1", "brm1", "brm1",
                                     "LD", "LD", "LD"))

coldata_brm1 <- data.frame(cells = colnames(merged_brm1),
                           state = c("brm1", "brm1", "brm1",
                                     "KR", "KR", "KR"))

coldata_brm5 <- data.frame(cells = colnames(merged_brm5),
                               state = c("brm5", "brm5", "brm5",
                                         "normal", "normal", "normal"))

coldata_crt <- data.frame(cells = colnames(merged_crt),
                           state = c("crt", "crt", "crt",
                                     "normal", "normal", "normal"))

coldata_kr <- data.frame(cells = colnames(merged_kr),
                           state = c("KR", "KR", "KR",
                                     "BG", "BG", "BG"))

coldata_ld <- data.frame(cells = colnames(merged_ld),
                           state = c("LD", "LD", "LD",
                                     "BG", "BG", "BG"))

coldata_ld <- data.frame(cells = colnames(merged_ld),
                         state = c("LD", "LD", "LD",
                                   "KR", "KR", "KR"))



rownames(coldata_bg) <- coldata_bg$cells
coldata_bg$cells <- NULL
all(colnames(merged_bg) == rownames(coldata_bg))

rownames(coldata_brm1) <- coldata_brm1$cells
coldata_brm1$cells <- NULL
all(colnames(merged_brm1) == rownames(coldata_brm1))

rownames(coldata_brm5) <- coldata_brm5$cells
coldata_brm5$cells <- NULL
all(colnames(merged_brm5) == rownames(coldata_brm5))

rownames(coldata_crt) <- coldata_crt$cells
coldata_crt$cells <- NULL
all(colnames(merged_crt) == rownames(coldata_crt))

rownames(coldata_kr) <- coldata_kr$cells
coldata_kr$cells <- NULL
all(colnames(merged_kr) == rownames(coldata_kr))

rownames(coldata_ld) <- coldata_ld$cells
coldata_ld$cells <- NULL
all(colnames(merged_ld) == rownames(coldata_ld))



### DESeq object ###
dds_bg <- DESeqDataSetFromMatrix(countData = merged_bg,
                                  colData = coldata_bg,
                                  design = ~ state)

dds_brm1 <- DESeqDataSetFromMatrix(countData = merged_brm1,
                                   colData = coldata_brm1,
                                   design = ~ state)

dds_brm5 <- DESeqDataSetFromMatrix(countData = merged_brm5,
                                       colData = coldata_brm5,
                                       design = ~ state)

dds_crt <- DESeqDataSetFromMatrix(countData = merged_crt,
                                   colData = coldata_crt,
                                   design = ~ state)

dds_kr <- DESeqDataSetFromMatrix(countData = merged_kr,
                                   colData = coldata_kr,
                                   design = ~ state)

dds_ld <- DESeqDataSetFromMatrix(countData = merged_ld,
                                   colData = coldata_ld,
                                   design = ~ state)

### filtering out low counts ###
keep_bg <- rowSums(counts(dds_bg)) >= 10
dds_bg <- dds_bg[keep_bg,]

keep_brm1 <- rowSums(counts(dds_brm1)) >= 10
dds_brm1 <- dds_brm1[keep_brm1,]

keep_brm5 <- rowSums(counts(dds_brm5)) >= 10
dds_brm5 <- dds_brm5[keep_brm5,]

keep_crt <- rowSums(counts(dds_crt)) >= 10
dds_crt <- dds_crt[keep_crt,]

keep_kr <- rowSums(counts(dds_kr)) >= 10
dds_kr <- dds_kr[keep_kr,]

keep_ld <- rowSums(counts(dds_ld)) >= 10
dds_ld <- dds_ld[keep_ld,]


### stting the reference ###
dds_bg$state <- relevel(dds_bg$state, ref = "Col")
dds_brm1$state <- relevel(dds_brm1$state, ref = "KR")
dds_brm5$state <- relevel(dds_brm5$state, ref = "normal")

dds_crt$state <- relevel(dds_crt$state, ref = "normal")
dds_kr$state <- relevel(dds_kr$state, ref = "BG")
dds_ld$state <- relevel(dds_ld$state, ref = "BG")
dds_ld$state <- relevel(dds_ld$state, ref = "KR")

### running the analysis
dds_bg <- DESeq(dds_bg)
res_bg <- results(dds_bg, alpha = alpha)

dds_brm1 <- DESeq(dds_brm1)
res_brm1 <- results(dds_brm1, alpha = alpha)

dds_brm5 <- DESeq(dds_brm5)
res_brm5 <- results(dds_brm5, alpha = alpha)

dds_crt <- DESeq(dds_crt)
res_crt <- results(dds_crt, alpha = alpha)

dds_kr <- DESeq(dds_kr)
res_kr <- results(dds_kr, alpha = alpha)

dds_ld <- DESeq(dds_ld)
res_ld <- results(dds_ld, alpha = alpha)

### visualizing the results ###
plotMA(res_bg)
rld <- rlog(dds_bg, blind=FALSE)
plotPCA(rld, intgroup=c("state"), ntop = 1000)

plotMA(res_brm1)
rld_brm1 <- rlog(dds_brm1, blind=FALSE)
plotPCA(rld_brm1, intgroup=c("state"), ntop = 1000)

plotMA(res_brm5)
rld_brm5 <- rlog(dds_brm5, blind = FALSE)
plotPCA(rld_brm5, intgroup=c("state"), ntop = 1000)

plotMA(res_crt)
rld_crt <- rlog(dds_crt, blind=FALSE)
plotPCA(rld_crt, intgroup=c("cells"), ntop = 1000)

plotMA(res_kr)
rld_kr <- rlog(dds_kr, blind=FALSE)
plotPCA(rld_kr, intgroup=c("state"), ntop = 1000)

plotMA(res_ld)
rld_ld <- rlog(dds_ld, blind = FALSE)
plotPCA(rld_ld, intgroup=c("state"), ntop = 1000)


### generating reports ###

report_bg <- DESeq2Report(dds_bg, project = "DESeq2 BG report", 
                           intgroup = "state", outdir = outdir, output = "index")

report_brm1 <- DESeq2Report(dds_brm1, project = "DESeq2 brm1 report", 
                            intgroup = "state", outdir = outdir, output = "index")

report_brm5 <- DESeq2Report(dds_brm5, project = "DESeq2 brm5 report", 
                                intgroup = "state", outdir = outdir, output = "index")


report_crt <- DESeq2Report(dds_crt, project = "DESeq2 crt report", 
                           intgroup = "cells", outdir = outdir, output = "index")

report_kr <- DESeq2Report(dds_kr, project = "DESeq2 KR report", 
                            intgroup = "state", outdir = outdir, output = "index")

report_ld <- DESeq2Report(dds_ld, project = "DESeq2 LD report", 
                                intgroup = "state", outdir = outdir, output = "index")


### saving
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\exp2"
setwd(outdir)

res_bg <- as.data.frame(res_bg)
write.csv(res_bg, 'res_df_bg_vs_col.csv')
sig_s <- res_bg[which(res_bg$padj <= 0.05),] 
write.csv(sig_s, 'sig_de_bg_vs_col.csv')


res_brm1 <- as.data.frame(res_brm1)
write.csv(res_brm1, 'res_df_brm1_vs_kr.csv')
sig_brm1 <- res_brm1[which(res_brm1$padj <= 0.05),] 
write.csv(sig_brm1, 'sig_de_brm1_vs_kr.csv')

res_brm5 <- as.data.frame(res_brm5)
write.csv(res_brm5, 'res_brm5.csv')
sig_brm5 <- res_brm5[which(res_brm5$padj <= 0.05),] 
write.csv(sig_brm5, 'sig_de_brm5.csv')

res_crt <- as.data.frame(res_crt)
write.csv(res_crt, 'res_df_crt.csv')
sig_crt <- res_crt[which(res_crt$padj <= 0.05),] 
write.csv(sig_crt, 'sig_de_crt.csv')

res_kr <- as.data.frame(res_kr)
write.csv(res_kr, 'res_df_kr_vs_bg.csv')
sig_kr <- res_kr[which(res_kr$padj <= 0.05),] 
write.csv(sig_kr, 'sig_de_kr_vs_bg.csv')

res_ld <- as.data.frame(res_ld)
write.csv(res_ld, 'res_ld_vs_bg.csv')
sig_ld <- res_ld[which(res_ld$padj <= 0.05),] 
write.csv(sig_ld, 'sig_de_ld_vs_bg.csv')

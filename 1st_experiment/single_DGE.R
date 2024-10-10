library(DESeq2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(regionReport)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

alpha <- 0.05

### reading in the data ###
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\counts_v2"
setwd(outdir)
samples <- c('3h1_2_s', '3h1_3_s', '3h1_4_s', 
             'Ihp1_1_s', 'Ihp1_3_s', 'Ihp1_4_s',
             'Ihp1_3h1_2_s', 'Ihp1_3h1_3_s', 'Ihp1_3h1_4_s',
             'wt_1_s', 'wt_3_s', 'wt_4_s')

for (sample in samples) {
  file_dir <- file.path(paste0(sample, ".txt"))
  df <- read.table(file_dir, header = TRUE, sep = "\t")
  
  df_long <- df %>%
    separate_rows(Chr, Start, End, Strand, sep = ";") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  name <- paste0("df_", sample)
  assign(name, df_long)
}

df_list <- c('df_3h1_2_s', 'df_3h1_3_s', 'df_3h1_4_s', 
             'df_Ihp1_1_s', 'df_Ihp1_3_s', 'df_Ihp1_4_s',
             'df_Ihp1_3h1_2_s', 'df_Ihp1_3h1_3_s', 'df_Ihp1_3h1_4_s',
             'df_wt_1_s', 'df_wt_3_s', 'df_wt_4_s')


for (df_name in df_list) {
  df <- get(df_name)
  
  # Remove columns 2 to 6
  df <- df[, -c(2:6)]
  
  # Remove duplicates based on the 'Geneid' column
  df <- df[!duplicated(df$Geneid), ]
  
  # Assign the modified data frame back to the original name
  assign(df_name, df)
}


### merging the data ###
merged_3h1 <- df_3h1_2_s
merged_3h1 <- merge(x = merged_3h1, y = df_wt_4_s)

merged_Ihp1 <- df_Ihp1_1_s
merged_Ihp1 <- merge(x = merged_Ihp1, y=df_wt_4_s)

merged_Ihp1_3h1 <- df_Ihp1_3h1_2_s
merged_Ihp1_3h1 <- merge(x=merged_Ihp1_3h1, y=df_wt_4_s)


### preparing the input ###
rownames(merged_3h1) <- merged_3h1$Geneid
merged_3h1$Geneid <- NULL

rownames(merged_Ihp1) <- merged_Ihp1$Geneid
merged_Ihp1$Geneid <- NULL

rownames(merged_Ihp1_3h1) <- merged_Ihp1_3h1$Geneid
merged_Ihp1_3h1$Geneid <- NULL


coldata_3h1 <- data.frame(cells = colnames(merged_3h1),
                           state = c("3h1", "3h1", "3h1",
                                     "normal", "normal", "normal"))

coldata_Ihp1 <- data.frame(cells = colnames(merged_Ihp1),
                           state = c("Ihp1", "Ihp1", "Ihp1",
                                     "normal", "normal", "normal"))

coldata_Ihp1_3h1 <- data.frame(cells = colnames(merged_Ihp1_3h1),
                               state = c("Ihp1_3h1", "Ihp1_3h1", "Ihp1_3h1",
                                         "normal", "normal", "normal"))

rownames(coldata_3h1) <- coldata_3h1$cells
coldata_3h1$cells <- NULL
all(colnames(merged_3h1) == rownames(coldata_3h1))

rownames(coldata_Ihp1) <- coldata_Ihp1$cells
coldata_Ihp1$cells <- NULL
all(colnames(merged_Ihp1) == rownames(coldata_Ihp1))

rownames(coldata_Ihp1_3h1) <- coldata_Ihp1_3h1$cells
coldata_Ihp1_3h1$cells <- NULL
all(colnames(merged_Ihp1_3h1) == rownames(coldata_Ihp1_3h1))


### DESeq object ###
dds_3h1 <- DESeqDataSetFromMatrix(countData = merged_3h1,
                                    colData = coldata_3h1,
                                    design = ~ state)

dds_Ihp1 <- DESeqDataSetFromMatrix(countData = merged_Ihp1,
                                  colData = coldata_Ihp1,
                                  design = ~ state)

dds_Ihp1_3h1 <- DESeqDataSetFromMatrix(countData = merged_Ihp1_3h1,
                                       colData = coldata_Ihp1_3h1,
                                       design = ~ state)

### filtering out low counts ###
keep_3h1 <- rowSums(counts(dds_3h1)) >= 10
dds_3h1 <- dds_3h1[keep_3h1,]

keep_Ihp1 <- rowSums(counts(dds_Ihp1)) >= 10
dds_Ihp1 <- dds_Ihp1[keep_Ihp1,]

keep_Ihp1_3h1 <- rowSums(counts(dds_Ihp1_3h1)) >= 10
dds_Ihp1_3h1 <- dds_Ihp1_3h1[keep_Ihp1_3h1,]

### stting the reference ###
dds_3h1$state <- relevel(dds_3h1$state, ref = "normal")
dds_Ihp1$state <- relevel(dds_Ihp1$state, ref = "normal")
dds_Ihp1_3h1$state <- relevel(dds_Ihp1_3h1$state, ref = "normal")

### running the analysis
dds_3h1 <- DESeq(dds_3h1)
res_3h1 <- results(dds_3h1, alpha = alpha)

dds_Ihp1 <- DESeq(dds_Ihp1)
res_Ihp1 <- results(dds_Ihp1, alpha = alpha)

dds_Ihp1_3h1 <- DESeq(dds_Ihp1_3h1)
res_Ihp1_3h1 <- results(dds_Ihp1_3h1, alpha = alpha)

### visualizing the results ###
plotMA(res_3h1)
rld <- rlog(dds_3h1, blind=FALSE)
plotPCA(rld, intgroup=c("state"), ntop = 1000)

plotMA(res_Ihp1)
rld_Ihp1 <- rlog(dds_Ihp1, blind=FALSE)
plotPCA(rld_Ihp1, intgroup=c("state"), ntop = 1000)

plotMA(res_Ihp1_3h1)
rld_Ihp1_3h1 <- rlog(dds_Ihp1_3h1, blind = FALSE)
plotPCA(rld_Ihp1_3h1, intgroup=c("state"), ntop = 1000)


### generating reports ###

report_3h1 <- DESeq2Report(dds_3h1, project = "DESeq2 3h1 report", 
                         intgroup = "state", outdir = outdir, output = "index")

report_Ihp1 <- DESeq2Report(dds_Ihp1, project = "DESeq2 Ihp1 report", 
                           intgroup = "state", outdir = outdir, output = "index")

report_Ihp1_3h1 <- DESeq2Report(dds_Ihp1_3h1, project = "DESeq2 Ihp1_3h1 report", 
                           intgroup = "state", outdir = outdir, output = "index")


### saving
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\deseq_v2\\short\\"
setwd(outdir)

res_df_3h1 <- as.data.frame(res_3h1)
write.csv(res_df_3h1, 'res_df_3h1.csv')
sig_s <- res_df_3h1[which(res_df_3h1$padj <= 0.05),] 
write.csv(sig_s, 'sig_de_3h1.csv')

res_df_Ihp1 <- as.data.frame(res_Ihp1)
write.csv(res_df_Ihp1, 'res_df_Ihp1.csv')
sig_Ihp1 <- res_df_Ihp1[which(res_df_Ihp1$padj <= 0.05),] 
write.csv(sig_Ihp1, 'sig_de_Ihp1.csv')

res_df_Ihp1_3h1 <- as.data.frame(res_Ihp1_3h1)
write.csv(res_df_Ihp1_3h1, 'res_df_Ihp1_3h1.csv')
sig_Ihp1_3h1 <- res_df_Ihp1_3h1[which(res_df_Ihp1_3h1$padj <= 0.05),] 
write.csv(sig_Ihp1_3h1, 'sig_de_Ihp1_3h1.csv')

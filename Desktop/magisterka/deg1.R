library(DESeq2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(regionReport)
library(ggplot2)

### reading in the data ###

DIR <- "C:\\Users\\aurin\\Desktop\\magisterka\\counts_v2"
samples <- c('3h1_2_l', '3h1_3_l', '3h1_4_l', 
             'Ihp1_1_l', 'Ihp1_3_l', 'Ihp1_4_l',
             'Ihp1_3h1_2_l', 'Ihp1_3h1_3_l', 'Ihp1_3h1_4_l',
             'wt_1_l', 'wt_3_l', 'wt_4_l',
             '3h1_2_s', '3h1_3_s', '3h1_4_s', 
             'Ihp1_1_s', 'Ihp1_3_s', 'Ihp1_4_s',
             'Ihp1_3h1_2_s', 'Ihp1_3h1_3_s', 'Ihp1_3h1_4_s',
             'wt_1_s', 'wt_3_s', 'wt_4_s')


for (sample in samples) {
  file_dir <- file.path(DIR, paste0(sample, ".txt"))
  df <- read.table(file_dir, header = TRUE, sep = "\t")
  
  df_long <- df %>%
    separate_rows(Chr, Start, End, Strand, sep = ";") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  name <- paste0("df_", sample)
  assign(name, df_long)
}

df_list <- list('df_3h1_2_l','df_3h1_3_l', 'df_3h1_4_l', 
                'df_Ihp1_1_l', 'df_Ihp1_3_l', 'df_Ihp1_4_l',
                'df_Ihp1_3h1_2_l', 'df_Ihp1_3h1_3_l', 'df_Ihp1_3h1_4_l',
                'df_wt_1_l', 'df_wt_3_l', 'df_wt_4_l',
                'df_3h1_2_s', 'df_3h1_3_s', 'df_3h1_4_s', 
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

# long
merged_df <- df_3h1_2_l
merged_df <- merge(x = merged_df, y = df_wt_4_l)
write.table(merged_df, file="merged_df_long.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#short
merged_df_s <- df_3h1_2_s
merged_df_s <- merge(x=merged_df_s, y = df_wt_4_s)


### preparing dataframes as input to Deseq ###
#long
rownames(merged_df) <- merged_df$Geneid
merged_df$Geneid <- NULL

coldata_long <- data.frame(cells = colnames(merged_df),
                           state = c("3h1", "3h1", "3h1",
                                     "Ihp1", "Ihp1", "Ihp1",
                                     "Ihp1_3h1", "Ihp1_3h1", "Ihp1_3h1",
                                     "normal", "normal", "normal"))

rownames(coldata_long) <- coldata_long$cells
coldata_long$cells <- NULL
all(colnames(merged_df) == rownames(coldata_long))


#short
rownames(merged_df_s) <- merged_df_s$Geneid
merged_df_s$Geneid <- NULL

coldata_short <- data.frame(cells = colnames(merged_df_s),
                           state = c("3h1", "3h1", "3h1",
                                     "Ihp1", "Ihp1", "Ihp1",
                                     "Ihp1_3h1", "Ihp1_3h1", "Ihp1_3h1",
                                     "normal", "normal", "normal"))

rownames(coldata_short) <- coldata_short$cells
coldata_short$cells <- NULL
all(colnames(merged_df_s) == rownames(coldata_short))


### Deseq object ###
#long
dds_long <- DESeqDataSetFromMatrix(countData = merged_df,
                       colData = coldata_long,
                       design = ~ state) #design factor


#short
dds_short <- DESeqDataSetFromMatrix(countData = merged_df_s,
                                   colData = coldata_short,
                                   design = ~ state)

### filtering out low counts ###
#long
keep_l <- rowSums(counts(dds_long)) >= 10
dds_long <- dds_long[keep_l,]

#short
keep_s <- rowSums(counts(dds_short)) >= 10
dds_short <- dds_short[keep_s,]

## normal = wt as a reference
# long
dds_long$state <- relevel(dds_long$state, ref = "normal")

#short
dds_short$state <- relevel(dds_short$state, ref = "normal")

### running the analysis ###
#long
dds_long <- DESeq(dds_long)
res_long <- results(dds_long)

#short
dds_short <- DESeq(dds_short)
res_short <- results(dds_short)

## converting res object to a df
# long
res_df_l <- as.data.frame(res_long)
write.csv(res_df_l, 'res_df_l.csv')
sig_l <- res_df_l[which(res_df_l$padj <= 0.05),] 
up_l <- sig_l[which(sig_l$log2FoldChange > 0),]
down_l <- sig_l[which(sig_l$log2FoldChange < 0),]

write.csv(sig_l, 'sig_de_long.csv')

#short
res_df_s <- as.data.frame(res_short)
write.csv(res_df_s, 'res_df_s.csv')
sig_s <- res_df_s[which(res_df_s$padj <= 0.05),] 
up_s <- sig_s[which(sig_s$log2FoldChange > 0),]
down_s <- sig_s[which(sig_s$log2FoldChange < 0),]

write.csv(sig_s, 'sig_de_short.csv')

### visualizing the results ###

summary(res_long)
resultsNames(dds_long)

summary(res_short)
resultsNames(dds_short)


plotMA(res_long)
plotMA(res_short)
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\deseq_v2\\short\\"
setwd(outdir)
report <- DESeq2Report(dds_long, project = "DESeq2 long report", 
                       intgroup = "state", outdir = outdir, output = "index")

report_s <- DESeq2Report(dds_short, project = "DESeq2 short report", 
                       intgroup = "state", outdir = outdir, output = "index")

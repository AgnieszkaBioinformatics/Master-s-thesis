library(tidyr)
library(dplyr)

############# preparing RNA samples for the pipeline
### converting featurecounts to tpm
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\3h1_lhp1\\counts_v2\\counts_gene\\"

fc_preprocess <- function(sample) {
  file_dir <- file.path(paste0(outdir, sample, ".txt"))
  df <- read.table(file_dir, header = TRUE, sep = "\t")
  
  df_long <- df %>%
    separate_rows(Chr, Start, End, Strand, sep = ";") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))
  #name <- paste0("df_", sample)
  #assign(name, df_long)
  
  df_long <- df_long[!duplicated(df_long$Geneid), ]
  return(df_long)
}


samples <- c('3h1_2_sorted_s', '3h1_3_sorted_s', '3h1_4_sorted_s', 
             'Ihp1_1_sorted_s', 'Ihp1_3_sorted_s', 'Ihp1_4_sorted_s',
             'Ihp1_3h1_2_sorted_s', 'Ihp1_3h1_3_sorted_s', 'Ihp1_3h1_4_sorted_s',
             'wt_1_sorted_s', 'wt_3_sorted_s', 'wt_4_sorted_s')

for (s in samples) {
  name <- paste0("X", s)
  assign(name, fc_preprocess(s))
}



calculate_tpm <- function(df) {
  bam_col <- df %>% 
    select(last_col()) %>%
    pull()
  
  df <- df %>%
    mutate(Length_kb = Length / 1000) %>%
    mutate(RPK = bam_col / Length_kb)
  
  scaling_factor <- sum(df$RPK) / 1e6
  df <- df %>% 
    mutate(tpm = RPK / scaling_factor)
  
  return(df)
}

df_samples <- c('X3h1_2_sorted_s', 'X3h1_3_sorted_s', 'X3h1_4_sorted_s', 
             'XIhp1_1_sorted_s', 'XIhp1_3_sorted_s', 'XIhp1_4_sorted_s',
             'XIhp1_3h1_2_sorted_s', 'XIhp1_3h1_3_sorted_s', 'XIhp1_3h1_4_sorted_s',
             'Xwt_1_sorted_s', 'Xwt_3_sorted_s', 'Xwt_4_sorted_s')

for (df in df_samples) {
  current_df <- get(df)
  updated_df <- calculate_tpm(current_df)
  assign(df, updated_df)
}

outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\3h1_lhp1\\counts_v2\\counts_gene\\tpm\\"

for (s in df_samples) {
  write.table(s, paste0(outdir, s), sep="\t")
}

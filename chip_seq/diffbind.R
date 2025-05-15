library(DiffBind)
library(tidyverse)


                      ################ DiffBind ################

s_dir <- "C:\\Users\\aurin\\Downloads\\metadata_f.csv"
samples <- read.csv(s_dir, stringsAsFactors = FALSE)

## Create DBA object
s_dba <- dba(sampleSheet = samples, scoreCol = 5)
plot(s_dba)
#dba.plotHeatmap(s_dba)

## counts
# binding matrix based on read counts
counts_paired <- dba.count(s_dba)

## normalizing
normalized <- dba.normalize(counts_paired)

dba.plotPCA(normalized)
## establishing a model design and a contrast
contrast <- dba.contrast(normalized, categories = DBA_FACTOR, minMembers=2)

# Set WT samples to FALSE
contrast$contrasts[[6]]$group1[grepl("^WT_", names(contrast$contrasts[[6]]$group1))] <- FALSE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[6]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[6]]$group1))] <- TRUE

# Set WT samples to FALSE
contrast$contrasts[[6]]$group2[grepl("^WT_", names(contrast$contrasts[[6]]$group2))] <- TRUE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[6]]$group2[grepl("^lhp1_3h1_", names(contrast$contrasts[[6]]$group2))] <- FALSE

contrast$contrasts[[6]]$name1 <- "lhp1_3h1"
contrast$contrasts[[6]]$name2 <- "WT"

results <- dba.analyze(contrast)
dba.show(results, bContrasts=TRUE)
dba.plotPCA(results)
dba.plotHeatmap(results)
dba.plotHeatmap(results, th=0.01)
#plot(results, contrast=4)

outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\diffbind\\"
setwd(outdir)

report1 <- dba.report(results, th=0.01, fold=1, contrast = 1)
report2 <- dba.report(results, th=0.01, fold=1, contrast = 2)
report3 <- dba.report(results, th=0.01, fold=1, contrast = 3)
report4 <- dba.report(results, th=0.01, fold=1, contrast = 4)
report5 <- dba.report(results, th=0.01, fold=1, contrast = 5)
report6 <- dba.report(results, th=0.01, fold=1, contrast = 6)

## saving 
dba.save(results, file="results_obj", 
         dir = outdir)


## stats
sum(report1$Fold>1) # 185   # które ważne? porównać z profilowaniem rnaseq
sum(report1$Fold<0) # 96

sum(report2$Fold>1) # 6
sum(report2$Fold<0) # 4

sum(report3$Fold>1) #137
sum(report3$Fold<0) # 76

sum(report4$Fold>1) # 56
sum(report4$Fold<0) # 158

sum(report5$Fold>1) # 0
sum(report5$Fold<0) # 0

sum(report6$Fold>1) # 87
sum(report6$Fold<0) # 301

## Venn diagrams
dba.plotVenn(results, contrast = c(3,6), th=0.01, bDB=TRUE)   
dba.plotVenn(results, contrast = 1, th=0.01, bGain=TRUE, 
             bAll=FALSE, bLoss=TRUE, bDB=TRUE)

## MA plot
dba.plotMA(results, th=0.01, fold=1)

## volcano plot
dba.plotVolcano(results, th=0.01, fold=1)

## pvalue boxplot
dba.plotBox(results, th=0.01, fold=1)



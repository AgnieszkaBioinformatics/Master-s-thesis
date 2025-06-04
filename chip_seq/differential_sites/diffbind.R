library(DiffBind)
library(tidyverse)
install.packages("parallelly", type = "binary")
BiocManager::install("DiffBind", force=TRUE)
update.packages(ask = FALSE, checkBuilt = TRUE)

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


# Set WT samples to FALSE
contrast$contrasts[[1]]$group1[grepl("^X3h1_", names(contrast$contrasts[[1]]$group1))] <- FALSE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[1]]$group1[grepl("^lhp1_", names(contrast$contrasts[[1]]$group1))] <- TRUE
contrast$contrasts[[1]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[1]]$group1))] <- FALSE
# Set WT samples to FALSE
contrast$contrasts[[1]]$group2[grepl("^X3h1_", names(contrast$contrasts[[1]]$group2))] <- TRUE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[1]]$group2[grepl("^lhp1_", names(contrast$contrasts[[1]]$group2))] <- FALSE

contrast$contrasts[[1]]$name1 <- "lhp1"
contrast$contrasts[[1]]$name2 <- "3h1"


# Set WT samples to FALSE
contrast$contrasts[[2]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[2]]$group1))] <- TRUE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[2]]$group1[grepl("^X3h1_", names(contrast$contrasts[[2]]$group1))] <- FALSE
# Set WT samples to FALSE
contrast$contrasts[[2]]$group2[grepl("^lhp1_3h1_", names(contrast$contrasts[[2]]$group2))] <- FALSE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[2]]$group2[grepl("^X3h1_", names(contrast$contrasts[[2]]$group2))] <- TRUE

contrast$contrasts[[2]]$name1 <- "lhp1_3h1"
contrast$contrasts[[2]]$name2 <- "3h1"



# Set WT samples to FALSE
contrast$contrasts[[4]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[4]]$group1))] <- TRUE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[4]]$group1[grepl("^lhp1_", names(contrast$contrasts[[4]]$group1))] <- FALSE
contrast$contrasts[[4]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[4]]$group1))] <- TRUE
# Set WT samples to FALSE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[4]]$group2[grepl("^lhp1_", names(contrast$contrasts[[4]]$group2))] <- TRUE
contrast$contrasts[[4]]$group2[grepl("^lhp1_3h1_", names(contrast$contrasts[[4]]$group2))] <- FALSE

contrast$contrasts[[4]]$name1 <- "lhp1_3h1"
contrast$contrasts[[4]]$name2 <- "lhp1"


# Set WT samples to FALSE
contrast$contrasts[[6]]$group1[grepl("^lhp1_3h1_", names(contrast$contrasts[[6]]$group1))] <- TRUE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[6]]$group1[grepl("^WT", names(contrast$contrasts[[4]]$group1))] <- FALSE

# Set WT samples to FALSE

# Set lhp1_3h1 samples to TRUE
contrast$contrasts[[6]]$group2[grepl("^lhp1_3h1_", names(contrast$contrasts[[6]]$group2))] <- FALSE
contrast$contrasts[[6]]$group2[grepl("^WT_", names(contrast$contrasts[[6]]$group2))] <- TRUE

contrast$contrasts[[6]]$name1 <- "lhp1_3h1"
contrast$contrasts[[6]]$name2 <- "WT"

results <- dba.analyze(contrast)
dba.show(results, bContrasts=TRUE)
dba.plotPCA(results)
dba.plotHeatmap(results)
dba.plotHeatmap(results, th=0.05)
#plot(results, contrast=4)

outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\diffbind\\"
setwd(outdir)

report1 <- dba.report(results, th=0.05, fold=1, contrast = 1)
report2 <- dba.report(results, th=0.05, fold=1, contrast = 2)
report3 <- dba.report(results, th=0.05, fold=1, contrast = 3)
report4 <- dba.report(results, th=0.05, fold=1, contrast = 4)
report5 <- dba.report(results, th=0.05, fold=0.5, contrast = 5)
report6 <- dba.report(results, th=0.05, fold=1, contrast = 6)

report_t <- dba.report(results, th=0.01, fold=1, contrast = 1:2, bDB=TRUE)

## saving 
dba.save(results, file="results_obj2", 
         dir = outdir)

results <- dba.load(file="results_obj2")

## stats
sum(report1$Fold>1) # 252
sum(report1$Fold<0) # 120

sum(report2$Fold>1) # 9
sum(report2$Fold<0) # 10

sum(report3$Fold>1) # 168
sum(report3$Fold<0) # 99

sum(report4$Fold>1) # 68
sum(report4$Fold<0) # 199

sum(report5$Fold>0.5) # 8
sum(report5$Fold<0) # 16

sum(report6$Fold>1) # 126
sum(report6$Fold<0) # 400

## Venn diagrams  -> to są nie odfiltrowane wyniki
dba.plotVenn(results, contrast=1, bDB=TRUE,
               bGain=TRUE, bLoss=TRUE, bAll=FALSE)


## MA plot
dba.plotMA(results, th=0.01, fold=1)

## volcano plot
dba.plotVolcano(results, th=0.01, fold=1)

## pvalue boxplot
dba.plotBox(results, th=0.05, fold=1)

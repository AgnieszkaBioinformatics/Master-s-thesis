library(DiffBind)
library(tidyverse)


################ DiffBind

s_dir <- "C:\\Users\\aurin\\Downloads\\metadata_f.csv"
samples <- read.csv(s_dir, stringsAsFactors = FALSE)

# Create DBA object
s_dba <- dba(sampleSheet = samples, scoreCol = 5)
plot(s_dba)
#dba.plotHeatmap(s_dba)

# counts
# binding matrix based on read counts
counts_paired <- dba.count(s_dba)

## normalizing
normalized <- dba.normalize(counts_paired)

dba.plotPCA(normalized)
## establishing a model design and a contast
contrast <- dba.contrast(normalized, categories = DBA_FACTOR, minMembers=2)
results <- dba.analyze(contrast)
dba.show(results, bContrasts=TRUE)
plot(results, contrast=4)

report <- dba.report(results, th=0.01)


## saving 
dba.save(results, file="results_obj", 
         dir = "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\diffbind")


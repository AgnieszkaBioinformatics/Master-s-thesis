library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(clusterProfiler)

diff_dir <- "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\diffbind\\"
setwd(diff_dir)
samplefiles <- list.files(diff_dir, pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)

tables <- list()

for (t in 1:length(samplefiles)) {
  name <- basename(samplefiles[[t]])
  rt <- read.table(samplefiles[[t]])
  rt[,1] <- sub("^chr", "", rt[,1])
  write.table(rt, file = name, 
              sep = "\t",
              row.names = FALSE,  # Removes row numbers
              col.names = FALSE,
              quote = FALSE)
}

txdb <- TxDb.Athaliana.BioMart.plantsmart28

lhp1_WT <- read.table(samplefiles[[7]])
# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")


X3h1_vs_WT <- as.data.frame(peakAnnoList[[1]]@anno)
lhp1_3h1_vs_3h1 <- as.data.frame(peakAnnoList[[2]]@anno)
lhp1_3h1_vs_lhp1 <- as.data.frame(peakAnnoList[[3]]@anno)
lhp1_3h1_vs_WT <- as.data.frame(peakAnnoList[[4]]@anno)
lhp1_vs_3h1 <- as.data.frame(peakAnnoList[[5]]@anno)
lhp1_vs_WT <- as.data.frame(peakAnnoList[[6]]@anno)
###### write annotations
outdir <- "C:\\Users\\aurin\\Desktop\\magisterka\\chipseq\\chipseeker\\"
setwd(outdir)
annotations_l <- list("X3h1_vs_WT" = X3h1_vs_WT, 
                      "lhp1_3h1_vs_3h1" = lhp1_3h1_vs_3h1, 
                      "lhp1_3h1_vs_lhp1" = lhp1_3h1_vs_lhp1,
                      "lhp1_3h1_vs_WT" = lhp1_3h1_vs_WT, 
                      "lhp1_vs_3h1" = lhp1_vs_3h1, 
                      "lhp1_vs_WT" = lhp1_vs_WT)
annotations_s <- c("X3h1_vs_WT", "lhp1_3h1_vs_3h1", "lhp1_3h1_vs_lhp1",
                   "lhp1_3h1_vs_WT", "lhp1_vs_3h1", "lhp1_vs_WT")

for (i in 1:length(annotations_l)) {
  name <- paste0(outdir, annotations_s[[i]], "_annotation.txt")
  write.table(annotations_l[[i]], file = name, sep="\t", 
              quote=FALSE, row.names=FALSE)
}

######## GO enrichment analysis

for (i in 1:length(annotations_l)) {
  name <- paste0(outdir, annotations_s[[i]], "_enrichment.csv")
  id <- annotations_l[[i]]$geneId %>% 
    as.character() %>% 
    unique()
  ego <- enrichGO(gene = id, 
                         keyType = "TAIR", 
                         OrgDb = org.At.tair.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)
  
  # Output results from GO analysis to a table
  summary <- data.frame(egoT)
  write.csv(summary, name)
}

## X3h1_vs_WT      -----------------> nic nie ma
entrezids_X3h1_vs_WT <- X3h1_vs_WT$geneId %>% 
  as.character() %>% 
  unique()

ego_X3h1_vs_WT <- enrichGO(gene = entrezids_X3h1_vs_WT, 
                keyType = "TAIR", 
                OrgDb = org.At.tair.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)



## lhp1_3h1_vs_3h1   ----------------> nic nie ma
entrezids_lhp1_3h1_vs_3h1 <- lhp1_3h1_vs_3h1$geneId %>% 
  as.character() %>% 
  unique()

ego_lhp1_3h1_vs_3h1 <- enrichGO(gene = entrezids_lhp1_3h1_vs_3h1, 
                           keyType = "TAIR", 
                           OrgDb = org.At.tair.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = TRUE)



## lhp1_3h1_vs_lhp1
entrezids_lhp1_3h1_vs_lhp1 <- lhp1_3h1_vs_lhp1$geneId %>% 
  as.character() %>% 
  unique()

ego_lhp1_3h1_vs_lhp1 <- enrichGO(gene = entrezids_lhp1_3h1_vs_lhp1, 
                           keyType = "TAIR", 
                           OrgDb = org.At.tair.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = TRUE)



## lhp1_3h1_vs_WT ----------------> nic nie ma
id_lhp1_3h1_vs_WT <- lhp1_3h1_vs_WT$geneId %>% 
  as.character() %>% 
  unique()

ego_lhp1_3h1_vs_WT <- enrichGO(gene = id_lhp1_3h1_vs_WT, 
                                 keyType = "TAIR", 
                                 OrgDb = org.At.tair.db, 
                                 ont = "BP", 
                                 pAdjustMethod = "BH", 
                                 qvalueCutoff = 0.05, 
                                 readable = TRUE)



## lhp1_vs_3h1
id_lhp1_vs_3h1 <- lhp1_vs_3h1$geneId %>% 
  as.character() %>% 
  unique()

ego_lhp1_vs_3h1 <- enrichGO(gene = id_lhp1_vs_3h1, 
                               keyType = "TAIR", 
                               OrgDb = org.At.tair.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = TRUE)



## lhp1_vs_WT
id_lhp1_vs_WT <- lhp1_vs_WT$geneId %>% 
  as.character() %>% 
  unique()

ego_lhp1_vs_WT <- enrichGO(gene = id_lhp1_vs_WT, 
                            keyType = "TAIR", 
                            OrgDb = org.At.tair.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)


# Dotplot visualization
dotplot(ego_X3h1_vs_WT, showCategory=50)


kegg_lhp1_3h1_vs_WT <- enrichKEGG(gene = id_lhp1_3h1_vs_WT,
                    organism = 'ath',
                    pvalueCutoff = 0.05)

dotplot(kegg_lhp1_3h1_vs_WT)  # ----------> też nic



### all samples comparison
# Create a list with genes from each sample
genes = lapply(annotations_l, function(i) as.data.frame(i)$geneId)

comp <- compareCluster(geneCluster = genes, 
                           fun = "enrichGO",
                       keyType = "TAIR",
                       OrgDb = org.At.tair.db,
                       ont = "BP",
                       qvalueCutoff = 0.05, 
                           pAdjustMethod = "BH")
dotplot(comp, showCategory = 20, title = "GO Pathway Enrichment Analysis")



## integration with lhp1 rna seq profile
rna_dir <- "C:\\Users\\aurin\\Desktop\\magisterka\\3h1_lhp1\\dge\\v3\\"
sig_lhp1 <- read.csv(paste0(rna_dir, "sig_de_Ihp1.csv"))
background <- read.csv(paste0(rna_dir, "res_df_Ihp1.csv"))
annotations_l <- list("X3h1_vs_WT" = X3h1_vs_WT, 
                      "lhp1_3h1_vs_3h1" = lhp1_3h1_vs_3h1, 
                      "lhp1_3h1_vs_lhp1" = lhp1_3h1_vs_lhp1,
                      "lhp1_3h1_vs_WT" = lhp1_3h1_vs_WT, 
                      "lhp1_vs_3h1" = lhp1_vs_3h1, 
                      "lhp1_vs_WT" = lhp1_vs_WT,
                      "lhp1_RNA-seq" = sig_lhp1)

genes = lapply(annotations_l, function(i) as.data.frame(i)$geneId)
genes$`lhp1_RNA-seq` <- sig_lhp1$X

comp <- compareCluster(geneCluster = genes, 
                       fun = "enrichGO",
                       keyType = "TAIR",
                       OrgDb = org.At.tair.db,
                       ont = "BP",
                       qvalueCutoff = 0.05, 
                       pAdjustMethod = "BH")



dotplot(comp, showCategory = 50, font.size = 5,
        title = "GO Pathway Enrichment Analysis \n lhp1 RNA-seq profile")


comp_rna <- compareCluster(geneCluster = list("lhp1_rna" = sig_lhp1$X), 
                           fun = "enrichGO",
                           keyType = "TAIR",
                           universe = background$X,
                           OrgDb = org.At.tair.db,
                           ont = "BP",
                           qvalueCutoff = 0.05, 
                           pAdjustMethod = "BH")

dotplot(comp_rna, showCategory = 50, title = "lhp1 RNA-seq profile")
emapplot(comp, showCategory = 50)

### intgracja rnaseq -> chipseq niezbyt informatywna
### chyba lepiej chipseq -> rnaseq
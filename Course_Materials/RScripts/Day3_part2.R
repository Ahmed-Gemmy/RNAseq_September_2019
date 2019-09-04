# Annotation

library(biomaRt)
library(DESeq2)
library(tidyverse)

load("Robjects/DE.RData")

listMarts()
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl) %>%
  filter(str_detect(description, "Mouse"))

ensembl  = useDataset("mmusculus_gene_ensembl", mart=ensembl)


#query

## we need to give it the following:
### 1. filter (gene Ids)
### 2. values (where our filters are)
### 3. attributes (what we want back)

listFilters(ensembl) %>%
  filter(str_detect(name, "ensembl"))

ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(resLvV)[1:1000]

listAttributes(ensembl) %>%
  head(20)
attributeNames <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name")

annot <- getBM(attributes = attributeNames, 
               filters = ourFilterType, 
               values = filterValues,
               mart=ensembl)

annot %>%
  add_count(ensembl_gene_id) %>%
  filter(n>1)

## Challenge
### 1. attributes: gene description biotype
### 2. all genes not just 1000
### 3. how many have multiple Entrez IDs?
### 4. how many don't have annotation?



listAttributes(ensembl) %>% 
  filter(str_detect(name, "description"))

listAttributes(ensembl) %>% 
  filter(str_detect(name, "biotype"))
attributeNames <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "description", "gene_biotype")
filterValues <- rownames(resLvV)

annot <- getBM(attributes = attributeNames, 
               filters = ourFilterType, 
               values = filterValues,
               mart=ensembl)



multiEntrez <-annot %>%
  filter(entrezgene_id != 'NA') %>%
  add_count(entrezgene_id) %>%
  filter(n>1)

dim(multiEntrez)

noAnnotation <- annot %>%
  filter(entrezgene_id == 'NA' & external_gene_name=='NA' & description == 'NA' & gene_biotype=='NA')

# In Class solution
load("Robjects/Full_annotation.RData")

dup<-annot %>%
  add_count(ensembl_gene_id) %>%
  filter(n>1) %>%
  distinct(ensembl_gene_id) %>%
  nrow()
### other way to do it
sum(duplicated(annot$ensembl_gene_id))

## To find missing genes
missingGenes <- !rownames(resLvV)%in%annot$ensembl_gene_id
rownames(resLvV)[missingGenes]



#### Now we need to combine annotation table and result table #####

load("Robjects/Ensembl_annotations.RData")

annotLvV <- as.data.frame(resLvV) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

write_tsv(annotLvV, "results/VirginVsLactating_Result_Annontated.txt")

annotLvV %>%
  arrange(FDR) %>%
  head(10)



### Visualization

ddsShrink <- lfcShrink(ddsObj, coef = "Status_lactate_vs_virgin")

shrinkLvL <- as.data.frame(ddsShrink) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

## pvalue hist

hist(shrinkLvL$pvalue)
#hist(shrinkLvL$FDR)


## MA plots

plotMA(ddsShrink, alpha=0.05)

cutoff <- sort(shrinkLvL$pvalue) [10]
shrinkLvL <- shrinkLvL %>%
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, Symbol, ""))

ggplot(shrinkLvL, 
       aes(x= log2(baseMean), y=logFC)) +
  geom_point(aes(colour=FDR<0.05), shape=20, size=0.5) +
  geom_text(aes(label=TopGeneLabel)) +
  labs(x="mean of normalised counts", y="log fold change")


## Challenge Volcano Plot

shrinkLvL <- shrinkLvL %>%
  mutate(log10=-log10(FDR))

ggplot(shrinkLvL, 
       aes(x= logFC, y=log10)) +
  geom_point(aes(colour=FDR<0.05), shape=20, size=0.5) +
  #geom_text(aes(label=TopGeneLabel)) +
  labs(x="LogFC", y="-Log10(FDR)")

# In Class solution

# Volcano Plot

shrinkLvV <- shrinkLvL %>%
  mutate(`-log10pvalue`= -log10(pvalue))

ggplot(shrinkLvV, 
       aes(x= logFC, y=`-log10pvalue`)) +
  geom_point(aes(colour=FDR<0.05), shape=20, size=0.5) +
  labs(x="LogFC", y="-log10pvalue")

## strip chart

topgene <- filter(shrinkLvV, Symbol=="Wap")
geneID <- topgene$GeneID

plotCounts(ddsObj, gene = geneID, intgroup = c("CellType", "Status"),
           returnData = T) %>% 
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Status), shape=21, size=2) +
  facet_wrap(~CellType) +
  expand_limits(y=0)


# Heatmaps
library(ComplexHeatmap)
library(circlize)

# getting signeficant genes

sigGenes <- as.data.frame(shrinkLvV) %>% 
  top_n(150, wt=-FDR) %>% 
  pull("GeneID")


plotDat <- vst(ddsObj)[sigGenes,] %>% 
  assay()  # just pulls out the counts
z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE)) 

# colour palette
myPalette <- c("red3", "ivory", "blue3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_names = FALSE,
        cluster_columns = FALSE)

# Better Heatmap

hcDat <- hclust(dist(z.mat))
cutGroups <- cutree(hcDat, h=4)

ha1 <- HeatmapAnnotation(df = colData(ddsObj)[, c("CellType", "Status")  ])
Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_names = FALSE,
        cluster_columns = FALSE,
        split = cutGroups,
        rect_gp = gpar(col = "darkgrey", lwd=0.5),
        top_annotation = ha1)



save(annotLvV, shrinkLvV, file="results/Annotated_Results_LvV.RData")












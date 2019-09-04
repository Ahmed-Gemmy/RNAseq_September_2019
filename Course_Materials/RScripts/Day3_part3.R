# Functional Analysis

## GSEA
library(tidyverse)
library(fgsea)


load("Robjects/Annotated_Results_LvV.RData")

### ranking
## -- a way of ranking our genes --

gseaDat <- filter(shrinkLvV, !is.na(Entrez))

rankData <- gseaDat$logFC
#then we need to add geneids to this vector

names(rankData) <- gseaDat$Entrez
head(rankData)


# Load Pathways

load("Robjects/mouse_H_v5.RData")
Mm.H[[1]]
names(Mm.H)[1]

pathwaysH <- Mm.H

## Conduct analysis

fgseaRes <- fgsea(pathwaysH, 
                  rankData,
                  minSize = 15, 
                  maxSize =  500,
                  nperm=1000)

head(fgseaRes)

### pathway        pval        padj         ES           NES            nMoreExtreme                                                      size                               leadingEdge
### pathway name   pval        FDR    enrichment score  Normalised(ES)          in the 1000 randomly rnaking our genes how extreme was it

### pvalue depends on the value of nperm .... 
  
plotEnrichment(pathwaysH[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], rankData)

# GO enrichment

## goseq
library(goseq)

# list of significant genes

sigData <- as.integer( shrinkLvV$FDR <0.01 & !is.na(shrinkLvV$FDR) )

names(sigData) <- shrinkLvV$GeneID

# fit probability weight function
pwf <- nullp(sigData, "mm10", "ensGene", bias.data = shrinkLvV$medianTxLength)

# run the GO analysis

goResults <- goseq(pwf, "mm10", "ensGene", test.cats = c("GO:BP"))

head(goResults)

# plot top 10 go terms

goResults %>%
  top_n(10, wt=-over_represented_pvalue) %>%
  mutate(histPerc=numDEInCat*100/numDEInCat) %>%
  ggplot(aes(x=histPerc, y=term, colour=over_represented_pvalue, size=numDEInCat)) +
  geom_point() + 
  expand_limits(x=0) +
  labs(x="hist %", y="Go Term", colour="p value", "size"="Count")


# KEGG pathway
library(clusterProfiler)

# significant gene list
sigGenes <- shrinkLvV %>% 
  filter(FDR < 0.05 & !is.na(FDR) & 
           abs(logFC) > 1 & 
           !is.na(Entrez)) %>% 
  pull(Entrez)

## run KEGG analysis

keggResults <- enrichKEGG(gene = sigGenes, organism = 'mmu')
head(keggResults, n=10)

browseKEGG(keggResults, 'mmu03320')

# Generate a figure

library(pathview)
























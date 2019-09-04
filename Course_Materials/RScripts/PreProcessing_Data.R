library(DESeq2)
library(tidyverse) # used to manipulate the TSV data

# read sample info
sampleinfo <- read_tsv("data/SampleInfo.txt")
sampleinfo

# read count matrix
seqdata <- read_tsv("data/GSE60450_Lactation.featureCounts", comment = "#")
seqdata

#dplyr

#make a copy of sampleinfo
newTable <-sampleinfo
newTable

#filter rows
newTable <- filter(newTable, CellType == "basal")
newTable

# select columns
newTable <- select(newTable, CellType, Group)
newTable

#rename CellType
newTable <- rename(newTable, Cell=CellType)
newTable

# using pipe (%>%)
newTable <- sampleinfo %>%
  filter(CellType=="basal") %>%
  select(CellType, Group) %>%
  rename(Cell=CellType)

newTable

# Just another way of doing it by '[]' filter
# newTable <- sampleinfo
# newTable <- newTable[newTable$CellType=="basal", ]
# newTable

# reformat count matrix
countdata <- seqdata %>%
  column_to_rownames("Geneid") %>%
  rename_all(str_remove,".bam") %>%
  select(sampleinfo$Sample) %>%
  as.matrix()

dim(countdata)
head(countdata)

# filter genes
keep <- rowSums(countdata) > 5
table(keep)

countdata <- countdata[keep,]
head(countdata)
dim(countdata)


# QA# library size


librarySizes <- colSums(countdata)

barplot(librarySizes,
        names=names(librarySizes),
        las=2,
        main="Barplot of library sizes"
        )
abline(h=20e6, lty=2)

# distribution pf read counts
logcounts <- log2(countdata +1)

# color
statusCol <- match(sampleinfo$Status, c("virgin", "pregnant", "lactate")) + 1
statusCol

#boxplot
boxplot(logcounts,
        xlab="",
        ylab="log2(Counts)",
        col=statusCol)
abline(h=median(as.matrix(logcounts)), col="blue")

# Challenge 1
# User rlog() to transform counts
# plot

rlogcount <- rlog(countdata)

# color
statusCol <- match(sampleinfo$Status, c("virgin", "pregnant", "lactate")) + 1
#statusCol <- as.numeric(factor(sampleinfo$Sample)) +1
#statusCol

#boxplot
boxplot(rlogcount,
        xlab="",
        ylab="rLog(Counts)",
        las=2,
        col=statusCol)
abline(h=median(as.matrix(rlogcount)), col="blue")


# PCA
library(ggfortify)

#run PCA

pcDat <- prcomp( t(rlogcount) )
class(pcDat)

# plot PCA

autoplot(pcDat)

# coloring points
autoplot(pcDat,
         data = sampleinfo,
         colour = "CellType",
         shape = "Status",
         size = 5)

# label samples with name
library(ggrepel)
autoplot(pcDat,
         data = sampleinfo, 
         colour = "CellType",
         shape="Status",
         size = 5) + 
  geom_text_repel( aes(x=PC1, y=PC2, label=Sample), box.padding = 0.8 )

# correct sample info

sampleinfo_copy <- sampleinfo
sampleinfo <- sampleinfo %>%
  mutate( CellType = ifelse(Sample == "MCL1.DG", "basal", CellType)) %>%
  mutate(CellType = ifelse(Sample == "MCL1.LA", "luminal", CellType)) %>%
  mutate(Group=str_c(CellType, ".", Status))
sampleinfo

autoplot(pcDat,
         data = sampleinfo, 
         colour = "CellType",
         shape="Status",
         size = 5) + 
  geom_text_repel( aes(x=PC1, y=PC2, label=Sample), box.padding = 0.8 )


# make DESeq object
sampleinfo <- read_tsv("data/SampleInfo.Corrected.txt")
sampleinfo

all( sampleinfo$Sample == colnames(countdata))


design <- as.formula( ~ CellType )

ddsObj <- DESeqDataSetFromMatrix(countData = countdata, 
                                 colData = sampleinfo,
                                 design = design)
class(ddsObj)

# Normalisation
ddsObj <- estimateSizeFactors(ddsObj)
sizeFactors(ddsObj)

# MA plot
library(limma)

# log counts non-normatilsed
par( mfrow=c(1,2) )
plotMA(logcounts, array=7)
abline(h=0, col="red")
plotMA(logcounts, array=11)
abline(h=0, col="red")

# normalized counts
normalisedCounts <- counts(ddsObj, normalized=TRUE)
class(normalisedCounts)
dim(normalisedCounts)
#log the normalised counts
logNormalisedCounts <-  log2(normalisedCounts+1)
par( mfrow=c(2,2) )
plotMA(logNormalisedCounts, array=7)
abline(h=0, col="red")
plotMA(logNormalisedCounts, array=11)

abline(h=0, col="red")
plotMA(logNormalisedCounts, array=10)
abline(h=0, col="red")
plotMA(logNormalisedCounts, array=9)
abline(h=0, col="red")

save(countdata, sampleinfo, file="results/preprocessing.RData")
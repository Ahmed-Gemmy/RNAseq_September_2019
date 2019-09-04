library(tidyverse)
library(DESeq2)

load("Robjects/preprocessing.RData")
dim(countdata)
sampleinfo

# Create a DESeq Object
## count data
## sample info
## Model (which are the designs for our experiments)

design <- as.formula(~CellType + Status ) ## ~telda means function of something.. what are the parameteres in the experiment that changing gene expression
#modelMatrix <- model.matrix(design, data=sampleinfo)
#modelMatrix

## It makes more sense to have 'virgin' in the intercept .. R takes them ordered alphbetically .. hence we will change this below

sampleinfo
sampleinfo$Status <- factor(sampleinfo$Status, levels = c("virgin", "pregnant", "lactate"))

sampleinfo$Status

modelMatrix <- model.matrix(design, data=sampleinfo)
modelMatrix


# Create DESeq2 object

ddsObj.row <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleinfo, design = design)

ddsObj.row


## PCA
vstcount <- vst(ddsObj.row, blind = TRUE)
plotPCA(vstcount, intgroup=c("Status", "CellType"))

# The DESeq2 workflow

####STEPS####

## 1. Estimate size factors

ddsObj <- estimateSizeFactors(ddsObj.row)
sizeFactors(ddsObj)

## 2. Estimate dispersion paramter

ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

## 3. Apply the negative bionomial model

ddsObj <- nbinomWaldTest(ddsObj)

#DESeq command
ddsObj <- DESeq(ddsObj.row) #One command doning all the analysis 

# Generate results table

res <- results(ddsObj, alpha = 0.05)
head(res)


resLvV <- res

# Getting other contrasts

##Main effect
resultsNames(ddsObj) #tells us what main effects are there .. default is always last one which we got in privious run ... LvV

resPvV <- results(ddsObj, name="Status_pregnant_vs_virgin", alpha = 0.05)
head(resPvV)

# Retrieving 'top genes' by pvalue
topGen <- as.data.frame(resPvV) %>%
  rownames_to_column("GeneID") %>%
  top_n(-padj, n=100)
topGen

# Chalenge
# Obtain result for contrast luminal v basal and find top 200 genes by padj (FDR)
resultsNames(ddsObj)

resLvB <- results(ddsObj, name = "CellType_luminal_vs_basal", alpha = 0.05 )
topGen2 <- as.data.frame(resPvV) %>%
  rownames_to_column("GeneID") %>%
  top_n(-padj, n=200)
topGen2

## obtaining unnamed contrast

resLvP <- results(ddsObj,
                  contrast = c("Status", "lactate", "pregnant"),  #we have to specify 1) column name, 2)the contrasts
                  alpha = 0.05)
head(resLvP)

# Compare tow models

designC <- as.formula(~CellType) #I think only cell type is important

# test the design

ddsObjC <- DESeq(ddsObj, test="LRT", reduced=designC) #Compare to this simpler design

resCvCS <- results(ddsObjC)
head(resCvCS)

### Now we got 2 models how to figure out which one is better?! #####
### You have to think about Biology and other factors you are studying. #####
### One way is look at FDR and see if the simpler model is segnifically different than the more complex model .. if it is not we can drop Status ###


table(resCvCS$padj <= 0.05)
## FALSE  TRUE 
## 13010  7932
### All trues indicate that complex model is more segnificant than the simpler one... in this case we have 13K doesn't matter.


### Challenge 2

# 1. fit model with interaction (replace '+' with '*')
# 2. Use the "LRT" test to compare models
# 3. Is the interaction more appropriate?
# Do we have enough replication for the interactive model?

designI <- as.formula(~CellType*Status)
ddsObj.rowI <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleinfo, design = designI)
ddsObjI <- DESeq(ddsObj.rowI)

resultsNames(ddsObjI)

ddsObjCM <- DESeq(ddsObjI, test="LRT", reduced = design)
resIvNI <- results(ddsObjCM)
table(resIvNI$padj <= 0.05)

##### In class solution

## create new DESeq object
designI <- as.formula(~CellType*Status)
### OR
designI <- as.formula(~CellType + Status + CellType:Status)

ddsObjI <-DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = designI)
ddsObjI <- DESeq(ddsObjI)

## Test the interactive model vs addiditive model

ddsObjIvA <- DESeq(ddsObjI, test = "LRT", reduced = design)
resIvA <- results(ddsObjIvA)
head(resIvA)

table(resIvA$padj <=0.05)
### for 9k genes it is the better model and biologically it makes sense


## pulling conrasts from the interactive model
resultsNames(ddsObjI)

### Contrast Luminal vs Basil in virgin mice
## our base line (intercept) is the virgin mice, so this is just the main effect.
resLvBV <- results(ddsObjI, name = "CellType_luminal_vs_basal", alpha = 0.05)

head(resLvBV)


### Contrast Luminal v Basal in pregnant mice
resLvBP <- results(ddsObjI, 
                   contrast = list(c("CellType_luminal_vs_basal", "CellTypeluminal.Statuspregnant")),
                   # Basically adding the interactive term to the main effect.
                   alpha = 0.05
                   )
head(resLvBP)



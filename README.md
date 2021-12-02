---
title: "Final_Project"
author: "Maricarmen Pachicano"
date: "11/30/2021"
output: html_document
---


# Setting up R environment and Dataset input:
## HTSeq-Count Input
I have imported a total of 80 HTSeq-count files into a directory, and is displayed in table sampleTable. There are n = 40 from each treatment group in my data file (20 Female, 20 Male for mulitple comparisons).

## Creating a variable called 'directory' which points to where my output files are located:
Please change path according to where you saved the folder 'files'. 
```{r}
directory <- "/Users/maricarmen/Desktop/files"
```
## Specifying files that contain "treated" in the filename, and using that to identify the condition status:
```{r}
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleType <- sub(".*treated","",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          type = sampleType)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$type <- factor(sampleTable$type)
```

## Building the DESeqDataSet using the ddsHTSeq function:
```{r}
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq
```
## Pre-filtering low count genes before running DESeq2 functions:
```{r}
keep <- rowSums(counts(ddsHTSeq)) >= 10
dds <- ddsHTSeq[keep,]
```

## Explicitly setting factor levels for comparison of treated vs. untreated:
```{r}
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
```

# Differential Expression Analysis: 
The differential expression analysis is done with ONE single function, the DESeq function. Results table is given with function 'res', and it's output includes log2 fold changes, p values and adjusted p values. The first few lines of the results table explains comparison for log2 fold changes between treated/untreated samples. 
```{r}
dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, name="condition_treated_vs_untreated")
res
```

## Log fold change shrinkage: 
This is a necessary step in order to properly visualize and analyze various genes. This is done with the function lfcShrink. Here, the apeglm method is used for shrinkage. We first set dds as the object, and then set the name of the coefficient we want to shrink.
```{r}
resultsNames(dds)
```
```{r}
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```

## Ordering results table (by p-values): 
We can order results by the smallest p-value and store in resOrdered.
```{r}
resOrdered <- res[order(res$pvalue),]
```
```{r}
summary(res)
```

## How many p-values are less than 0.1? There are approximately 201 p-values that are less than 0.1.
```{r}
sum(res$padj < 0.1, na.rm=TRUE)
```

## Adjusting p-value to 0.5
```{r}
res05 <- results(dds, alpha=0.05)
summary(res05)
```

## How many p-values are less than 0.5?
```{r}
sum(res05$padj < 0.05, na.rm=TRUE)
```

# Exploring and exporting results
## MAPlot for shrunkin log2 fold changes:
This MAplot shows the log2 fold changes for each variable over the mean of normalized counts for all samples in the dataset. Red dots indicate p values less than 0.1, and open triangles either pointing up or down indicate points that fall out of window. 
```{r}
plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
```
## Alternative shrinkage estimators: 
Instead of setting coefficient by its name, we can alternitevely set the coefficient by the order it appears in resultsNames(dds) by setting coef = 2
```{r}
resultsNames(dds)
```

## Setting 'coef = 2' with the 'ashr' method (this may take awhile, grab a coffee):
```{r}
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```
## Plotting MAplots for apeglm, normal, and ashr methods.
```{r}
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```
## Plot counts: 
By using the plotCounts function, we can analyze counts of reads for one gene across both treated and untreated groups. plotCounts uses normalized counts and adds a psuedocount of 1/2 for log scale plotting.The gene ploted below, ENSG0000175426.9, had the smallest p-value from the results table. 
```{r}
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```
## Plotting using ggplot: using the function returnData to call for a ggplot.
```{r}
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```
# Multi-factor design
Our sample population has a condition of interest (treated for liver cancer or untreated), but we also have the Sex of each patient under "type". By using the colData function, we see that we have to chop/clean up our 'type' names. This is not the most elegant way to grab 'type' data, but it's the only way I could figure out how to add a multiple factor design and it works!
## With colData, you can see condition and type as the headers. Because of the long names of each samples, we need to do some clean up. The following steps cleans up the "Type" so we can use it in our multi-factor analysis. 
```{r}
colData(dds)
```

First, I create a copy of the DDS so we can run a multi-factor analysis.
```{r}
ddsMF <- dds
```

Then, we change how "Type" is displayed for each sample participant.
```{r}
levels(ddsMF$type)
```

We will chop/clean up "Type" so that it only specifies F for female patient or M for male patient. 
```{r}
levels(ddsMF$type) <- sub("_.*", "", levels(ddsMF$type))
levels(ddsMF$type)
```

As condition is the variable of interest, we put it at the end of our formula. We then re-run DeSeq analysis:
```{r}
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
```

We access results using the results function:
```{r}
resMF <- results(ddsMF)
head(resMF)
```

Now we use the contrast function to extract log2 fold changes while specifying Type
```{r}
resMFType <- results(ddsMF,
                     contrast=c("type", "F", "M"))
head(resMFType)
```

## Info on results columns: 
More information regarding which tests were used can be used via the mcols function.
```{r}
mcols(res)$description
```

## Exporting to CSV file: 
Results of differentially expressed genes can be exported using the base R function: write.csv. The output will be located on the Desktop or your relative path.  Although it does not output the HUGO ID, it does provide the Ensembl.I will address this as a Known Issue. 
```{r}
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
```

Exporting results with ONLY  genes that pass an adjusted p-value threshold can be done using the subset function, followed by the write.csv function.
```{r}
resSig <- subset(resOrdered, padj < 0.1)
resSig
```

```{r}
write.csv(as.data.frame(resSig), 
          file="condition_treated_results.csv")
```


# Data transformations and visualizations
## Extracting transformed values: 
The expression of variance for the variance stabilizing transformation is used by the 'vst' function. On the other hand, the regularized log transformation is used by the 'rlog' function (note: this takes a VERY long time to run, take a break go get a coffee).  
```{r}
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 1)
```

## Effects on transformation on variance for the normal logarithm transformation (ntd):
SdPlot is used to plot the expression of different transformations. The SdPlots plot the standard deviation transformation data across all samples, against the mean, using the shifted log transformation, the regularized log transformation, and the variance stabilizing transformation. 
```{r}
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```
## Effects on transformation on variance (Standard Deviation) for the variance stabilizing transformation (vsd):
```{r}
meanSdPlot(assay(vsd))
```
## Effects on transformation on variance (Standard Deviation) for the regularized log transformation (rld):
```{r}
meanSdPlot(assay(rld))
```

# Heatmap of count matrix
Known issue: I had to remove colnames because they were too complicated, and the heatmap wasn't being displayed. Will need to come back to this. 
## This is a heatmap of the multi-factor design (type + condition):
The heatmaps below are of the various transformed data. 
```{r}
library("pheatmap")
select <- order(rowMeans(counts(ddsMF,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsMF)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
## Heatmap using ntd data described above:
```{r}
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsMF)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
## Heatmap using vsd data described above:
```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

## Heatmap using rld data described above:
```{r}
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
# Sample clustering: 
Another way transformation data can be utilized is by sample clustering. 
```{r}
sampleDists <- dist(t(assay(vsd)))
```

## Creating a heatmap of the sample-to-sample distances:
This heatmap shows  the similarities and dissimiliarities between each sample. 
```{r}
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
# Principle Component plots
## PCA plotting:
This PCA plot shows the samples in a 2D plane spanned  by their first two princial components. These plots are utilized for visualizing and analyzing the overall effect of sample variance. 
This is a PCA plot using vsd transformed data
```{r}
plotPCA(vsd, intgroup="condition", ntop = 500, returnData = FALSE)
```

# Variations to the standard workflow
## Performing Wald test on dataset: 
The Wald test is a hypothesis test where the estimated standard error of a log2 fold change is used in order to see if it is equal to 0.  
```{r}
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
```

## Performing Likelihood ratio test on dataset: 
The likelihood ratio test is another type of hypothesis test, but it examines two models for the counts: a full model and reduced model. This hypothesis test determines if the increased liklihood of the data using the extra terms in the full model is more than expected if those extra terms are 0. 
```{r}
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
```

## False sign or small rate in MA-plot: 
This uses s-values as opposed to p-values (Stephens 2016). The plot is of log fold change for s-values against normalized counts. 
```{r}
resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
```
# Dispersion plot and fitting alternatives
## Dispersion plot:
plotDispEsts is useful for plotting dispersion estimates. As you can more or less see, the final estimates shrunk from the gene-wise estimates towards the fitted estimates. But it is difficult to see since there are so many samples. Typically, with less samples, you are able to see the amount of shrinkage, the number of coefficients, and the row mean and variability of gene-wise estimates. 
```{r}
plotDispEsts(dds)
```
## Custom dispersion fit: 
Using the dispersionFunction we can replace dispersion estimates above a threshold, and  using the estimateDispersionGeneEst function we can provide fitted values for dispersion estimation. 
```{r}
ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                     na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)
```

## Independent filtering of results: 
The DESeq vignette performs independent filtering using the mean of normalized counts and a threshold of significance level alpha. 
```{r}
metadata(res)$alpha
```
```{r}
metadata(res)$filterThreshold
```

## Supplying a custom dispersion fit:
```{r}
ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                     na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)
```
## Independent filtering of results:
```{r}
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
```
## Turn off independent filtering:
```{r}
resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))
```

## 4 possible alternative hypotheses in MA-plot:
The four possible values are greater than absolute value (for two-tail tests), less than absolute value (for p-values that are the maximum of the upper and lower tests), greater than and less than. Depending on what you want to test and analyze, you will chose accordingly.Make sure to turn off independent filtering and reload/rerun Wald test before running. Note: takes long time to run, go grab a coffee. 
```{r}
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
```
## Access to all calculated values: 
All row-wise calculated values (intermediate dispersion calculation, coefficients, standard errors, etc) are stored in the DESeqDataSet object (dds). You can access these values via the mcols function. 
```{r}
mcols(dds,use.names=TRUE)[1:4,1:4]
```
## The function substr is only used for display purposes here: 
```{r}
substr(names(mcols(dds)),1,10) 
```
```{r}
mcols(mcols(dds), use.names=TRUE)[1:4,]
```

## Mean values and Cook's distances: 
Mean values and Cook's distances for each gene and sample are stored  as matrices in asssays slot. 
```{r}
head(assays(dds)[["mu"]])
```

```{r}
head(assays(dds)[["cooks"]])
```

# Dispersion
##Dispersions can be accessed with the 'dispersion' function. 
```{r}
head(dispersions(dds))
```
```{r}
head(mcols(dds)$dispersion)
```

# Size factors: 
##Size factors are available with the sizeFactors function:
```{r}
sizeFactors(dds)
```

# extracting the matrix [βir] for all genes i and model coefficients r: the coef function can be used for extracting the matrix for all genes and model coefficients. 
```{r}
head(coef(dds))
```

## B prior variance: 
The beta prior variance is also stored in the vignette. 
```{r}
attr(dds, "betaPriorVar")
```

# General Information (packages, info on log fold change shrinkage)
## General information using the 'prior' function:
```{r}
priorInfo(resLFC)
```
```{r}
priorInfo(resNorm)
```

```{r}
priorInfo(resAsh)
```

## Dispersion prior variance: 
Dispersion prior variance is also stored in the vignette. 
```{r}
dispersionFunction(dds)
```
```{r}
attr(dispersionFunction(dds), "dispPriorVar")
```

## The version of DESeq2 is stored here. 
```{r}
metadata(dds)[["version"]]
```

# Sample/Gene-dependent Nomralization factors: 
Here the vignette provides a matrix with row-wise geometric means of 1. In order to do this, we divide out the currernt row geometric means. 
```{r}
normFactors <- matrix(runif(nrow(dds)*ncol(dds),0.5,1.5),
                      ncol=ncol(dds),nrow=nrow(dds),
                      dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors
```

# Count outlier detection: 
Plotting the maximun value of Cook's distance for each row over the rank of the test statistic. 
```{r}
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
```
# Filtering criteria: 
This can be done by graphing the mean of normalized counts regardless of condition
```{r}
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
```
# Histogram p-values for all tests: 
The area shaded in blue indicates the portion of those who pass the filtering, and those in yellow who do not pass the filtering. 
```{r}
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
```
```{r}
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
```
# Section 2: Data
The data file is a zipped file located in the shared Google Drive named XXXX. Please unzip and add it to your Desktop. If you save the file elsewhere, please change the path in the "Input Data" portion of the vignette. 

# Section 3: Known issues
1) When I download my "cart" from https://portal.gdc.cancer.gov/, I have to unzip each individual file myself. I think using linux there is a faster way, but I am not sure how to do this yet. For now, I have unzipped each HTSeq file one by one and placed them in one folder. 
2) I had to re-name the file names so that the Condition and Type for each sample was labeled. This is not the most elegant way to do this, but it works. 
3) The file names are very long because each sample has a unique ID. Because of this, I needed to cut Sample ID's out of my heatmaps because it was not loading the graph successfully. I will need to figure out a easy way to rename all the sample ID's. 
4) My CSV file only includes Ensembl ID's. I was not able to figure out code to replace Ensembl with Hugo. I will need to come back to this. 

---
title: "Final_Project"
author: "Maricarmen Pachicano"
date: "12/3/2021"
output: html_document
---

# Differential gene expression in treated and untreated patients with Early Stage (I/II/III) Hepatocellular carcinoma Liver Cancer
Some of the biggest lifestyle risk factors for Liver Cancer is alcohol abuse, smoking and obesity. Although this disease highly impacts substantial populations in the United States, majority of the time it goes untreated. I will be looking at differences in gene expression in early-stage Hepatocellular carcinoma liver cancer patients who are either treated or remain untreated. This analysis will utilize the package DESeq2 vignette. 


# Section 1: DESeq2 Vignette 

## HTSeq-Count Input
I have imported a total of **80 HTSeq-count files** into a dataframe. There are n = 40 from each treatment group in my data file (20 Female, 20 Male for mulitple comparisons)

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
![](https://user-images.githubusercontent.com/89544326/144116365-ba6415f1-63e8-4082-b749-503100b86618.png)

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
![](https://user-images.githubusercontent.com/89544326/144496138-dd7179eb-25ad-41ed-a222-40975449c874.png)

## Log fold change shrinkage: 
This is a necessary step in order to properly visualize and analyze various genes. This is done with the function lfcShrink. Here, the apeglm method is used for shrinkage. We first set dds as the object, and then set the name of the coefficient we want to shrink.
```{r}
resultsNames(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496145-3715f368-e812-4899-8b2b-5df065b39a09.png)

```{r}
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```
![](https://user-images.githubusercontent.com/89544326/144496147-60e38c44-89ae-4ea1-bee8-c74309a00e9e.png)

## Ordering results table (by p-values): 
We can order results by the smallest p-value and store in resOrdered.
```{r}
resOrdered <- res[order(res$pvalue),]
```
```{r}
summary(res)
```
![](https://user-images.githubusercontent.com/89544326/144496152-f23a60bc-2b45-42f9-a56f-438504da79ee.png)

## How many p-values are less than 0.1? There are approximately 201 p-values that are less than 0.1.
```{r}
sum(res$padj < 0.1, na.rm=TRUE)
```
![](https://user-images.githubusercontent.com/89544326/144496153-aa024758-c157-49bc-ad6b-476c001cd486.png)

## Adjusting p-value to 0.5
```{r}
res05 <- results(dds, alpha=0.05)
summary(res05)
```
![](https://user-images.githubusercontent.com/89544326/144496155-596b4282-e459-413f-9649-a98f8ab40207.png)

## How many p-values are less than 0.5?
```{r}
sum(res05$padj < 0.05, na.rm=TRUE)
```
![](https://user-images.githubusercontent.com/89544326/144496156-ec4b1fbe-87e0-4c6b-bd72-83bf38c8bb53.png)

# Exploring and exporting results
## MAPlot for shrunkin log2 fold changes:
This MAplot shows the log2 fold changes for each variable over the mean of normalized counts for all samples in the dataset. Red dots indicate p values less than 0.1, and open triangles either pointing up or down indicate points that fall out of window. 
```{r}
plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
```
![](https://user-images.githubusercontent.com/89544326/144496159-1a2645da-6e58-4594-a7a8-9fa8ef9830d3.png)

## Alternative shrinkage estimators: 
Instead of setting coefficient by its name, we can alternitevely set the coefficient by the order it appears in resultsNames(dds) by setting coef = 2
```{r}
resultsNames(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496163-b337e0fd-df1e-459b-882e-be72ad50e4ef.png)

## Setting 'coef = 2' with the 'ashr' method (this may take awhile, grab a coffee):
```{r}
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```
![](https://user-images.githubusercontent.com/89544326/144496164-ef5f4d24-e694-49b6-ba5c-ff24dc31ee64.png)

## Plotting MAplots for apeglm, normal, and ashr methods.
```{r}
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```
![](https://user-images.githubusercontent.com/89544326/144496165-7e42d228-984b-4785-99f1-5926ee32be2b.png)

## Plot counts: 
By using the plotCounts function, we can analyze counts of reads for one gene across both treated and untreated groups. plotCounts uses normalized counts and adds a psuedocount of 1/2 for log scale plotting.The gene ploted below, ENSG0000175426.9, had the smallest p-value from the results table. 
```{r}
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```
![](https://user-images.githubusercontent.com/89544326/144496167-3581e095-817c-4aea-b19f-5591b6ebe097.png)

## Plotting using ggplot: using the function returnData to call for a ggplot.
```{r}
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```
![](https://user-images.githubusercontent.com/89544326/144496169-616add01-c63f-4965-a632-79d3e548b18d.png)

# Multi-factor design
Our sample population has a condition of interest (treated for liver cancer or untreated), but we also have the Sex of each patient under "type". By using the colData function, we see that we have to chop/clean up our 'type' names. This is not the most elegant way to grab 'type' data, but it's the only way I could figure out how to add a multiple factor design and it works!
## With colData, you can see condition and type as the headers. Because of the long names of each samples, we need to do some clean up. The following steps cleans up the "Type" so we can use it in our multi-factor analysis. 
```{r}
colData(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496172-3eac4158-77f5-45d9-ba6c-4f160886c16d.png)

First, I create a copy of the DDS so we can run a multi-factor analysis.
```{r}
ddsMF <- dds
```

Then, we change how "Type" is displayed for each sample participant.
```{r}
levels(ddsMF$type)
```
![](https://user-images.githubusercontent.com/89544326/144496174-182752f7-f439-48ee-96a5-494f1e9f2c16.png)

We will chop/clean up "Type" so that it only specifies F for female patient or M for male patient. 
```{r}
levels(ddsMF$type) <- sub("_.*", "", levels(ddsMF$type))
levels(ddsMF$type)
```
![](https://user-images.githubusercontent.com/89544326/144496176-b7f93863-e0bf-49f2-98fe-a9f508840a9b.png)

As condition is the variable of interest, we put it at the end of our formula. We then re-run DeSeq analysis:
```{r}
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
```
![](https://user-images.githubusercontent.com/89544326/144496178-34fd743c-e782-410b-9789-05ff232d7e37.png)

We access results using the results function:
```{r}
resMF <- results(ddsMF)
head(resMF)
```
![](https://user-images.githubusercontent.com/89544326/144496180-55cd6a9b-e18e-4a98-a17d-77363305bf24.png)

Now we use the contrast function to extract log2 fold changes while specifying Type
```{r}
resMFType <- results(ddsMF,
                     contrast=c("type", "F", "M"))
head(resMFType)
```
![](https://user-images.githubusercontent.com/89544326/144496185-8bd85942-44ea-403f-8ac9-de96b1cde97a.png)

## Info on results columns: 
More information regarding which tests were used can be used via the mcols function.
```{r}
mcols(res)$description
```
![](https://user-images.githubusercontent.com/89544326/144496186-7439c5b3-d852-41e9-91f9-445047e4c1c1.png)

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
![](https://user-images.githubusercontent.com/89544326/144496188-db1d4d6f-fc1d-4b9f-9357-becb17d9d1b5.png)

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
![](https://user-images.githubusercontent.com/89544326/144496190-d1441ca6-667d-410c-bb21-5e8eeadb02e8.png)

## Effects on transformation on variance for the normal logarithm transformation (ntd):
SdPlot is used to plot the expression of different transformations. The SdPlots plot the standard deviation transformation data across all samples, against the mean, using the shifted log transformation, the regularized log transformation, and the variance stabilizing transformation. 
```{r}
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```
![](https://user-images.githubusercontent.com/89544326/144496192-17b18709-fc6d-459a-affe-89c496693880.png)

## Effects on transformation on variance (Standard Deviation) for the variance stabilizing transformation (vsd):
```{r}
meanSdPlot(assay(vsd))
```
![](https://user-images.githubusercontent.com/89544326/144496195-aeecf209-bb30-44e6-a718-e9a8fcdc425d.png)


## Effects on transformation on variance (Standard Deviation) for the regularized log transformation (rld):
```{r}
meanSdPlot(assay(rld))
```
![](https://user-images.githubusercontent.com/89544326/144496197-5ab4049b-aab2-4ce0-a3a6-49168dc5ae8b.png)

# Heatmap of count matrix 
Overall, this analysis (shown in the following heatmaps below) illustrate differential gene expression between groups. As you can see in the figures, it appears as though the "Treated" group has an increase differential gene expression compared to the "Untreated" groups. In addition, there is also less differential genetic expression in the Male groups. 

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
![](https://user-images.githubusercontent.com/89544326/144496198-a05014e6-780d-48c9-892b-d63cc021e682.png)

## Heatmap using ntd data described above:
```{r}
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsMF)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

![](https://user-images.githubusercontent.com/89544326/144496199-24ef2435-a736-44ff-84e5-f3bf392a90cc.png)

## Heatmap using vsd data described above:
```{r}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![](https://user-images.githubusercontent.com/89544326/144496200-ffbdf234-8aec-4d22-bbc4-91e336fcc32b.png)

## Heatmap using rld data described above:
```{r}
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![](https://user-images.githubusercontent.com/89544326/144496203-dd42fe3a-594e-4260-9d5d-bc45770e0ded.png)

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
![](https://user-images.githubusercontent.com/89544326/144496204-f3f77827-e6c8-42b9-bf72-d2479c53e7bd.png)

# Principle Component plots
## PCA plotting:
This PCA plot shows the samples in a 2D plane spanned  by their first two princial components. These plots are utilized for visualizing and analyzing the overall effect of sample variance. Each dot on the PCA plot below represent one of the 80 samples included in this analysis. As you can see, the Treated and Untreated groups are 'mixed' along the projected plane. This indicates that there is not much difference between the two groups. 
This is a PCA plot using vsd transformed data
```{r}
plotPCA(vsd, intgroup="condition", ntop = 500, returnData = FALSE)
```
![](https://user-images.githubusercontent.com/89544326/144496205-69bc29ec-2146-42a0-bad0-a6bf4df0968a.png)

However, if you add in Sex for multiple comparisons, you can see that there is a difference between male and female, treated and untreated patients. This is indicated by less mixing, and more grouping of same groups, specifically the Untreated group (shown in red). 

```{r}
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```
![](https://user-images.githubusercontent.com/89544326/144516408-ba16ed61-cd2a-4c6c-b5bf-c2e17f49dc43.png)


# Variations to the standard workflow
## Performing Wald test on dataset: 
The Wald test is a hypothesis test where the estimated standard error of a log2 fold change is used in order to see if it is equal to 0.  
```{r}
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496206-681dd97d-e492-408f-82cd-41c3bd7fde0c.png)

## Performing Likelihood ratio test on dataset: 
The likelihood ratio test is another type of hypothesis test, but it examines two models for the counts: a full model and reduced model. This hypothesis test determines if the increased liklihood of the data using the extra terms in the full model is more than expected if those extra terms are 0. 
```{r}
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496208-547e0271-27d2-40aa-b71e-68d0e7d05975.png)

## False sign or small rate in MA-plot: 
This uses s-values as opposed to p-values (Stephens 2016). The plot is of log fold change for s-values against normalized counts. 
```{r}
resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
```
![](https://user-images.githubusercontent.com/89544326/144496209-a94dcf19-16be-45d6-a439-958e126655ab.png)

# Dispersion plot and fitting alternatives
## Dispersion plot:
plotDispEsts is useful for plotting dispersion estimates. As you can more or less see, the final estimates shrunk from the gene-wise estimates towards the fitted estimates. But it is difficult to see since there are so many samples. Typically, with less samples, you are able to see the amount of shrinkage, the number of coefficients, and the row mean and variability of gene-wise estimates. 
```{r}
plotDispEsts(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496213-0f3459dc-20f8-4d88-b9a3-a8f4471e37b6.png)

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

![](https://user-images.githubusercontent.com/89544326/144496215-95ef1982-d83a-4740-b8b5-eae62ab45177.png)

## Independent filtering of results: 
The DESeq vignette performs independent filtering using the mean of normalized counts and a threshold of significance level alpha. 
```{r}
metadata(res)$alpha
```
![](https://user-images.githubusercontent.com/89544326/144496217-7a1c3631-906a-44c5-858f-b4a53836500c.png)

```{r}
metadata(res)$filterThreshold
```
![](https://user-images.githubusercontent.com/89544326/144496218-2ec64c58-b0f8-4b18-8fe0-dfec3124b877.png)

## Supplying a custom dispersion fit:
```{r}
ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],
                     na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)
```
![](https://user-images.githubusercontent.com/89544326/144496219-244f5ee8-745e-4ed2-a5e8-ef113a4cab8d.png)

## Independent filtering of results:
```{r}
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
```
![](https://user-images.githubusercontent.com/89544326/144496222-6358351b-e29a-4df2-8b89-4a6eed8dfc46.png)

## Turn off independent filtering:
```{r}
resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))
```
![](https://user-images.githubusercontent.com/89544326/144496223-01283e7b-658a-486c-b731-9400d45e1e18.png)

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
![](https://user-images.githubusercontent.com/89544326/144496224-f0fdfd9b-a5cf-436c-a245-d18f40cb2f9b.png)

## Access to all calculated values: 
All row-wise calculated values (intermediate dispersion calculation, coefficients, standard errors, etc) are stored in the DESeqDataSet object (dds). You can access these values via the mcols function. 
```{r}
mcols(dds,use.names=TRUE)[1:4,1:4]
```
![](https://user-images.githubusercontent.com/89544326/144496226-80f2256f-06b5-408d-91d3-4ac3e1ea02ea.png)

## The function substr is only used for display purposes here: 
```{r}
substr(names(mcols(dds)),1,10) 
```
![](https://user-images.githubusercontent.com/89544326/144496227-c262db89-e772-42ec-ae08-e601d0cde306.png)

```{r}
mcols(mcols(dds), use.names=TRUE)[1:4,]
```
![](https://user-images.githubusercontent.com/89544326/144496228-752d0bd0-ba73-41eb-8c87-bd153bc335d5.png)

## Mean values and Cook's distances: 
Mean values and Cook's distances for each gene and sample are stored  as matrices in asssays slot. 
```{r}
head(assays(dds)[["mu"]])
```
![](https://user-images.githubusercontent.com/89544326/144496229-3f96257d-27c6-4b82-bb26-679daed99c37.png)

```{r}
head(assays(dds)[["cooks"]])
```
![](https://user-images.githubusercontent.com/89544326/144496230-632e0070-bf86-4ffd-9667-034d877bd465.png)

# Dispersion
## Dispersions can be accessed with the 'dispersion' function. 
```{r}
head(dispersions(dds))
```
![](https://user-images.githubusercontent.com/89544326/144496231-9b192c7e-21d4-408c-bc90-ff9ddf960674.png)

# Size factors: 
##Size factors are available with the sizeFactors function:
```{r}
sizeFactors(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496233-53e9fc81-32ea-4765-935f-6d4b7163525b.png)

# Extracting the matrix [βir] for all genes i and model coefficients r: the coef function can be used for extracting the matrix for all genes and model coefficients. 
```{r}
head(coef(dds))
```
![](https://user-images.githubusercontent.com/89544326/144496236-8ca32b69-e9c2-4f51-83d8-006f6efff0d4.png)

## B prior variance: 
The beta prior variance is also stored in the vignette. 
```{r}
attr(dds, "betaPriorVar")
```
![](https://user-images.githubusercontent.com/89544326/144496238-4ff975d5-0d4d-4a99-b7ad-2977e102b307.png)

# General Information (packages, info on log fold change shrinkage)
## General information using the 'prior' function:
```{r}
priorInfo(resLFC)
```
![](https://user-images.githubusercontent.com/89544326/144496242-696fb877-da6f-4c73-8831-5a2819ea670e.png)

```{r}
priorInfo(resNorm)
```
![](https://user-images.githubusercontent.com/89544326/144496243-6859dcac-17aa-43d5-b025-037f76e8bbf5.png)

```{r}
priorInfo(resAsh)
```
![](https://user-images.githubusercontent.com/89544326/144496244-2203431e-e286-4f47-9373-08c1eebb1f10.png)


## Dispersion prior variance: 
Dispersion prior variance is also stored in the vignette. 
```{r}
dispersionFunction(dds)
```
![](https://user-images.githubusercontent.com/89544326/144496245-8e4a9275-faae-467d-9f2b-8e2e6ad198d3.png)

```{r}
attr(dispersionFunction(dds), "dispPriorVar")
```
![](https://user-images.githubusercontent.com/89544326/144496246-03cf37a6-c0c5-47ae-9cae-c1b6a802fdfd.png)

## The version of DESeq2 is stored here. 
```{r}
metadata(dds)[["version"]]
```
![](https://user-images.githubusercontent.com/89544326/144496248-55156017-7cc9-4038-b2b3-e16043d5858f.png)

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
![](https://user-images.githubusercontent.com/89544326/144496249-b86ffe19-2411-42df-a0dc-7490cdd3c9bd.png)

# Filtering criteria: 
This can be done by graphing the mean of normalized counts regardless of condition
```{r}
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
```
![](https://user-images.githubusercontent.com/89544326/144496252-eef8a4ce-b4ac-40ca-8971-647fc60c699e.png)

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
![](https://user-images.githubusercontent.com/89544326/144496254-67407dda-8575-401b-a611-70b6913fc857.png)

# Section 2: Conclusions


# Section 3: Data
The data file is a zipped file located in the shared Google Drive named XXXX. Please unzip and add it to your Desktop. If you save the file elsewhere, please change the path in the "Input Data" portion of the vignette. 

# Section 4: Known issues
1) When I download my "cart" from https://portal.gdc.cancer.gov/, I have to unzip each individual file myself. I think using linux there is a faster way, but I am not sure how to do this yet. For now, I have unzipped each HTSeq file one by one and placed them in one folder. 
2) I had to re-name the file names so that the Condition and Type for each sample was labeled. This is not the most elegant way to do this, but it works. 
3) The file names are very long because each sample has a unique ID. Because of this, I needed to cut Sample ID's out of my heatmaps because it was not loading the graph successfully. I will need to figure out a easy way to rename all the sample ID's. 
4) My CSV file only includes Ensembl ID's. I was not able to figure out code to replace Ensembl with Hugo. I will need to come back to this. 

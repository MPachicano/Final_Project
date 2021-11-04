# Final Project Outline

## Title
*Differential gene expression in Treated and Untreated patients with Early Stage (I/II/III) Hepatocellular carcinoma Liver Cancer*

## Author
Maricarmen Pachicano

## Overview of Project
Some of the biggest lifestyle risk factors for Liver Cancer is alcohol abuse, smoking and obesity. Although this disease highly impacts substantial populations in the United States, majority of the time it goes untreated. I will be looking at differences in gene expression in early-stage Hepatocellular carcinoma liver cancer patients who are either treated or remain untreated. This analysis will utilize the package DESeq2 and follow this vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 

My analysis will use the TCGA cohort, and I have identified about 404 cases HT-Seq Counts files, with 289 liver cancer patients who received no treatment, and 33 patients who did receive treatment. I have confirmed that these patients have HT-Seq Counts files that are able to be downloaded, and that is what the DESeq vignette requires. 

## Data
I will use data from https://portal.gdc.cancer.gov/repository. Out of the 322 Liver samples in TCGA, I have identified a total of 289 who did not receive treatment and 33 patients who did receive treatment. Patients have downloadable HT-Seq Counts files, and I have ensured that I am able to download and open them via the R application. 

## Milestone 1
I plan on fully loading the data into vignette through HT-Seq steps, and attempt first round of analyses.

## Milestone 2
I will complete the entire first round of analysis in the vignette, and will compose written analyses alongside the vignette.

## Deliverable 
I will deliver a complete repository of differential gene expression in Treated and Untreated Liver patients along with a description of my analysis and results.

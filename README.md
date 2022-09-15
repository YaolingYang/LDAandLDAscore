# Linkage Disequilibrium of Ancestry (LDA)

This is the location for the LDA and LLDA score tools that were used in Barrie et al. 2022 (in prep).

## Introduction
Linkage Disequilibrium of Ancestry (LDA) quantifies the correlations between the ancestry of two SNPs, measuring the proportion of individuals who have experienced a recombination leading to a change in ancestry, relative to the genome-wide baseline.

LDA score is the total amount of genome in LDA with each SNP (measured in recombination map distance), which is useful for detecting the signal of selection.

A detailed description can be found in "**LDA and LDA score.pdf**".

This work is done by Yaoling Yang (<yaoling.yang@bristol.ac.uk>) under the supervision of Dr. Daniel Lawson (<dan.lawson@bristol.ac.uk>).

## Download R package "LDAandLDAS"
R Package "LDAandLDAS" includes two functions "LDA" and "LDAS" for calculating LDA and LDA score, respectively.

You can download the R package by typing the following codes in R:
```
devtools::install_github("YaolingYang/LDAandLDAscore")
```

## A simple example:
```
library(LDAandLDAS)

#combine the painting data for two ancestries
data=list(painting_p1[,-1],painting_p2[,-1])

#calculate the pairwise LDA of SNPs
LDA_result <- LDA(data,SNPlimit=1200)

#map is the data containing the physical position and recombination distance of the SNPs
#calculate the LDA score for the SNPs
LDA_score <- LDAS(LDA_result,map,window=10)

#visualise the LDA score (pd is the physical distance of a SNP)
plot(x=LDA_score$pd,y=LDA_score$LDAS)
```

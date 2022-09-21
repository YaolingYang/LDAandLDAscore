# Linkage Disequilibrium of Ancestry (LDA)

This is the location for the LDA and LDA score tools that were used in Barrie et al. 2022 (submitted).

* Author: Yaoling Yang (yaoling.yang@bristol.ac.uk)   

* Supervisor: Daniel Lawson (dan.lawson@bristol.ac.uk)  

* License: GPL-3.0-or-later

## Introduction
Linkage Disequilibrium of Ancestry (LDA) quantifies the correlations between the ancestry of two SNPs, measuring the proportion of individuals who have experienced a recombination leading to a change in ancestry, relative to the genome-wide baseline.

LDA score is the total amount of genome in LDA with each SNP (measured in recombination map distance), which is useful for detecting the signal of selection.

A detailed description is available here: [LDA and LDA score.pdf](https://github.com/YaolingYang/LDAandLDAscore/blob/master/LDA%20and%20LDA%20score.pdf).

## Install R package "LDAandLDAS"
R Package "LDAandLDAS" includes two functions "LDA" and "LDAS" for calculating LDA and LDA score, respectively.

You can download the R package by typing the following codes in R:
```
devtools::install_github("YaolingYang/LDAandLDAscore")
```

## Examples
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

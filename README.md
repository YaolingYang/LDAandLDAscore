# Linkage Disequilibrium of Ancestry (LDA) and LDA score

This is the location for the LDA and LDA score tools that were used in [Barrie W, Yang
Y, Attfield K E, et al. Elevated genetic risk for multiple sclerosis originated in Steppe Pastoralist populations. Nature
(2023)](https://www.biorxiv.org/content/10.1101/2022.09.23.509097v2).

-   Authors:  
    Yaoling Yang (<yaoling.yang@bristol.ac.uk>)  
    Daniel Lawson (<dan.lawson@bristol.ac.uk>)

-   License: GPL-3

## Introduction
Linkage Disequilibrium of Ancestry (LDA) quantifies the correlations between the ancestry of two SNPs, measuring the proportion of individuals who have experienced a recombination leading to a change in ancestry, relative to the genome-wide baseline.

LDA score is the total amount of genome in LDA with each SNP (measured in recombination map distance), which is useful for detecting the signal of selection.

## Install R package "LDAandLDAS" from CRAN
```
install.packages("LDAandLDAS")
```

## Examples
```
library(LDAandLDAS)

# visualize the painting data
# Painting data are the average probabilities of different populations
head(LDAandLDAS::example_painting_p1[1:5,],10)

# combine the painting data for two ancestries as a list
# to make to input data for function 'LDA'.
paintings=list(LDAandLDAS::example_painting_p1,
          LDAandLDAS::example_painting_p2)
          
# calculate the pairwise LDA of SNPs
LDA_result <- LDA(paintings)

# map is the data containing two columns
# The first column is the physical position (unit: b) (decreasing order)
# The second column is the recombination distance (unit: cM) of the SNPs
head(LDAandLDAS::example_map,10)

# calculate the LDA score for the SNPs
LDA_score <- LDAS(LDA_result,LDAandLDAS::example_map,window=10)

#visualize the LDA scores
plot(x=LDA_score$SNP,y=LDA_score$LDAS)
```

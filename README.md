Linkage Disequilibrium of Ancestry (LDA) quantifies the correlations between the ancestry of two SNPs, measuring the proportion of individuals who have experienced a recombination leading to a change in ancestry, relative to the genome-wide baseline.

LDA score is the total amount of genome in LDA with each SNP (measured in recombination map distance), which is useful for detecting the signal of selection.

A detailed description can be found in LDA_and_LDA_score.pdf.

R Package "LDAandLDAS" includes two functions "LDA" and "LDAS" for calculating LDA and LDA score, respectively.

You can download the R package by typing the following codes in R:
```
devtools::install_github("YaolingYang/LDAandLDAscore/LDAandLDAS")
```

A simple example:
```
library(LDAandLDAS)

#painting data for two ancestries

data1 <- read.csv("https://raw.githubusercontent.com/YaolingYang/LDAandLDAscore/master/painting_p1.csv")

data2 <- read.csv("https://raw.githubusercontent.com/YaolingYang/LDAandLDAscore/master/painting_p2.csv")

data=cbind(data1[,-1],data2[,-1])

#calculate the pairwise LDA of SNPs

LDA_result <- LDA(data,ancestry=2)

#A file with physical position and recombination distance of the SNPs

map <- read.csv("https://raw.githubusercontent.com/YaolingYang/LDAandLDAscore/master/map.csv")

#calculate the LDA score for the SNPs

LDA_score <- LDAS(LDA_result,map,window=4)
```

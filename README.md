
# galgo

<!-- badges: start -->
<!-- badges: end -->

The goal of galgo is  the Identification and Study of Prognostic Gene Expression Signatures in Cancer

## Installation

You can install the released version of galgo using devtools with:

``` r
devtools::install_github("https://github.com/harpomaxx/galgo")
```

## Example

This is a basic example which shows how to find gene Expression signatures for the LUAD Dataset

``` r

library(galgo)

## basic example code
load("/home/galgo/galgo/Data/RNA_LUAD_all.rda")
downl=names(esets)

trainSet=esets[["TCGA"]]
prob_matrix= exprs(trainSet)
clinical=pData(trainSet)
OS=survival::Surv(time=clinical$time,event=clinical$status)
chrom_length= nrow(prob_matrix)   #length of chromosome

galgo::search_ges(generations = 10, population = 30,prob_matrix = prob_matrix, chrom_length = nrow(prob_matrix),OS=OS)

```


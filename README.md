
galgoR <img src="vignettes/GalgoR.png" align="right" alt="" width="120" />
================================================================================
<!-- badges: start -->
<!-- badges: end -->

A multi-objective optimization algorithm for disease subtype discovery based on a  non-dominated sorting genetic algorithm. The galgo framework combines the advantages of clustering algorithms for grouping heterogeneous omics data and the searching properties of genetic algorithms for feature selection and optimal number of clusters determination to find features that maximize the survival difference between subtypes while keeping cluster consistency high.

Installation
-------------
You can install the released version of galgo using devtools with:

``` r
devtools::install_github("https://github.com/harpomaxx/galgo")
```

Example
-------
This is a basic example which shows how to find gene Expression signatures for the LUAD Dataset

``` r

library(galgoR)

## basic example code
rna_luad<-use_rna_luad()

prm <- rna_luad$prm 
clinical <- rna_luad$clinical
OS <- survival::Surv(time=clinical$time,event=clinical$status)

galgoR::galgo(generations = 10, population = 30,prob_matrix = prm, OS=OS)

```
## GPU support

By default galgo runs some portions of its code in GPU, provided by the gpuR package. Before installing gpuR, the opencl backend should be configured. 

In linux systems install lastest nvidia cuda drivers and the opencl backend.

```
       apt-get install nvidia-418 nvidia-opencl-icd-418 libcuda1-418
       apt-get install opencl-headers  ocl-icd-opencl-dev
       
```

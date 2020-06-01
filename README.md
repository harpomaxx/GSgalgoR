
galgoR <img src="inst/extdata/GalgoR.png" align="right" alt="" width="120" />
================================================================================
<!-- badges: start -->
<!-- badges: end -->

A multi-objective optimization algorithm for disease subtype discovery based on a  non-dominated sorting genetic algorithm. The galgo framework combines the advantages of clustering algorithms for grouping heterogeneous omics data and the searching properties of genetic algorithms for feature selection and optimal number of clusters determination to find features that maximize the survival difference between subtypes while keeping cluster consistency high.

Installation
-------------
You can install the released version of galgoR using devtools with:

```
devtools::install_github("https://github.com/harpomaxx/galgo")

library(galgoR)
```

Examples datasets
------------------

To standardize the structure of genomic data, we use the [ExpressionSet](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html) structure for the examples given in this guide. The `ExpressionSet` objects can hold different types of data in a single structure but in this case we opted for using a simplified format to facilitate the example to those not familiar with the [Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html) package. The `ExpressionSet` objects are formed mainly by:

- A matrix of genetic expression, usually derived from microarray or RNAseq experiments. 
- Phenotypic data, where we find information on the samples (condition, status, treatment, survival, and other covariates). 
- Finally, these objects can also contain Annotations and feature Meta-data.

To start testing galgoR, the package contains two reduced lung adenocarcinoma gene expression datasets ([TCGA](https://gdc.cancer.gov/about-data/publications/luad_2014) and [GSE68465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465)), that can be download using the function use_rna_luad(). Additionally, It also contains the [Wilkerson's centroids](http://cancer.unc.edu/nhayes/publications/adenocarcinoma.2012/wilkerson.2012.LAD.predictor.centroids.csv.zip) to perform lung adenocarcinoma sample classification. 
 
```
 rna_luad<- use_rna_luad()
 TCGA<- rna_luad$TCGA #Access TCGA dataset
 GSE68465<- rna_luad$GSE68465 #Access GSE68465 dataset

 #To access gene expression data
 TCGA_expr<- TCGA$expression_data

 #To access feature data
 TCGA_features<- TCGA$feature_data

 #To access clinical data
 TCGA_clinic <- TCGA$pheno_data

 #To get wilkerson centroids
 WilkCentroids <- rna_luad$WilkCentroids
```
# An Example

## Loading data 


```
rna_luad <- use_rna_luad()
```

We will also load the `survival` library

```
library(survival)
```
## Run galgo()

The main function in this package is galgo(). It accepts an expression matrix and survival object to find robust gene expression signatures related to a given outcome.
This function contains some parameters that can be modified, according to the characteristics of the analysis to be performed.

### Setting parameters

The principal parameters are:

- population: a number indicating the number of solutions in the population of solutions that will be evolved
- generations: a number indicating the number of iterations of the galgo algorithm
- nCV: number of cross-validation sets
- usegpu: logical default to FALSE, set to TRUE if you wish to use gpu computing (gpuR package must be properly installed and loaded)
- distancetype: character, it can be 'pearson' (centered pearson), 'uncentered' (uncentered pearson), 'spearman' or 'euclidean'
- TournamentSize: a number indicating the size of the tournaments for the selection procedure
- period: a number indicating the outcome period to evaluate the RMST

```
population <- 30 # For testing reasons it is set to a low number but ideally should be above 100            
generations <-15 # For testing reasons it is set to a low number but ideally should be above 150            
nCV <- 5                      
distancetype <- "pearson"     
TournamentSize <- 2
period <- 1825
```

### Expression matrix

Create an expression matrix for the TCGA_LUAD data example.

```{r TCGA_expr}
TCGA_expr <- rna_luad$TCGA$expression_matrix
```

### Survival Object

The 'OS' object is created by the Surv() function of the survival package. This uses phenotypic data that are contained in the TCGA dataset.  

```{r Surv}
TCGA_clinic <- rna_luad$TCGA$pheno_data

OS <- survival::Surv(time=TCGA_clinic$time,event=TCGA_clinic$status)
```

### Run Galgo algorithm

```
output <- galgoR::galgo(generations = generations, 
                        population = population, 
                        prob_matrix = TCGA_expr, 
                        OS = OS,
                        nCV = nCV, 
                        distancetype = distancetype,
                        TournamentSize = TournamentSize, 
                        period = period)

```


### Galgo Object

The output of the galgo() function is an object of type 'galgo.Obj' that has two slots with the elements:

- Solutions 
- ParetoFront.

#### Solutions 

Is a l x (n + 5) matrix where n is the number of features evaluated and l is the number of solutions obtained. 

- The submatrix l x n is a binary matrix where each row represents the chromosome of an evolved solution from the solution population, where each feature can be present (1) or absent (0) in the solution. 
- Column n+1 represent the k number of clusters for each solutions 
- Column n+2 shows the SC Fitness 
- Column n+3 represent Survival Fitness values
- Column n+4 shows the solution rank
- Column n+5 represent the crowding distance of the solution in the final pareto front

#### ParetoFront

Is a list of length equal to the number of generations run in the algorithm. Each element is a l x 2 matrix where l is the number of solutions obtained and the columns are the SC Fitness and the Survival Fitness values respectively.


For easier interpretation of the 'galgo.Obj', the output can be transformed to a List or to a DataFrame objects.

## toList() function

This function restructurates a galgo.Obj to a more easy to understand an use list. This output is particularly useful if one wants to select a given solution and use its outputs in a new classifier. The output of type list has a length equals to the number of solutions obtained by the galgo algorithm.

Basically this output is a list of lists, where each element of the output is named after the solution's name (solution.n, where n is the number assigned to that solution), and inside of it, it has all the constituents for that given solution with the following structure:

- solution.n$Genes: A vector of the features included in the solution
- solution.n$k: The number of partitions found in that solution
- solution.n$SC.Fit: The average silhouette coefficient of the partitions found
- solution.n$Surv.Fit: The survival fitnes value
- solution.n$Rank: The solution rank
- CrowD: The solution crowding distance related to the rest of the solutions

```
outputList <- toList(output)
head(names(outputList))
```

To evaluate the structure of the first solution we can run:

```
outputList[["Solution.1"]]
```

## toDataFrame() function

The current function restructurates a galgo.Obj to a more easy to understand an use data.frame. The output data.frame has m x n dimensions, were the rownames (m) are the solutions obtained by the galgo algorithm. The columns has the following structure:

- Genes: The features included in each solution in form of a list
- k: The number of partitions found in that solution
- SC.Fit: The average silhouette coefficient of the partitions found
- Surv.Fit: The survival fitnes value
- Rank: The solution rank
- CrowD: The solution crowding distance related to the rest of the solutions

```
outputDF <- toDataFrame(output)
head(outputDF)
```

## plot_pareto()

Once we obtain the `galgo.obj` from the output of `galgo()` we can plot the obtained Pareto front and see how it evolved trough the tested number of generations

```
plot_pareto(output)
```
![](./inst/extdata/images/pareto2.jpg)


GPU support
---
Galgo is able to fasten its execution using GPU computing utilizing the [gpuR package](https://cran.r-project.org/package=gpuR "gpuR R package") . Before installing [gpuR](https://cran.r-project.org/package=gpuR "gpuR R package"), the opencl backend should be configured. 

In linux systems install latest nvidia cuda drivers and the opencl backend.

```
       apt-get install nvidia-418 nvidia-opencl-icd-418 libcuda1-418
       apt-get install opencl-headers  ocl-icd-opencl-dev
       
```

For installing [gpuR](https://github.com/cdeterman/gpuR/wiki) and enable GPU computing in different operating systems, follow [gpuR installation guide](https://github.com/cdeterman/gpuR/wiki "gpuR installation guide").

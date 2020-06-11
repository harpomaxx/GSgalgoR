#' Download a toy dataset to test galgoR
#'
#' The current function downloads two reduced lung adenocarcinoma gene expression datasets (TCGA lung adenocarcinoma project and GEO GSE68465 datasets) to test \code{\link[galgoR:galgo]{galgoR}}. It also contains the Wilkerson's centroids to perform lung adenocarcinoma sample classification.
#' @param userdir the name of the folder to install the datasets. A temporary file will be used if left blank.
#' @return The dataset gets downloaded in a temporary folder or the \code{userdir} folder specified by the user. The output of the function is a \code{list} with two \code{lists} with the scaled gene expression and clinical data from TCGA lung adenocarcinoma and \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465}{GEO GSE68465} datasets.It also contains the Wilkerson's centroids extracted from \url{http://cancer.unc.edu/nhayes/publications/adenocarcinoma.2012/wilkerson.2012.LAD.predictor.centroids.csv.zip} to perform sample classification into "Bronchoid", "Magnoid" and "Squamoid" subtypes.
#' @export
#' @references
#' \itemize{
#'   \item Shedden, K., Taylor, J., Enkemann, S. et al. Gene expression–based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. Nat Med 14, 822-827 (2008). \url{https://doi.org/10.1038/nm.1790}
#'   \item Collisson, E., Campbell, J., Brooks, A. et al. Comprehensive molecular profiling of lung adenocarcinoma. Nature 511, 543-550 (2014). \url{https://doi.org/10.1038/nature13385}
#'   \item Matthew D. Wilkerson, Xiaoying Yin, Katherine A. Hoadley, Yufeng Liu, Michele C. Hayward, Christopher R. Cabanski, Kenneth Muldrew, C. Ryan Miller, Scott H. Randell, Mark A. Socinski, Alden M. Parsons, William K. Funkhouser, Carrie B. Lee, Patrick J. Roberts, Leigh Thorne, Philip S. Bernard, Charles M. Perou and D. Neil Hayes. Clin Cancer Res October 1 2010 (16) (19) 4864-4875 \url{https://doi.org/10.1158/1078-0432.CCR-10-0199}
#'   \item Wilkerson MD, Yin X, Walter V, Zhao N, Cabanski CR, Hayward MC, et al. (2012) Differential Pathogenesis of Lung Adenocarcinoma Subtypes Involving Sequence Mutations, Copy Number, Chromosomal Instability, and Methylation. PLoS ONE 7(5): e36530. \url{https://doi.org/10.1371/journal.pone.0036530}
#' }
#' @examples
#' \dontrun{
#' rna_luad <- use_rna_luad()
#' TCGA <- rna_luad$TCGA # Access TCGA dataset
#' GSE68465 <- rna_luad$GSE68465 # Access GSE68465 dataset
#'
#' # To access gene expression data
#' TCGA_expr <- TCGA$expression_data
#'
#' # To access feature data
#' TCGA_features <- TCGA$feature_data
#'
#' # To access clinical data
#' TCGA_clinic <- TCGA$pheno_data
#'
#' # To get wilkerson centroids
#' WilkCentroids <- rna_luad$WilkCentroids
#' }
#' @importFrom tools md5sum

use_rna_luad <- function(userdir = "") {
  rna_luad_digest <- "900a74e7c4fdd0dcb7a3f2ddb44bb680"
  rna_luad_url <- "https://bit.ly/luad_data_galgo"
  if (userdir == "") {
    dest_dir <- paste0(tempdir(), "/")
  } else {
    dest_dir <- userdir
  }
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }

  dest_file <- paste(dest_dir, "luad_data.rds", sep = "")
  message("Trying to load dataset from ", dest_file)
  download_fail <- 0
  if (!file.exists(dest_file)) {
    message(dest_file, " not found.")
    tryCatch(utils::download.file(rna_luad_url,
      dest_file,
      mode = "wb"
    ),
    warning = function(e) {
      message("rna_luad dataset is not available at [", rna_luad_url, "]. Are you connected to the internet? ")
      download_fail <<- sub(".+HTTP status was ", "", e)
    }
    )
  }
  # message(download_fail)
  if (download_fail == 0 && tools::md5sum(dest_file) == rna_luad_digest) {
    data <- readRDS(dest_file)
  }
  else {
    message("Something goes wrong. Could not down(load) dataset. Please try again later or download it by your own from ", rna_luad_url)
    if (file.exists(dest_file)) {
      file.remove(dest_file)
    }
  }
}


#' Title
#'
#' @param eset 
#' @param map 
#'
#' @return
#' @noRd 
#'
#' @examples
drop_duplicates <- function(eset, map = "Gene.symbol") {

  # Drop NA's
  drop <- which(is.na(Biobase::fData(eset)[, map]))
  eset <- eset[-drop, ]

  # Drop duplicates
  drop <- NULL
  Dup <- as.character(unique(Biobase::fData(eset)[which(duplicated(Biobase::fData(eset)[, map])), map]))
  Var <- apply(Biobase::exprs(eset), 1, stats::var)
  for (j in Dup) {
    pos <- which(Biobase::fData(eset)[, map] == j)
    drop <- c(drop, pos[-which.max(Var[pos])])
  }
  eset <- eset[-drop, ]

  Biobase::featureNames(eset) <- Biobase::fData(eset)[, map]
  return(eset)
}

# Custom function to expand probesets mapping to multiple genes
#' Title
#'
#' @param eset 
#' @param sep 
#' @param map 
#'
#' @return
#' @noRd
#'
#' @examples
expand_probesets <- function (eset, sep = "///", map="Gene.symbol")
{
  x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(Biobase::fData(eset)[,map]), function(x) strsplit(x, sep))
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(Biobase::fData(eset)[,map]), function(x) strsplit(x, sep))
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  idy <- unlist(sapply(1:length(y), function(i) rep(i, length(y[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  idy <- idy[xx]
  x <- unlist(x)[xx]
  y <- unlist(y)[xx]
  
  eset <- eset[idx, ]
  Biobase::featureNames(eset) <- x
  Biobase::fData(eset)[,map] <- x
  Biobase::fData(eset)$gene <- y
  return(eset)
  
}

#' Title
#'
#' @return
#' @noRd
#'
#' @examples
get_upp_dataset<- function(){
  utils::data("upp", package = "breastCancerUPP")
  get("upp")
}

#' Title
#'
#' @return
#' @noRd
#'
#' @examples

get_transbig_dataset<-function(){
  utils::data("transbig", package = "breastCancerTRANSBIG")
  get("transbig")
}

#' Title
#'
#' @return
#' @noRd
#'
#' @examples
#' 

use_BRCA <- function() {

  # To access gene expression data
  train_expr <- Biobase::exprs(get_upp_dataset())
  test_expr <- Biobase::exprs(get_transbig_dataset())
  # To access feature data
  train_features <- Biobase::fData(get_upp_dataset())
  test_features <- Biobase::fData(get_transbig_dataset())
  # To access clinical data
  train_clinic <- Biobase::pData(get_upp_dataset())
  test_clinic <- Biobase::pData(get_transbig_dataset())

  train <- drop_duplicates(get_upp_dataset())
  train <- expand_probesets(train)
  train <- train[, !is.na(survival::Surv(time = train_clinic$t.rfs, event = train_clinic$e.rfs))] # Drop NAs in survival
  test <- drop_duplicates(get_transbig_dataset())
  test <- expand_probesets(test)
  test <- test[, !is.na(survival::Surv(time = test_clinic$t.rfs, event = test_clinic$e.rfs))] # Drop NAs in survival

  # Determine common probes (Genes)
  int <- intersect(rownames(train), rownames(test))
  train <- train[int, ]
  test <- test[int, ]
  identical(rownames(train), rownames(test))
  # First we will get PAM50 centroids from genefu package
  PAM50Centroids <- genefu::pam50$centroids
  PAM50Genes <- genefu::pam50$centroids.map$probe
  PAM50Genes <- Biobase::featureNames(train)[Biobase::featureNames(train) %in% PAM50Genes]
  # Now we sample 200 random genes from expression matrix
  Non_PAM50Genes <- Biobase::featureNames(train)[!Biobase::featureNames(train) %in% PAM50Genes]
  Non_PAM50Genes <- sample(Non_PAM50Genes, 200, replace = FALSE)
  reduced_set <- c(PAM50Genes, Non_PAM50Genes)
  # Now we get the reduced training and test sets
  train <- train[reduced_set, ]
  test <- test[reduced_set, ]
  Biobase::exprs(train) <- t(apply(Biobase::exprs(train), 1, genefu::rescale, na.rm = TRUE, q = 0.05))
  Biobase::exprs(test) <- t(apply(Biobase::exprs(test), 1, genefu::rescale, na.rm = TRUE, q = 0.05))
  train_expr <- Biobase::exprs(train)
  test_expr <- Biobase::exprs(test)
  # Survival objects
  train_surv <- survival::Surv(time = train_clinic$t.rfs, event = train_clinic$e.rfs)
  test_surv <- survival::Surv(time = test_clinic$t.rfs, event = test_clinic$e.rfs)
  
  return(list(train = list(expression_data = train_expr, phenotype_data = train_surv),
              test =  list(expression_data = test_expr, phenotype_data = test_surv)
              )
         )
}


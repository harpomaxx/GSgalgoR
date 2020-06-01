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
#' rna_luad<- use_rna_luad()
#' TCGA<- rna_luad$TCGA #Access TCGA dataset
#' GSE68465<- rna_luad$GSE68465 #Access GSE68465 dataset
#'
#' #To access gene expression data
#' TCGA_expr<- TCGA$expression_data
#'
#' #To access feature data
#' TCGA_features<- TCGA$feature_data
#'
#' #To access clinical data
#' TCGA_clinic <- TCGA$pheno_data
#'
#' #To get wilkerson centroids
#' WilkCentroids <- rna_luad$WilkCentroids
#' }

use_rna_luad <- function(userdir=""){
  if(userdir == ""){
    dest_dir <- paste0(tempfile(),"/")
  }else{
    dest_dir <- userdir
  }
  if (!dir.exists(dest_dir))
      dir.create(dest_dir,recursive=TRUE)
  dest_file <- paste(dest_dir,"luad_data.rds",sep="")
  if (! file.exists(dest_file))
      utils::download.file("https://bit.ly/luad_data_galgo",dest_file, mode="wb")
  data<-readRDS(dest_file)
}

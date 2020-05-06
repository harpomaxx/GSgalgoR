#' calculate_distance
#'
#' Functions to calculate distance matrices using cpu or gpu computing
#'
#' @param distancetype a \code{character} that can be either \code{'pearson'}, \code{'uncentered'}, \code{'spearman'} or \code{'euclidean'}
#' @param usegpu \code{logical} \code{TRUE} or \code{FALSE}
#' @param x an expression matrix with features as rows and samples as columns
#'
#' @return
#' \code{select_distance(distancetype, usegpu)} assigns global function calculate_distance according to the parameters specified
#'
#' \code{calculate_distance_pearson_gpu(x)} returns columnwise pearson distance calculated using the GPU
#'
#' \code{calculate_distance_pearson_cpu(x)} returns columnwise pearson distance calculated using the CPU
#'
#' \code{calculate_distance_uncentered_gpu(x)} returns columnwise uncentered pearson distance calculated using the GPU
#'
#' \code{calculate_distance_uncentered_cpu(x)} returns columnwise uncentered pearson distance calculated using the CPU
#'
#' \code{calculate_distance_spearman_gpu(x)} returns columnwise spearman distance calculated using the GPU
#'
#' \code{calculate_distance_spearman_cpu(x)} returns columnwise spearman distance calculated using the CPU
#'
#' \code{calculate_distance_euclidean_gpu(x)} returns columnwise euclidean distance calculated using the GPU
#'
#' \code{calculate_distance_euclidean_cpu(x)} returns columnwise euclidean distance calculated using the CPU
#'
#' @author Martin E Guerrero-Gimenez, \email{mguerrero@mendoza-conicet.gob.ar}
#' @name calculate_distance
#' @examples
#' rna_luad<-use_rna_luad()
#' library(Biobase)
#' prm <- exprs(rna_luad$TCGA)
#' select_distance(distancetype= "pearson", usegpu=TRUE) #Will select the calculate_distance_pearson_gpu() function to calculate the distance matrix
#' calculate_distance(prm)
#' calculate_distance_pearson_gpu(prm)
#' calculate_distance_pearson_cpu(prm)
#' calculate_distance_uncentered_gpu(prm)
#' calculate_distance_uncentered_cpu(prm)
#' calculate_distance_spearman_gpu(prm)
#' calculate_distance_spearman_cpu(prm)
#' calculate_distance_euclidean_gpu(prm)
#' calculate_distance_euclidean_cpu(prm)
NULL

#' @rdname calculate_distance
calculate_distance_pearson_gpu <- function(x) {
        x <- gpuR::vclMatrix(x)
        mx <- gpuR::colMeans(x)
        x2 <- gpuR::colMeans(x^2)
        mx2 <- mx^2
        sx <- sqrt(x2 - mx2)
        pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
        # pDist=gpuR::cov(x)/sx%o%sx
        pDist <- as.matrix(pDist)
        as.dist(1 - pDist)
 }

#' @rdname calculate_distance
calculate_distance_spearman_gpu <- function(x) {
        x <- gpuR::vclMatrix(x)
        mx <- gpuR::colMeans(x)
        x2 <- gpuR::colMeans(x^2)
        mx2 <- mx^2
        sx <- sqrt(x2 - mx2)
        pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
        pDist <- as.matrix(pDist)
        # pDist=gpuR::cov(x)/sx%o%sx
        as.dist(1 - pDist)
}

#' @rdname calculate_distance
calculate_distance_uncentered_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        d2 <- gpuR::vclMatrix(1, ncol = 1, nrow = dim(mgpu)[2])
        a1 <- gpuR::tcrossprod(mgpu, mgpu)
        a2 <- mgpu^2 %*% d2
        pDist <- a1 / sqrt(gpuR::tcrossprod(a2, a2))
        pDist <- as.matrix(pDist)
        as.dist(1 - pDist)
}

#' @rdname calculate_distance
calculate_distance_euclidean_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        pDist <- suppressWarnings(gpuR::distance(mgpu, mgpu, method = "euclidean"))
        pDist <- as.matrix(pDist)
        as.dist(pDist)
      }

#' @rdname calculate_distance
calculate_distance_pearson_cpu <-function(x){
    print("pearson....")
    mx=base::colMeans(x)
    x2=base::colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((base::crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    as.dist(1-pDist)
}

#' @rdname calculate_distance
calculate_distance_spearman_cpu <- function(x){
    x=apply(x,2,rank)
    mx=colMeans(x)
    x2=colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    as.dist(1-pDist)
}

#' @rdname calculate_distance
calculate_distance_uncentered_cpu <- function(x){
    mgpu=t(x)
    d2= matrix(1,ncol=1,nrow=dim(mgpu)[2])
    a1=tcrossprod(mgpu,mgpu)
    a2= mgpu^2 %*% d2
    pDist= a1/ sqrt(tcrossprod(a2,a2))
    as.dist(1-pDist)
}

#' @rdname calculate_distance
calculate_distance_euclidean_cpu<-function(x){
    d=stats::dist(t(x),method="euclidean")
}

#' @rdname calculate_distance
select_distance <- function(distancetype = "pearson", usegpu = TRUE ) {
  distancetype <- match.arg(distancetype, c("pearson", "uncentered", "euclidean", "spearman"))

  if (usegpu == FALSE ) {
    computingtype <- "Using CPU for computing"
    # Centered pearson distance
    if (distancetype == "pearson") {
      calculate_distance <- calculate_distance_pearson_cpu
    }
    # Spearman distance
    if (distancetype == "spearman") {
      calculate_distance <- calculate_distance_spearman_cpu
    }

    if (distancetype == "uncentered") {
      calculate_distance <- calculate_distance_uncentered_cpu
    }
    if (distancetype == "euclidean") {
      calculate_distance <- calculate_distance_euclidean_cpu
    }
  } else {
    computingtype <- "Using GPU for computing"
    # Centered pearson distance
    if (distancetype == "pearson") {
     calculate_distance <- calculate_distance_pearson_gpu
    }
    # Spearman distance
    if (distancetype == "spearman") {
     calculate_distance <- calculate_distance_spearman_gpu
    }

    if (distancetype == "uncentered") {
     calculate_distance <- calculate_distance_uncentered_gpu
    }
    if (distancetype == "euclidean") {
     calculate_distance <- calculate_distance_euclidean_gpu
    }
  }
  print(paste(computingtype,distancetype,"distance"))
  calculate_distance
}


#' Distance to centroid classifier
#'
#' Given an \eqn{n x m} matrix of centroids, where \eqn{m} are the prototypic centroids with \eqn{n} features, classify new samples according to the distance to the centroids.
#'
#' @param data a \code{data.frame} of dimensions \eqn{n x p} with the samples to classify, were \eqn{n} are the same set of features as in the centroids
#' @param centroid a \code{data.frame} of dimensions \eqn{n x m}, where each column is a prototypic centroid to classify the samples
#' @param method Character string indicating which method to use to calculate distance to centroid. Options are \code{"pearson"} (default), \code{"kendall"}, or \code{"spearman"}
#'
#' @return Returns a numeric vector of length \eqn{p} with the class assigned to each sample according to the shortest distance to centroid
#' @export
#'
#' @examples
#' rna_luad<-use_rna_luad()
#' library(Biobase)
#'
#' #The expression of the toy datasets are already scaled
#' prm <- exprs(rna_luad$TCGA)
#'
#' #We change the rownames to be gene Symbol insted of Gene Id.
#' rownames(prm)<- fData(rna_luad$TCGA)$gene
#'
#' #Wilkerson's centroids
#' centroids<- rna_luad$WilkCentroids
#'
#' #Extract features from both data.frames
#' inBoth<- Reduce(intersect, list(rownames(prm),rownames(centroids)))
#'
#' #Classify samples
#' Wilk.Class<- classify(prm[inBoth,],centroids[inBoth,])
#' table(Wilk.Class)

classify <- function(data, centroid, method = "pearson") {
    R <- stats::cor(data, centroid, method = method)
    scores <- apply(R, 1, which.max)
}


#' Wrapper function to perform partition around medioids (PAM) for GalgoR
#'
#' In \code{galgoR}, the partition around medioids (PAM) algotithm is the default clustering process used under the evolutionary process.
#'
#' @param c a dissimilarity matix object of type \code{'dist'}
#' @param k positive integer specifying the number of clusters, less than the number of observations
#'
#' @return Returns a \code{'list'} with the value \code{'$cluster'} which contains the cluster assignment of each of the samples evaluated
#' @details The function runs the \code{\link[cluster:pam]{pam}} function of the \code{'cluster'} package with options \code{cluster.only =TRUE}, \code{diss = TRUE}, \code{do.swap=TRUE}, \code{keep.diss=FALSE}, \code{keep.data = FALSE}, \code{pamonce= 2}
#' @references
#' \itemize{
#'   \item Reynolds, A., Richards, G., de la Iglesia, B. and Rayward-Smith, V. (1992) Clustering rules: A comparison of partitioning and hierarchical clustering algorithms; Journal of Mathematical Modelling and Algorithms 5, 475--504. 10.1007/s10852-005-9022-1.
#'   \item Erich Schubert and Peter J. Rousseeuw (2019) Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms; Preprint, \url{(https://arxiv.org/abs/1810.05691)}.
#' }
#' @export
#' @noRd
#' @examples
#' rna_luad<-use_rna_luad()
#'library(Biobase)
#'prm <- exprs(rna_luad$TCGA)
#'Dist= calculate_distance_pearson_cpu(prm)
#'k=4
#'Pam= cluster_algorithm(Dist,k)
#'table(Pam$cluster)

cluster_algorithm <- function(c, k) {
  return(list(cluster = cluster::pam(c, k, cluster.only = TRUE, diss = TRUE, do.swap = TRUE, keep.diss = FALSE, keep.data = FALSE, pamonce = 2)))
}

#' cosine similarity
#'
#' Cosine similarity is a metric of similarity between two non-zero vectors of an inner product space that measures the cosine of the angle between them.
#'Two vectors with the same orientation have a cosine similarity of 1, if they are perpendicular they have a similarity of 0, and if they have oposing directions the cosine similarity is -1, independent of their magnitude.
#'One advantage of cosine similarity is its low-complexity, especially for sparse vectors where only the non-zero dimensions need to be considered, which is a common case in \code{galgoR}.
#'Other names of cosine similarity are Otuska-Orchini similarity when it is applied to binary data, which is the case for \code{galgoR}, where individual solutions represented as strings of 0 and 1 are compared with this metric.
#'
#' @param a,b A string of numbers with equal length. It can also be two binary strings of 0's and 1's
#'
#' @return In practice, the function can return numeric values from -1 to 1 according the vector orientations, where a cosine similarity of 1 implies same orientation of the vectors while -1 imply vector of oposing directions. In the binary application, values range from 0 to 1, where 0 are totally discordant vectors while 1 are identical binary vectors.
#' @export
#' @examples
#' solution1= c(1,0,0,1,0,0,1)
#' solution2= solution1
#' r=cosine(solution1, solution2)
#' #the cosine similarity (r) equals 1
#' solution2= abs(solution1-1)
#' r2=cosine(solution1, solution2)
#' #the cosine similarity (r2) equals 0

cosine <- function(a, b) {
  a %*% b / sqrt(a %*% a * b %*% b)
} # Calculates the cosine distance between two solutions to estimate the mutational rate

#' Function to calculate the centroids of different groups (classes)
#'
#' This function calculates the mean value for each feature of each class to calculate the prototypic centroids of the different groups
#' @param data a scaled gene expression \code{matrix} or \code{data.frame} with samples as columns and features as rows
#' @param class a vector with the samples classes
#'
#' @return returns a \code{data.frame} with the estimated prototypic centroids for each class with the features names as rownames
#' @export
#'
#' @examples
#' rna_luad<-use_rna_luad()
#' library(Biobase)
#' prm <- exprs(rna_luad$TCGA)
#' Dist= calculate_distance_pearson_cpu(prm)
#' k=4
#' Pam= cluster_algorithm(Dist,k)$cluster
#' centroids <- kcentroid(prm,Pam)

kcentroid <- function(data, class) {
  L <- list()
  c <- unique(unlist(class))
  for (i in c) {
    if (sum(unlist(class) == i) > 1) {
      x <- rowMeans(data[, unlist(class) == i])
      L[[i]] <- x
    } else {
      L[[i]] <- data[, unlist(class) == i]
    }
  }
  L <- t(do.call(rbind, L))
}


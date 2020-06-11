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
#' # load example dataset
#' require(iC10TrainingData)
#' require(pamr)
#'
#' data(train.Exp)
#'
#' calculate_distance <- select_distance(distancetype = "pearson", usegpu = FALSE)
#' Dist <- calculate_distance(train.Exp)
#' k <- 4
#' Pam <- cluster_algorithm(Dist, k)
#' table(Pam$cluster)
NULL

#' @rdname calculate_distance
#' @export
calculate_distance_pearson_gpu <- function(x) {
  if (requireNamespace("gpuR", quietly = TRUE)) {
    x <- gpuR::vclMatrix(x)
    mx <- gpuR::colMeans(x)
    x2 <- gpuR::colMeans(x^2)
    mx2 <- mx^2
    sx <- sqrt(x2 - mx2)
    pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
    # pDist=gpuR::cov(x)/sx%o%sx
    pDist <- as.matrix(pDist)
    stats::as.dist(1 - pDist)
  } else {
    message("Package gpuR not available in your platform. Fallback to CPU.")
    calculate_distance_pearson_cpu(x)
  }
}

#' @rdname calculate_distance
#'
calculate_distance_spearman_gpu <- function(x) {
  if (requireNamespace("gpuR", quietly = TRUE)) {
    x <- gpuR::vclMatrix(x)
    mx <- gpuR::colMeans(x)
    x2 <- gpuR::colMeans(x^2)
    mx2 <- mx^2
    sx <- sqrt(x2 - mx2)
    pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
    pDist <- as.matrix(pDist)
    # pDist=gpuR::cov(x)/sx%o%sx
    stats::as.dist(1 - pDist)
  } else {
    message("Package gpuR not available in your platform. Fallback to CPU.")
    calculate_distance_spearman_cpu(x)
  }
}

#' @rdname calculate_distance
#' @export
calculate_distance_uncentered_gpu <- function(x) {
  if (requireNamespace("gpuR", quietly = TRUE)) {
    mgpu <- gpuR::vclMatrix(t(x))
    d2 <- gpuR::vclMatrix(1, ncol = 1, nrow = dim(mgpu)[2])
    a1 <- gpuR::tcrossprod(mgpu, mgpu)
    a2 <- mgpu^2 %*% d2
    pDist <- a1 / sqrt(gpuR::tcrossprod(a2, a2))
    pDist <- as.matrix(pDist)
    stats::as.dist(1 - pDist)
  } else {
    message("Package gpuR not available in your platform. Fallback to CPU.")
    calculate_distance_uncentered_cpu(x)
  }
}

#' @rdname calculate_distance
#' @export
calculate_distance_euclidean_gpu <- function(x) {
  if (requireNamespace("gpuR", quietly = TRUE)) {
    mgpu <- gpuR::vclMatrix(t(x))
    pDist <- suppressWarnings(gpuR::distance(mgpu, mgpu, method = "euclidean"))
    pDist <- as.matrix(pDist)
    stats::as.dist(pDist)
  } else {
    message("Package gpuR not available in your platform. Fallback to CPU.")
    calculate_distance_euclidean_cpu(x)
  }
}
#' @rdname calculate_distance
#' @export
calculate_distance_pearson_cpu <- function(x) {
  # print("pearson....")
  mx <- base::colMeans(x)
  x2 <- base::colMeans(x^2)
  mx2 <- mx^2
  sx <- sqrt(x2 - mx2)
  pDist <- ((base::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
  stats::as.dist(1 - pDist)
}

#' @rdname calculate_distance
#' @export
calculate_distance_spearman_cpu <- function(x) {
  x <- apply(x, 2, rank)
  mx <- colMeans(x)
  x2 <- colMeans(x^2)
  mx2 <- mx^2
  sx <- sqrt(x2 - mx2)
  pDist <- ((crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
  stats::as.dist(1 - pDist)
}

#' @rdname calculate_distance
#' @export
calculate_distance_uncentered_cpu <- function(x) {
  mgpu <- t(x)
  d2 <- matrix(1, ncol = 1, nrow = dim(mgpu)[2])
  a1 <- tcrossprod(mgpu, mgpu)
  a2 <- mgpu^2 %*% d2
  pDist <- a1 / sqrt(tcrossprod(a2, a2))
  stats::as.dist(1 - pDist)
}

#' @rdname calculate_distance
#' @export
calculate_distance_euclidean_cpu <- function(x) {
  d <- stats::dist(t(x), method = "euclidean")
}

#' @rdname calculate_distance
#' @export
select_distance <- function(distancetype = "pearson", usegpu = TRUE) {
  distancetype <- match.arg(distancetype, c("pearson", "uncentered", "euclidean", "spearman"))

  if (usegpu == FALSE) {
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
  message(paste(computingtype, distancetype, "distance"))
  # assign("calculate_distance",calculate_distance,envir=.GlobalEnv)
  # return(invisible(calculate_distance))
  return(calculate_distance)
}

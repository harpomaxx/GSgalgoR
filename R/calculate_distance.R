# Distance type used in the algorithm
# distancetype <- "pearson" # Options are: "pearson","uncentered" (cosine),"spearman","euclidean"
calculate_distance_pearson_gpu <- function(x) {
        x <- gpuR::vclMatrix(x)
        mx <- gpuR::colMeans(x)
        x2 <- gpuR::colMeans(x^2)
        mx2 <- mx^2
        sx <- sqrt(x2 - mx2)
        pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
        # pDist=gpuR::cov(x)/sx%o%sx
        as.dist(1 - pDist)
 }
calculate_distance_spearman_gpu <- function(x) {
        x <- gpuR::vclMatrix(x)
        mx <- gpuR::colMeans(x)
        x2 <- gpuR::colMeans(x^2)
        mx2 <- mx^2
        sx <- sqrt(x2 - mx2)
        pDist <- ((gpuR::crossprod(x, x) / nrow(x)) - (mx %o% mx)) / sx %o% sx
        # pDist=gpuR::cov(x)/sx%o%sx
        as.dist(1 - pDist)
}

calculate_distance_uncentered_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        d2 <- gpuR::vclMatrix(1, ncol = 1, nrow = dim(mgpu)[2])
        a1 <- gpuR::tcrossprod(mgpu, mgpu)
        a2 <- mgpu^2 %*% d2
        pDist <- a1 / sqrt(gpuR::tcrossprod(a2, a2))
        as.dist(1 - pDist)
}
calculate_distance_euclidean_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        d <- suppressWarnings(gpuR::distance(mgpu, mgpu, method = "euclidean"))
        as.dist(d)
      }

calculate_distance_pearson_cpu<- function(x) {
      amap::Dist(x, method = "correlation")
    }

select_distance <- function(distancetype = "pearson", GPU = TRUE ) {
  distancetype <- match.arg(distancetype, c("pearson", "uncentered", "euclidean", "spearman"))

  if (GPU == FALSE ) {
    computingtype <- "using CPU computing"
    print(computingtype)
    calculate_distance <- calculate_distance_pearson_cpu

  } else {
    computingtype <- "using GPU computing"
    print(computingtype)
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
  calculate_distance
  }
}


# Given the distance to the centroids classify the samples
classify <- function(data, centroid, method = "pearson") {
    R <- stats::cor(data, centroid, method = method)
    scores <- apply(R, 1, which.max)
}


# To use partition around medioids (PAM), This is the default for the the current aplication
cluster_algorithm <- function(c, k) {
  return(list(cluster = cluster::pam(c, k, cluster.only = TRUE, diss = TRUE, do.swap = TRUE, keep.diss = FALSE, keep.data = FALSE, pamonce = 2)))
}


# cosine similarity
cosine <- function(a, b) {
  a %*% b / sqrt(a %*% a * b %*% b)
} # Calculates the cosine distance between two solutions to estimate the mutational rate

# Function to calculate the centroids of different groups (classes)
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


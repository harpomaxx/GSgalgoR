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
        pDist <- as.matrix(pDist)
        as.dist(1 - pDist)
 }
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

calculate_distance_uncentered_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        d2 <- gpuR::vclMatrix(1, ncol = 1, nrow = dim(mgpu)[2])
        a1 <- gpuR::tcrossprod(mgpu, mgpu)
        a2 <- mgpu^2 %*% d2
        pDist <- a1 / sqrt(gpuR::tcrossprod(a2, a2))
        pDist <- as.matrix(pDist)
        as.dist(1 - pDist)
}
calculate_distance_euclidean_gpu <- function(x) {
        mgpu <- gpuR::vclMatrix(t(x))
        pDist <- suppressWarnings(gpuR::distance(mgpu, mgpu, method = "euclidean"))
        pDist <- as.matrix(pDist)
        as.dist(pDist)
      }


calculate_distance_pearson_cpu <-function(x){
    print("pearson....")
    mx=base::colMeans(x)
    x2=base::colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((base::crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    as.dist(1-pDist)
}


calculate_distance_spearman_cpu <- function(x){
    x=apply(x,2,rank)
    mx=colMeans(x)
    x2=colMeans(x^2)
    mx2=mx^2
    sx=sqrt(x2-mx2)
    pDist=((crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
    as.dist(1-pDist)
}


calculate_distance_uncentered_cpu <- function(x){
    mgpu=t(x)
    d2= matrix(1,ncol=1,nrow=dim(mgpu)[2])
    a1=tcrossprod(mgpu,mgpu)
    a2= mgpu^2 %*% d2
    pDist= a1/ sqrt(tcrossprod(a2,a2))
    as.dist(1-pDist)
}

calculate_distance_euclidean_cpu<-function(x){
    d=stats::dist(t(x),method="euclidean")
}



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


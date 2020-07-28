#' Distance to centroid classifier function
#'
#' Given an \eqn{n x m} matrix of centroids, where \eqn{m} are the prototypic 
#' centroids with \eqn{n} features, classify new samples according to the 
#' distance to the centroids.
#'
#' @param data a \code{data.frame} of dimensions \eqn{n x p} with the samples 
#' to classify, were \eqn{n} are the same set of features as in the centroids
#' @param centroid a \code{data.frame} of dimensions \eqn{n x m}, where each 
#' column is a prototypic centroid to classify the samples
#' @param method Character string indicating which method to use to calculate 
#' distance to centroid. Options are \code{"pearson"} (default), 
#' \code{"kendall"}, or \code{"spearman"}
#'
#' @return Returns a numeric vector of length \eqn{p} with the class assigned 
#' to each sample according to the shortest distance to centroid
#' @export
#'
#' @examples
#'
#' # load example dataset
#' require(iC10TrainingData)
#' require(pamr)
#'
#' data(train.Exp)
#' data(IntClustMemb)
#' TrainData <- list(x = train.Exp, y = IntClustMemb)
#'
#' # Create prototypic centroids
#' pam <- pamr.train(TrainData)
#' centroids <- pam$centroids
#'
#' Class <- cluster_classify(train.Exp, centroids)
#' table(Class, IntClustMemb)
cluster_classify <- function(data, centroid, method = "pearson") {
    R <- stats::cor(data, centroid, method = method)
    scores <- apply(R, 1, which.max)
}


#' Wrapper function to perform partition around medioids (PAM) for GalgoR
#'
#' In \code{galgoR}, the partition around medioids (PAM) algorithm is the 
#' default clustering process used under the evolutionary process.
#'
#' @param c a dissimilarity matrix object of type \code{'dist'}
#' @param k positive integer specifying the number of clusters, less than the 
#' number of observations
#'
#' @return Returns a \code{'list'} with the value \code{'$cluster'} which 
#' contains the cluster assignment of each of the samples evaluated
#' @details The function runs the \code{\link[cluster:pam]{pam}} function of 
#' the \code{'cluster'} package 
#' with options \code{cluster.only =TRUE}, 
#' \code{diss = TRUE}, \code{do.swap=TRUE}, 
#' \code{keep.diss=FALSE}, \code{keep.data = FALSE}, 
#' \code{pamonce= 2}
#' @references
#' \itemize{
#' \item Reynolds, A., Richards, G., de la Iglesia, B. and Rayward-Smith, V. 
#' (1992) Clustering rules: A comparison of partitioning and 
#' hierarchical clustering algorithms; Journal of Mathematical Modelling and 
#' Algorithms 5, 475--504. 10.1007/s10852-005-9022-1.
#' \item Erich Schubert and Peter J. Rousseeuw (2019) Faster k-Medoids 
#' Clustering: Improving the PAM, CLARA, and CLARANS Algorithms; Preprint, 
#' (\url{https://arxiv.org/abs/1810.05691}).
#' }
#'
#' @export
#' @examples
#' # load example dataset
#' require(iC10TrainingData)
#' require(pamr)
#' data(train.Exp)
#'
#' calculate_distance <- select_distance(distancetype = "pearson", 
#' usegpu = FALSE)
#' Dist <- calculate_distance(train.Exp)
#' k <- 4
#' Pam <- cluster_algorithm(Dist, k)
#' table(Pam$cluster)
cluster_algorithm <- function(c, k) {
    return(list(
        cluster = cluster::pam(
            c,
            k,
            cluster.only = TRUE,
            diss = TRUE,
            do.swap = TRUE,
            keep.diss = FALSE,
            keep.data = FALSE,
            pamonce = 2
        )
    ))
}

#' Function for calculating the cosine similarity
#'
#' Cosine similarity is a metric of similarity between two non-zero vectors 
#' of an inner product space that measures the cosine of the angle between them.
#' Two vectors with the same orientation have a cosine similarity of 1, if 
#' they are perpendicular they have a similarity of 0, and if they have 
#' opposing directions the cosine similarity is -1, independent of their
#' magnitude.
#' One advantage of cosine similarity is its low-complexity, especially for 
#' sparse vectors where only the non-zero dimensions need to be considered, 
#' which is a common case in \code{galgoR}.
#' Other names of cosine similarity are Otuska-Orchini similarity when it is 
#' applied to binary data, which is the case for \code{galgoR}, where 
#' individual solutions represented as strings of 0 and 1 are compared with t
#' his metric.
#'
#' @param a,b A string of numbers with equal length. It can also be two 
#' binary strings of 0's and 1's
#'
#' @return In practice, the function can return numeric values from -1 to 1 
#' according the vector orientations, where a cosine similarity of 1 implies 
#' same orientation of the vectors while -1 imply vector of opposing directions.
#' In the binary application, values range from 0 to 1, where 0 are totally 
#' discordant vectors while 1 are identical binary vectors.
#' @export
#' @examples
#' solution1 <- c(1, 0, 0, 1, 0, 0, 1)
#' solution2 <- solution1
#' r <- cosine_similarity(solution1, solution2)
#' # the cosine similarity (r) equals 1
#' solution2 <- abs(solution1 - 1)
#' r2 <- cosine_similarity(solution1, solution2)
#' # the cosine similarity (r2) equals 0
cosine_similarity <- function(a, b) {
    a %*% b / sqrt(a %*% a * b %*% b)
}   # Calculates the cosine distance between two solutions to estimate the 
    #mutational rate

#' Function to calculate the centroids of different groups (classes)
#'
#' This function calculates the mean value for each feature of each class to 
#' calculate the prototypic centroids of the different groups
#' @param data a scaled gene expression \code{matrix} or \code{data.frame} 
#' with samples as columns and features as rows
#' @param class a vector with the samples classes
#'
#' @return returns a \code{data.frame} with the estimated prototypic centroids 
#' for each class with the features names as rownames
#' @export
#'
#' @examples
#' # load example dataset
#' require(iC10TrainingData)
#' require(pamr)
#'
#' data(train.Exp)
#'
#' calculate_distance <- select_distance(distancetype = "pearson", 
#' usegpu = FALSE)
#' Dist <- calculate_distance(train.Exp)
#' k <- 4
#' Pam <- cluster_algorithm(Dist, k)
#' table(Pam$cluster)
#' centroids <- k_centroids(train.Exp, Pam)
k_centroids <- function(data, class) {
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

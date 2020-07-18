# GalgoR results functions
#' Summary of the non dominated solutions
#'
#' The function uses a \code{'galgo.Obj'} as input an the training dataset to evaluate the non-dominated solutions found by GalgoR
#'
#' @param output An object of class \code{galgo.Obj}
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param OS a \code{survival} object (see \code{\link[survival]{Surv}} function from the \code{\link{survival}} package)
#' @param distancetype a \code{character} that can be either \code{'pearson'}, \code{'uncentered'}, \code{'spearman'} or \code{'euclidean'}
#' @param usegpu \code{logical} \code{TRUE} or \code{FALSE}
#'
#' @return Returns a \code{data.frame} with 5 columns and a number of rows equals to the non-dominated solutions found by GalgoR.
#' The first column has the name of the non-dominated solutions, the second the number of partitions found for each solution \code{(k)}, the third, the number of genes, the fourth the mean silhouette coefficient of the solution and the last columns has the estimated C.Index for each one.
#' @export
#'
#' @examples
#' # load example dataset
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#'
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#'
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 5, population = 15, prob_matrix = expression, OS = OS)
#' non_dominated_summary(
#'     output = output,
#'     OS = OS,
#'     prob_matrix = expression,
#'     distancetype = "pearson",
#'     usegpu = FALSE
#' )
non_dominated_summary <-
    function(output,
             prob_matrix,
             OS,
             distancetype = "pearson",
             usegpu = FALSE) {
        if (!methods::is(output, "galgo.Obj")) {
            stop("object must be of class 'galgo.Obj'")
        }
        output_df <- to_dataframe(output)
        NonDom_solutions <- output_df[output_df$Rank == 1, ]
        calculate_distance <- select_distance(distancetype = distancetype, usegpu = usegpu)
        RESULT <- data.frame(
            solution = as.character(),
            k = as.numeric(),
            ngenes = as.numeric(),
            mean.Silhouette = as.numeric(),
            C.Index = as.numeric(),
            stringsAsFactors = FALSE
        )
        for (i in seq_len(nrow(NonDom_solutions))) {
            name <- rownames(NonDom_solutions)[i]
            genes <- NonDom_solutions[i, "Genes"][[1]]
            k <- NonDom_solutions[i, "k"]
            Sub_matrix <- prob_matrix[genes, ]
            D <- calculate_distance(Sub_matrix)
            true_class <- cluster_algorithm(D, k)
            Centroids <- k_centroids(Sub_matrix, true_class$cluster)
            predicted_class <- cluster_classify(Sub_matrix, Centroids, method = distancetype)
            predicted_class <- as.factor(predicted_class)
            predicted_classdf <- as.data.frame(predicted_class)
            surv_formula <- stats::as.formula("OS~ predicted_class")
            tumortotal <- survival::survfit(surv_formula)
            totalsdf <- survival::survdiff(surv_formula)
            tumortotalpval <- 1 - stats::pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
            tumortotalpval <- format(tumortotalpval, digits = 4)
            coxsimple <- survival::coxph(surv_formula, data = predicted_classdf)
            CI <- survcomp::concordance.index(
                stats::predict(coxsimple, new = predicted_classdf),
                surv.time = OS[, 1],
                surv.event = OS[, 2],
                outx = FALSE
            )$c.index
            # CI <- intsurv::cIndex(risk_score = stats::predict(coxsimple,new=predicted_classdf), time = OS[, 1], event = OS[, 2])["index"]
            mean_Sil <- mean(cluster::silhouette(as.numeric(predicted_class), D)[, 3])
            row <- c(name, k, length(genes), mean_Sil, CI)
            RESULT[nrow(RESULT) + 1, ] <- row
        }
        return(RESULT)
    }


#' Create Centroids
#'
#' This functions create the signature centroids estimated from the GalgoR output and the expression matrix of the training sets.
#'
#' @param output @param output An object of class \code{galgo.Obj}
#' @param solution_names A \code{character} vector with the names of the solutions for which the centroids are to be calculated
#' @param trainset a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param distancetype a \code{character} that can be either \code{'pearson'}, \code{'uncentered'}, \code{'spearman'} or \code{'euclidean'}
#' @param usegpu \code{logical} \code{TRUE} or \code{FALSE}
#'
#' @return Returns a list with the centroid matrix for each of the solutions in \code{solution_names}, where each column represents the prototypic centroid of a subtype and each row the constituents features of the solution signature
#' @export
#'
#' @examples
#' # load example dataset
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#'
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#'
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 5, population = 15, prob_matrix = expression, OS = OS)
#' outputDF <- to_dataframe(output)
#' outputList <- to_list(output)
#'
#' RESULTS <- non_dominated_summary(
#'     output = output, OS = OS,
#'     prob_matrix = expression,
#'     distancetype = "pearson",
#'     usegpu = FALSE
#' )
#' CentroidsList <- create_centroids(output, RESULTS$solution, trainset = expression)
create_centroids <-
    function(output,
             solution_names,
             trainset,
             distancetype = "pearson",
             usegpu = FALSE) {
        calculate_distance <-
            select_distance(distancetype = distancetype, usegpu = usegpu)

        CentroidsList <- list()
        output_df <- to_dataframe(output)
        for (j in solution_names) {
            genes <- output_df[j, "Genes"][[1]]
            k <- output_df[j, "k"]
            name <- j
            Sub_matrix <- trainset[genes, ]

            D <- calculate_distance(Sub_matrix)

            true_class <- cluster_algorithm(D, k)

            Centroids <- k_centroids(Sub_matrix, true_class$cluster)
            CentroidsList[[name]] <- Centroids
        }
        return(CentroidsList)
    }


#' Classify samples from multiple centroids
#'
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param centroid_list a\code{list} with the centroid matrix for each of the signatures to evaluate, where each column represents the prototypic centroid of a subtype and each row the constituents features of the solution signature. The output of \code{\link[galgoR:create_centroids]{create_centroids}} can be used.
#' @param distancetype  a \code{character} that can be either \code{'pearson'} (default), \code{'spearman'} or \code{'kendall'}.
#'
#' @return Returns a \code{data.frame} with the classes assigned to each sample in each signature, were samples are a rows and signatures in columns
#' @export
#'
#' @examples
#' # load example dataset
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#'
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#'
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 5, population = 15, prob_matrix = expression, OS = OS)
#' outputDF <- to_dataframe(output)
#' outputList <- to_list(output)
#'
#' RESULTS <- non_dominated_summary(
#'     output = output, OS = OS,
#'     prob_matrix = expression,
#'     distancetype = "pearson",
#'     usegpu = FALSE
#' )
#' CentroidsList <- create_centroids(output, RESULTS$solution, trainset = expression)
#' classes <- classify_multiple(prob_matrix = expression, centroid_list = CentroidsList)
classify_multiple <-
    function(prob_matrix,
             centroid_list,
             distancetype = "pearson") {
        classes <-
            matrix(rep(NA, ncol(prob_matrix) * length(centroid_list)), ncol = length(centroid_list))
        as.data.frame <- classes
        colnames(classes) <- names(centroid_list)
        rownames(classes) <- colnames(prob_matrix)

        for (i in seq_len(length(centroid_list))) {
            centroids <- centroid_list[[i]]
            name <- names(centroid_list)[i]
            genes <- rownames(centroids)
            k <- ncol(centroids)
            Sub_matrix <- prob_matrix[genes, ]

            predicted_class <-
                cluster_classify(Sub_matrix, centroids, method = distancetype)

            predicted_class <- as.factor(predicted_class)
            classes[, name] <- predicted_class
        }

        return(classes)
    }

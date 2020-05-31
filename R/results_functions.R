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
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 2, population = 3, prob_matrix = TCGA_expr, OS = OS)
#' non_dominated_summary(
#'   output = output, OS = OS,
#'   prob_matrix = TCGA_expr,
#'   distancetype = "pearson",
#'   usegpu = FALSE
#' )
non_dominated_summary <- function(output, prob_matrix, OS, distancetype = "pearson", usegpu = FALSE) {
  if (!methods::is(output, "galgo.Obj")) {
    stop("object must be of class 'galgo.Obj'")
  }
  output_df <- toDataFrame(output)
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

  for (i in 1:nrow(NonDom_solutions)) {
    name <- rownames(NonDom_solutions)[i]
    genes <- NonDom_solutions[i, "Genes"][[1]]
    k <- NonDom_solutions[i, "k"]
    Sub_matrix <- prob_matrix[genes, ]

    D <- calculate_distance(Sub_matrix)

    true_class <- cluster_algorithm(D, k)

    Centroids <- kcentroid(Sub_matrix, true_class$cluster)

    predicted_class <- classify(Sub_matrix, Centroids, method = distancetype)

    predicted_class <- as.factor(predicted_class)
    predicted_classdf <- as.data.frame(predicted_class)


    surv_formula <- stats::as.formula("OS~ predicted_class")
    tumortotal <- survival::survfit(surv_formula)
    totalsdf <- survival::survdiff(surv_formula)
    tumortotalpval <- 1 - stats::pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
    tumortotalpval <- format(tumortotalpval, digits = 4)

    coxsimple <- survival::coxph(surv_formula, data = predicted_classdf)

    # CI=intsurv::cIndex(stats::predict(coxsimple),surv.time=OS[,1],surv.event=OS[,2],outx=FALSE)$c.index
    CI <- intsurv::cIndex(risk_score = stats::predict(coxsimple), time = OS[, 1], event = OS[, 2])["index"]
    mean_Sil <- mean(cluster::silhouette(as.numeric(predicted_class), D)[, 3])

    row <- c(name, k, length(genes), mean_Sil, CI)
    RESULT[nrow(RESULT) + 1, ] <- row
    # print(row)
  }
  return(RESULT)
}


#' Create Centroids
#'
#' This functions create the signature centroids estimated from the GalgoR output and the expression matrix of the training sets.
#'
#' @param output @param output An object of class \code{galgo.Obj}
#' @param solution.names A \code{character} vector with the names of the solutions for which the centroids are to be calculated
#' @param train.set a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param distancetype a \code{character} that can be either \code{'pearson'}, \code{'uncentered'}, \code{'spearman'} or \code{'euclidean'}
#' @param usegpu \code{logical} \code{TRUE} or \code{FALSE}
#'
#' @return Returns a list with the centroid matrix for each of the solutions in \code{solution.names}, where each column represents the prototypic centroid of a subtype and each row the constituents features of the solution signature
#' @export
#'
#' @examples
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 2, population = 3, prob_matrix = TCGA_expr, OS = OS)
#' RESULTS <- non_dominated_summary(
#'   output = output, OS = OS,
#'   prob_matrix = TCGA_expr,
#'   distancetype = "pearson",
#'   usegpu = FALSE
#' )
#' CentroidsList <- create_centroids(output, RESULTS$solution, train.set = TCGA_expr)
create_centroids <- function(output, solution.names, train.set, distancetype = "pearson", usegpu = FALSE) {
  calculate_distance <- select_distance(distancetype = distancetype, usegpu = usegpu)

  CentroidsList <- list()
  output_df <- toDataFrame(output)
  for (j in solution.names) {
    genes <- output_df[j, "Genes"][[1]]
    k <- output_df[j, "k"]
    name <- j
    Sub_matrix <- train.set[genes, ]

    D <- calculate_distance(Sub_matrix)

    true_class <- cluster_algorithm(D, k)

    Centroids <- kcentroid(Sub_matrix, true_class$cluster)
    CentroidsList[[name]] <- Centroids
  }
  return(CentroidsList)
}


#' Classify samples from multiple centroids
#'
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param centroid._list a\code{list} with the centroid matrix for each of the signatures to evaluate, where each column represents the prototypic centroid of a subtype and each row the constituents features of the solution signature. The output of \code{\link[galgoR:create_centroids]{create_centroids}} can be used.
#' @param distancetype  a \code{character} that can be either \code{'pearson'} (default), \code{'spearman'} or \code{'kendall'}.
#'
#' @return Returns a \code{data.frame} with the classes assigned to each sample in each signature, were samples are a rows and signatures in columns
#' @export
#'
#' @examples
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 2, population = 3, prob_matrix = TCGA_expr, OS = OS)
#' RESULTS <- non_dominated_summary(
#'   output = output, OS = OS,
#'   prob_matrix = TCGA_expr,
#'   distancetype = "pearson",
#'   usegpu = FALSE
#' )
#' CentroidsList <- create_centroids(output, RESULTS$solution, train.set = TCGA_expr)
#' TCGA_classes <- classify_multiple(prob_matrix = TCGA_expr, centroid._list = CentroidsList)
classify_multiple <- function(prob_matrix, centroid._list, distancetype = "pearson") {
  classes <- matrix(rep(NA, ncol(prob_matrix) * length(centroid._list)), ncol = length(centroid._list))
  as.data.frame <- classes
  colnames(classes) <- names(centroid._list)
  rownames(classes) <- colnames(prob_matrix)

  for (i in 1:length(centroid._list)) {
    centroids <- centroid._list[[i]]
    name <- names(centroid._list)[i]
    genes <- rownames(centroids)
    k <- ncol(centroids)
    Sub_matrix <- prob_matrix[genes, ]

    predicted_class <- classify(Sub_matrix, centroids, method = distancetype)

    predicted_class <- as.factor(predicted_class)
    classes[, name] <- predicted_class
  }

  return(classes)
}

#' Plot pareto front
#'
#' @param output An object of class \code{galgo.Obj}
#'
#' @return This function returns a scatterplot showing the solutions found by Galgo accross all generations in the solution space, where the Silhouette Fitness is in the x-axis and the survival fitness in the y-axis.
#' A line is drawn over all non-dominated solutions showing the estimated Pareto front
#' @export
#'
#' @examples
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 2, population = 3, prob_matrix = TCGA_expr, OS = OS)
#' plot_pareto(output)
plot_pareto <- function(output) {
  SC.Fit <- Surv.Fit <- Gen <- NULL
  PARETO <- output@ParetoFront

  for (i in 1:length(PARETO)) {
    PARETO[[i]] <- as.data.frame(PARETO[[i]])
    PARETO[[i]]$Gen <- i
  }
  PARETO <- do.call(rbind, PARETO)
  colnames(PARETO) <- c("SC.Fit", "Surv.Fit", "Gen")

  output_df <- toDataFrame(output) # Transform output to dataframe
  output_df <- output_df[output_df$Rank == 1, ]
  output_df <- output_df[order(output_df$SC.Fit), ]

  PlotPareto <- ggplot2::ggplot(PARETO, ggplot2::aes(x = SC.Fit, y = Surv.Fit, colour = Gen)) + ggplot2::geom_point() + ggplot2::theme_bw()
  PlotPareto <- PlotPareto + ggplot2::geom_point(ggplot2::aes(x = SC.Fit, y = Surv.Fit), colour = "black", size = 2, output_df) + ggplot2::geom_line(ggplot2::aes(x = SC.Fit, y = Surv.Fit), colour = "black", output_df)
  PlotPareto <- PlotPareto + ggplot2::ggtitle("Galgo run Pareto front")
  print(PlotPareto)
}


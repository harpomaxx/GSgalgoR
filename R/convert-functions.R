#' Convert galgo.Obj to list
#'
#' The current function transforms a \code{galgo.Obj} to a \code{list}
#'
#' @param output An object of class \code{galgo.Obj}
#'
#' @return The current function restructurates a \code{galgo.Obj} to a more easy to understand an use \code{list}. This output is particularly useful if one wants to select a given solution and use its outputs in a new classifier. The output of type \code{list} has a length equals to the number of solutions obtained by the \code{\link[galgoR:galgo]{galgo}} algorithm.
#'
#' Basically this output is a list of lists, where each element of the output is named after the solution's name (\code{solution.n}, where \code{n} is the number assigned to that solution), and inside of it, it has all the constituents for that given solution with the following structure:
#' \enumerate{
#' \item \strong{output$solution.n$Genes}: A vector of the features included in the solution
#' \item \strong{output$solution.n$k}: The number of partitions found in that solution
#' \item \strong{output$solution.n$SC.Fit}: The average silhouette coefficient of the partitions found
#' \item \strong{output$solution.n$Surv.Fit}: The survival fitness value
#' \item \strong{output$solution.n$Rank}: The solution rank
#' \item \strong{CrowD}: The solution crowding distance related to the rest of the solutions
#' }
#' @export
#' @author Martin E Guerrero-Gimenez, \email{mguerrero@mendoza-conicet.gob.ar}
#'
#' @examples
#' \dontrun{
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 10, population = 30, prob_matrix = TCGA_expr, OS = OS)
#' outputList <- to_list(output)
#' }
#'
to_list <- function(output) {
  if (!methods::is(output, "galgo.Obj")) {
    stop("object must be of class 'galgo.Obj'")
  }
  OUTPUT <- list()
  Genes <- colnames(output@Solutions)[1:(ncol(output@Solutions) - 5)]
  for (i in 1:nrow(output@Solutions)) {
    Sol <- paste("Solution", i, sep = ".")
    OUTPUT[[Sol]] <- list()
    OUTPUT[[Sol]][["Genes"]] <- Genes[as.logical(output@Solutions[i, 1:length(Genes)])]
    OUTPUT[[Sol]]["k"] <- output@Solutions[i, "k"]
    OUTPUT[[Sol]]["SC.Fit"] <- output@Solutions[i, length(Genes) + 2]
    OUTPUT[[Sol]]["Surv.Fit"] <- output@Solutions[i, length(Genes) + 3]
    OUTPUT[[Sol]]["rank"] <- output@Solutions[i, "rnkIndex"]
    OUTPUT[[Sol]]["CrowD"] <- output@Solutions[i, "CrowD"]
  }
  return(OUTPUT)
}

#' Convert galgo.Obj to data.frame
#'
#' The current function transforms a \code{galgo.Obj} to a \code{data.frame}
#'
#' @param output An object of class \code{galgo.Obj}
#'
#' @return The current function restructurates a \code{galgo.Obj} to a more easy to understand an use \code{data.frame}. The output \code{data.frame} has \eqn{ m x n} dimensions, were the rownames (\eqn{m}) are the solutions obtained by the \code{\link[galgoR:galgo]{galgo}} algorithm. The columns has the following structure:
#' \enumerate{
#'  \item \strong{Genes}: The features included in each solution in form of a \code{list}
#'  \item \strong{k}: The number of partitions found in that solution
#'  \item \strong{SC.Fit}: The average silhouette coefficient of the partitions found
#'  \item \strong{Surv.Fit}: The survival fitness value
#'  \item \strong{Rank}: The solution rank
#'  \item \strong{CrowD}: The solution crowding distance related to the rest of the solutions
#' }
#' @export
#'
#' @author Martin E Guerrero-Gimenez, \email{mguerrero@mendoza-conicet.gob.ar}
#'
#' @examples
#' \dontrun{
#' # Load data
#' rna_luad <- use_rna_luad()
#' TCGA_expr <- rna_luad$TCGA$expression_matrix
#' TCGA_clinic <- rna_luad$TCGA$pheno_data
#' OS <- survival::Surv(time = TCGA_clinic$time, event = TCGA_clinic$status)
#'
#' # Run galgo
#' output <- galgoR::galgo(generations = 10, population = 30, prob_matrix = TCGA_expr, OS = OS)
#' outputDF <- to_dataframe(output)
#' }
to_dataframe <- function(output) {
  if (!methods::is(output, "galgo.Obj")) {
    stop("object must be of class 'galgo.Obj'")
  }
  Genes <- colnames(output@Solutions)[1:(ncol(output@Solutions) - 5)]
  ListGenes <- list()
  for (i in 1:nrow(output@Solutions)) {
    ListGenes[[i]] <- list()
    ListGenes[[i]] <- Genes[as.logical(output@Solutions[i, 1:length(Genes)])]
  }
  
  OUTPUT <- data.frame(Genes = I(ListGenes), k = output@Solutions[, "k"], SC.Fit = output@Solutions[, ncol(output@Solutions) - 3], Surv.Fit = output@Solutions[, ncol(output@Solutions) - 2], Rank = output@Solutions[, "rnkIndex"], CrowD = output@Solutions[, "CrowD"])
  rownames(OUTPUT) <- paste("Solutions", 1:nrow(output@Solutions), sep = ".")
  return(OUTPUT)
}
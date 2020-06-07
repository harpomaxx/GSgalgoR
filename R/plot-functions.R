
#' Plot pareto front from an galgo.Obj
#'
#' @param output An object of class \code{galgo.Obj}
#'
#' @return This function returns a scatterplot showing the solutions found by Galgo accross all generations in the solution space, where the Silhouette Fitness is in the x-axis and the survival fitness in the y-axis.
#' A line is drawn over all non-dominated solutions showing the estimated Pareto front
#' @export
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
#' plot_pareto(output)
#' }
plot_pareto <- function(output) {
  SC.Fit <- Surv.Fit <- Gen <- NULL
  PARETO <- output@ParetoFront
  
  for (i in 1:length(PARETO)) {
    PARETO[[i]] <- as.data.frame(PARETO[[i]])
    PARETO[[i]]$Gen <- i
  }
  PARETO <- do.call(rbind, PARETO)
  colnames(PARETO) <- c("SC.Fit", "Surv.Fit", "Gen")
  
  output_df <- to_dataframe(output) # Transform output to dataframe
  output_df <- output_df[output_df$Rank == 1, ]
  output_df <- output_df[order(output_df$SC.Fit), ]
  
  PlotPareto <- ggplot2::ggplot(PARETO, ggplot2::aes(x = SC.Fit, y = Surv.Fit, colour = Gen)) + ggplot2::geom_point() + ggplot2::theme_bw()
  PlotPareto <- PlotPareto + ggplot2::geom_point(ggplot2::aes(x = SC.Fit, y = Surv.Fit), colour = "black", size = 2, output_df) + ggplot2::geom_line(ggplot2::aes(x = SC.Fit, y = Surv.Fit), colour = "black", output_df)
  PlotPareto <- PlotPareto + ggplot2::ggtitle("Galgo run Pareto front")
  print(PlotPareto)
}

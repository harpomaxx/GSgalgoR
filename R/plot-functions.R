#' Plot pareto front from an galgo.Obj
#'
#' @param output An object of class \code{galgo.Obj}
#'
#' @return This function returns a scatterplot showing the solutions found by Galgo accross all generations in the solution space, where the Silhouette Fitness is in the x-axis and the survival fitness in the y-axis.
#' A line is drawn over all non-dominated solutions showing the estimated Pareto front
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
#' plot_pareto(output)
plot_pareto <- function(output) {
    SC.Fit <- Surv.Fit <- Gen <- NULL
    PARETO <- output@ParetoFront

    for (i in seq_len(length(PARETO))) {
        PARETO[[i]] <- as.data.frame(PARETO[[i]])
        PARETO[[i]]$Gen <- i
    }
    PARETO <- do.call(rbind, PARETO)
    colnames(PARETO) <- c("SC.Fit", "Surv.Fit", "Gen")

    output_df <- to_dataframe(output) # Transform output to dataframe
    output_df <- output_df[output_df$Rank == 1, ]
    output_df <- output_df[order(output_df$SC.Fit), ]

    PlotPareto <-
        ggplot2::ggplot(PARETO, ggplot2::aes(x = SC.Fit, y = Surv.Fit, colour = Gen)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw()
    PlotPareto <-
        PlotPareto + ggplot2::geom_point(
            ggplot2::aes(x = SC.Fit, y = Surv.Fit),
            colour = "black",
            size = 2,
            output_df
        ) + ggplot2::geom_line(ggplot2::aes(x = SC.Fit, y = Surv.Fit), colour = "black", output_df)
    PlotPareto <-
        PlotPareto + ggplot2::ggtitle("Galgo run Pareto front")
    print(PlotPareto)
}

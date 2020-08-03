context("plot-functions")

test_that("plot_pareto return a Grob Object", {
    
    library(breastCancerTRANSBIG)
    data(transbig)
    Train <- transbig
    expression <- Biobase::exprs(Train)
    clinical <- Biobase::pData(Train)
    OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
    
    expression <- expression[sample(seq_len(nrow(expression)), 100), ]
    expression <- t(scale(t(expression)))
    output <- GSgalgoR::galgo(generations = 2, population = 3, prob_matrix = expression, OS = OS, 
                            distancetype = "pearson", verbose = 1)
    plot_output<-plot_pareto(output)
    expect_is(plot_output,"ggplot")
    
})

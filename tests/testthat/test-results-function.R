context("distance-functions")

test_that("non_dominated_summary returns a dataframe", {
    
    library(breastCancerTRANSBIG)
    data(transbig)
    Train <- transbig
    expression <- Biobase::exprs(Train)
    clinical <- Biobase::pData(Train)
    OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)

    expression <- expression[sample(seq_len(nrow(expression)), 100), ]
    expression <- t(scale(t(expression)))
    output <- galgoR::galgo(generations = 2, population = 3, prob_matrix = expression, OS = OS, 
                            distancetype = "pearson", verbose = 1)
    summary_results<-non_dominated_summary(output = output,OS = OS,
                                           prob_matrix = expression, 
                                           distancetype = "pearson")
    expect_is(summary_results,"data.frame")
    
})

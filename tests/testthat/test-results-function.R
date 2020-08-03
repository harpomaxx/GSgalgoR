context("results-functions")

test_that("classifiy_multiple returns a list", {
    
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
    summary_results<-non_dominated_summary(output = output,OS = OS,
                                           prob_matrix = expression, 
                                           distancetype = "pearson")
    centroids_list <- create_centroids(output, summary_results$solution, trainset = expression)
    classes <- classify_multiple(prob_matrix = expression, centroid_list = centroids_list)
    expect_is(classes,"matrix")
})

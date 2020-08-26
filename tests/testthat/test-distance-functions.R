context("distance-functions")

test_that("pearson distance works", {
  library(breastCancerTRANSBIG)
  data(transbig)
  Train <- transbig
  expression <- Biobase::exprs(Train)
  clinical <- Biobase::pData(Train)
  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
  expression <- t(scale(t(expression)))
  output <- GSgalgoR::galgo(generations = 2, population = 5, prob_matrix = expression, OS = OS, 
                          distancetype = "pearson", verbose =1 )
  expect_is(output,"galgo.Obj")
  
})


test_that("spearman distance works", {
  
  library(breastCancerTRANSBIG)
  data(transbig)
  Train <- transbig
  expression <- Biobase::exprs(Train)
  clinical <- Biobase::pData(Train)
  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
  expression <- t(scale(t(expression)))
  output <- GSgalgoR::galgo(generations = 2, population = 3, prob_matrix = expression, OS = OS, 
                          distancetype = "spearman", verbose = 1)
  expect_is(output,"galgo.Obj")
  
})


test_that("euclidean distance works", {
  
  library(breastCancerTRANSBIG)
  data(transbig)
  Train <- transbig
  expression <- Biobase::exprs(Train)
  clinical <- Biobase::pData(Train)
  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
  expression <- t(scale(t(expression)))
  output <- GSgalgoR::galgo(generations = 2, population = 5, prob_matrix = expression, OS = OS, 
                          distancetype = "euclidean", verbose = 1)
  expect_is(output,"galgo.Obj")
  
})

test_that("uncentered distance works", {
  
  library(breastCancerTRANSBIG)
  data(transbig)
  Train <- transbig
  expression <- Biobase::exprs(Train)
  clinical <- Biobase::pData(Train)
  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
  
  # We will use a reduced dataset for the example
  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
  
  # Now we scale the expression matrix
  expression <- t(scale(t(expression)))
  
  # Run galgo
  output <- GSgalgoR::galgo(generations = 2, population = 5, prob_matrix = expression, OS = OS, 
                          distancetype = "uncentered", verbose = 1)
  expect_is(output,"galgo.Obj")
  
})

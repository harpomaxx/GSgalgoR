context("convert-functions")

test_that("galgo.Obj can be converted to list", {
  
  # load example dataset
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
  output <- galgoR::galgo(generations = 2, population = 5, prob_matrix = expression, OS = OS, 
                          distancetype = "uncentered")
  output_dataframe <- to_list(output)
  expect_is(output_dataframe,"list")
  
})



test_that("galgo.Obj can be converted to data.frame", {
  
  # load example dataset
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
  output <- galgoR::galgo(generations = 3, population = 5, prob_matrix = expression, OS = OS, 
                          distancetype = "uncentered")
  output_dataframe <- to_dataframe(output)
  expect_is(output_dataframe,"data.frame")
  
})
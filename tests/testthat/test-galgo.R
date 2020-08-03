context("galgo function")

test_that("a valid galgo.Obj is generated", {
  
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
   output <- GSgalgoR::galgo(generations = 3, population = 5, prob_matrix = expression, OS = OS)
   expect_is(output,"galgo.Obj")
  
})


test_that("a valid galgo.Obj is generated using verbose 1", {
   
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
   output <- GSgalgoR::galgo(generations = 3, 
                           population = 5, 
                           prob_matrix = expression,
                           OS = OS,
                           verbose = 1)
   expect_is(output,"galgo.Obj")
   
})

test_that("a valid galgo.Obj is generated using verbose 0", {
   
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
   output <- GSgalgoR::galgo(generations = 3, 
                           population = 5, 
                           prob_matrix = expression,
                           OS = OS,
                           verbose = 0)
   expect_is(output,"galgo.Obj")
   
})

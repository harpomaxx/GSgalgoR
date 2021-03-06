context("callback functions")


test_that("final population is saved in /tmp/", {
  set.seed(29042010)
  library(breastCancerTRANSBIG)
  data(transbig)
  Train <- transbig
  expression <- Biobase::exprs(Train)
  clinical <- Biobase::pData(Train)
  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
  expression <- t(scale(t(expression)))
  tmpdir<-tempdir()
  output <- GSgalgoR::galgo(res_dir = tmpdir,
                          generations = 2, 
                          population = 3, 
                          prob_matrix = expression, 
                          OS = OS,
                          verbose = 0,
                          end_galgo_callback = callback_base_save_pop_final)
  expect_true(file.exists(file.path(tmpdir,"final.rda")))
})

#test_that("partial population is saved in /tmp/", {
#  # load example dataset
#  library(breastCancerTRANSBIG)
#  data(transbig)
#  Train <- transbig
#  expression <- Biobase::exprs(Train)
#  clinical <- Biobase::pData(Train)
#  OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#  expression <- expression[sample(seq_len(nrow(expression)), 100), ]
#  expression <- t(scale(t(expression)))
#  output <- GSgalgoR::galgo(res_dir = "/tmp/",
#                          generations = 3, 
#                          population = 3, 
#                          prob_matrix = expression, 
#                          OS = OS,
#                          verbose = 0,
#                          end_gen_callback = callback_base_save_pop_partial )
#  expect_true(file.exists("/tmp/2.rda"))
#})



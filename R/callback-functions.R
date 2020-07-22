#' A base callback function that returns a galgo.Obj
#' @param userdir the default directory used by `galgo()` to store files
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo across all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#'
#' @return an object of class galgo
#' @export 
#' @usage  callback_base_return_pop (userdir = "",generation, pop_pool, 
#' pareto, prob_matrix, current_time)
#' @examples
#' # load example dataset
#'
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo with base_return_pop_callback assigned to the end_galgo_callback
#' # hook-point
#' # By using this callback galgo() return a `galgo,Obj` object.
#' output <- galgoR::galgo(generations = 5,
#' population = 15,
#' prob_matrix = expression,
#' OS = OS,
#' end_galgo_callback = callback_base_return_pop
#' )
#'
#'
callback_base_return_pop <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        colnames(pop_pool)[seq_len(ncol(pop_pool) - 5)] <-
            rownames(prob_matrix)
        output <- list(Solutions = pop_pool, ParetoFront = pareto)
        output <-
            galgo.Obj(
                Solutions = output$Solutions,
                ParetoFront = output$ParetoFront
            )
        return(output)
    }

#' A base callback function that saves galgo.Obj
#'
#'
#'
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param prefix a prefix used for the file name
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo across all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#'
#' @return an object of class galgo
#' 
#' @usage  callback_base_save_pop (userdir = "",generation, pop_pool, 
#' pareto, prob_matrix, current_time)
#'
#' @examples
#' @noRd
callback_base_save_pop <-
    function(userdir = "",
            prefix,
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        if (userdir == "") {
            directory <- paste0(tempdir(), "/")
        } else {
            directory <- userdir
        }
        if (!dir.exists(directory)) {
            dir.create(directory, recursive = TRUE)
        }
        colnames(pop_pool)[seq_len(ncol(pop_pool) - 5)] <-
            rownames(prob_matrix)
        output <- list(Solutions = pop_pool, ParetoFront = pareto)
        output <-
            galgo.Obj(Solutions = output$Solutions,
                    ParetoFront = output$ParetoFront)
        filename <- paste0(directory, prefix, ".rda")
        save(file = filename, output)
        return(output)
    }

#' A callback for saving partial galgo.Obj (every 2 generations)
#'
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo accross all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#'
#' @return an object of class galgo
#'
#' @usage  callback_base_save_pop_partial (userdir = "",generation,
#' pop_pool, pareto, prob_matrix, current_time)
#'
#' @examples
#' @noRd
callback_base_save_pop_partial <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        if (userdir == "") {
            directory <- paste0(tempdir(), "/")
        } else {
            directory <- userdir
        }
        if (!dir.exists(directory)) {
            dir.create(directory, recursive = TRUE)
        }
        if (generation %% 2 == 0) {
            callback_base_save_pop(
                userdir = directory,
                prefix = generation,
                generation,
                pop_pool,
                pareto,
                prob_matrix,
                current_time
            )
        }
    }

#' A callback for saving final galgo.Obj
#'
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo accross all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#'
#' @return
#'
#' @usage  callback_base_save_pop_final <- function(userdir = "",generation, 
#' pop_pool, pareto, prob_matrix, current_time)
#'
#' @examples
#' @noRd
callback_base_save_pop_final <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        if (userdir == "") {
            directory <- paste0(tempdir(), "/")
        } else {
            directory <- userdir
        }
        if (!dir.exists(directory)) {
            dir.create(directory, recursive = TRUE)
        }
        output <-
            callback_base_save_pop(
                userdir = directory,
                prefix = "final",
                generation,
                pop_pool,
                pareto,
                prob_matrix,
                current_time
            )
        message(paste0("final population saved in final.rda in ", directory))
        return(output)
    }


#' Print basic info per generation
#'
#' @param userdir the default directory used by `galgo()` to store files
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo across all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#' @return Nothing.
#'
#' @export
#' @usage  callback_base_report (userdir, generation, pop_pool, 
#' pareto, prob_matrix, current_time)
#' @examples
#' # load example dataset
#'
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo with base_report_callback assigned to the report_callback
#' # hook-point
#' galgoR::galgo(generations = 5,
#' population = 15,
#' prob_matrix = expression,
#' OS = OS,
#' report_callback = callback_base_report
#' )


callback_base_report <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        chrom_length <- nrow(prob_matrix)
        message(paste0("Generation ", generation, " Non-dominated solutions:"))
        print(pop_pool[pop_pool[, "rnkIndex"] == 1,
                    (chrom_length + 1):(chrom_length + 5)])
        # print(Sys.time()- start_time)
    }

#' Print minimal information to the user about galgo execution.
#'
#' The main idea behind this callback function is to provide some feedback to the user about galgo execution.
#' No other relevant information is shown
#'
#' @param userdir the default directory used by `galgo()` to store files
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo across all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#' @return Nothing.
#'
#' @export
#' 
#' @usage  callback_no_report (userdir = "",generation, pop_pool, 
#' pareto, prob_matrix, current_time)
#' @examples
#' # load example dataset
#'
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo with no_report_callback assigned to the report_callback
#' # hook-point
#' galgoR::galgo(generations = 5,
#' population = 15,
#' prob_matrix = expression,
#' OS = OS,
#' report_callback = callback_no_report
#' )
callback_no_report <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
        if (generation %% 5 == 0) {
            cat("*")
        } else {
            cat(".")
        }
    }


#' A default call_back function that does nothing.
#'
#' @param userdir the default directory used by \code{galgo()} to store files
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo across all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#'
#' @return Nothing
#' 
#' @usage  callback_default (userdir = "",generation, pop_pool, 
#' pareto, prob_matrix, current_time)
#' @export
#' @examples
#' # load example dataset
#'
#' library(breastCancerTRANSBIG)
#' data(transbig)
#' Train <- transbig
#' rm(transbig)
#' expression <- Biobase::exprs(Train)
#' clinical <- Biobase::pData(Train)
#' OS <- survival::Surv(time = clinical$t.rfs, event = clinical$e.rfs)
#' # We will use a reduced dataset for the example
#' expression <- expression[sample(1:nrow(expression), 100), ]
#'
#' # Now we scale the expression matrix
#' expression <- t(scale(t(expression)))
#'
#' # Run galgo with default_callback assigned to all the hook-points
#'
#' galgoR::galgo(generations = 5,
#' population = 15,
#' prob_matrix = expression,
#' OS = OS,
#' start_galgo_callback = callback_default, # When Galgo is about to start.
#' end_galgo_callback = callback_default,   # When Galgo is about to finish.
#' start_gen_callback = callback_default,   # At the beginning of each generation/iteration.
#' end_gen_callback = callback_default,     # At the end of each generation/iteration.
# 'report_callback = callback_default,      # In the middle of the generation
#' )
#'
#'
callback_default <-
    function(userdir = "",
            generation,
            pop_pool,
            pareto,
            prob_matrix,
            current_time) {
    }

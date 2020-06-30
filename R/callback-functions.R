#' A base callback function that returns a galgo.Obj
#' @param userdir
#' @param prefix
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return an objet of class galgo
#'
#'
#' @examples
#' @noRd
base_return_pop_callback <-
    function(userdir = "",
             prefix,
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
#' @param userdir
#' @param prefix
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return an objet of class galgo
#'
#'
#' @examples
#' @noRd
base_save_pop_callback <-
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
            galgo.Obj(
                Solutions = output$Solutions,
                ParetoFront = output$ParetoFront
            )
        filename <- paste0(directory, prefix, ".rda")
        save(file = filename, output)
        return(output)
    }

#' A callback for daving partial galgo.Obj (every 2 generations)
#'
#' @param userdir
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
base_save_pop_partial_callback <-
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
            base_save_pop_callback(
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
#' @param userdir
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
base_save_pop_final_callback <-
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
            base_save_pop_callback(
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

#' Title
#' Print basic info per generation
#'
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
base_report_callback <-
    function(userdir = "",
             generation,
             pop_pool,
             pareto,
             prob_matrix,
             current_time) {
        chrom_length <- nrow(prob_matrix)
        message(paste0("Generation ", generation, " Non-dominated solutions:"))
        print(pop_pool[pop_pool[, "rnkIndex"] == 1, (chrom_length + 1):(chrom_length + 5)])
        # print(Sys.time()- start_time)
    }

#' Title
#'
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
no_report_callback <-
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

#' Title
#'
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
base_start_gen_callback <-
    function(userdir = "",
             generation,
             pop_pool,
             pareto,
             prob_matrix,
             current_time) {
        # start_time <- Sys.time()
        
    }

#' Title
#'
#' @param generation
#' @param pop_pool
#' @param pareto
#' @param prob_matrix
#' @param current_time
#'
#' @return
#'
#'
#' @examples
#' @noRd
base_end_gen_callback <-
    function(userdir = "",
             generation,
             pop_pool,
             pareto,
             prob_matrix,
             current_time) {
        # print(Sys.time()- start_time)
    }

#' A default call_back function that does nothing.
#'
#' @param generation a number indicating the number of iterations of the galgo algorithm
#' @param pop_pool a \code{data.frame} with the solution vectors, number of clusters and their ranking.
#' @param pareto the solutions found by Galgo accross all generations in the solution space
#' @param prob_matrix a \code{matrix} or \code{data.frame}. Must be an expression matrix with features in rows and samples in columns
#' @param current_time an \code{POSIXct} object
#'
#' @return
#'
#' 
#' @examples
#' @noRd
default_callback <-
    function(userdir = "",
             generation,
             pop_pool,
             pareto,
             prob_matrix,
             current_time) {

    }

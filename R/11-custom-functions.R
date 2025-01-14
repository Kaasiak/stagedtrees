#' TO DO: Add examples
#'

#' Full search of an ordered variable.
#'
#' Full search on one level of a staged event tree with
#' iterative joining of "adjacent" stage for an ordered variable.
#' 
#' @param object an object of class \code{sevt} with fitted probabilities and 
#' data, as returned by \code{full} or \code{sevt_fit}.
#' @param variable name of a variables that should be considered for the optimization.
#' @param n_bins final number of stages.
#' @param n_init initial number of stages per subset considered for the optimization.
#' @param score the score function to be maximized.
#' @param ignore vector of stages which will be ignored and left untouched,
#'               by default the name of the unobserved stages stored in
#'               `object$name_unobserved`.
#' @details For the given variable the algorithm separates the stages into bins 
#' and joins all stages in one bin. It searches all possible splittings of the 
#' stages into the bins and returns the best scoring model.
#' @return The final staged event tree obtained.
#' @importFrom stats  BIC
#' @export
stages_full_ordered_search <-
  function(object,
           variable = NULL,
           n_init = NULL,
           n_bins = 3,
           score = function(x) {
             return(-BIC(x))
           },
           ignore = object$name_unobserved) {
    check_sevt_fit(object)
    v <- variable
    if (is.null(n_init)) {
      stop("Initial number of stages per subset must be specified")
    }
    if (is.null(variable)) {
      stop("The variable to be optimised must be specified")
    }
    if (length(object$stages[[v]]) %% n_init != 0) {
      stop("Starting number of stages must be a multiple of n_init")
    }
    n_subtrees <- length(object$stages[[v]]) / n_init
    partitions <- as.matrix(gtools::combinations(n = n_init - 1, r = n_bins - 1, repeats.allowed = FALSE))
    maxScore <- -Inf
    best_part <- NULL
    best_object <- NULL
    for (i in seq(1, nrow(partitions))) {  # try every partition
      try_object <- object
      part <- partitions[i,]
      part <- c(0, part, n_init)
      for (k in seq(1, n_subtrees)) {  # merge stages of every subtree
        stages <- as.character(seq((k - 1) * n_init + 1, k * n_init))
        try_object <- partition_stages(try_object, v, stages, part)
        try_score <- score(try_object)
      }
      if (try_score > maxScore) {
        maxScore <- try_score
        best_part <- part
        best_object <- try_object
      }
    }
    print(best_part)
    return(best_object)
}

#' Backward hill-climbing for ordered variables
#'
#' Greedy search on one level of a staged event tree with
#' iterative joining of "adjacent" stages for an ordererd variable.
#' 
#' @param object an object of class \code{sevt} with fitted probabilities and 
#' data, as returned by \code{full} or \code{sevt_fit}.
#' @param variable name of a variables that should be considered for the optimization.
#' @param n_init initial number of stages per subset considered for the optimization.
#' if FALSE in first step merges all stages associated with the same outcome into the same stage.
#' @param score the score function to be maximized.
#' @param max_iter the maximum number of iterations per variable.
#' @param ignore vector of stages which will be ignored and left untouched,
#'               by default the name of the unobserved stages stored in
#'               `object$name_unobserved`.
#' @details For the given variable the algorithm tries to join "adjacent" stages
#' and moves to the best model that increases the score. When no
#' increase is possible it terminates and returns the best scoring model.
#' @return The final staged event tree obtained.
#' @importFrom stats  BIC
#' @export
stages_ordered_bhc <-
  function(object,
           variable = NULL,
           n_init = NULL,
           score = function(x) {
             return(-BIC(x))
           },
           max_iter = Inf,
           ignore = object$name_unobserved) {
    check_sevt_fit(object)
    v <- variable
    if (is.null(n_init)) {
      stop("Initial number of stages per subset must be specified")
    }
    if (is.null(variable)) {
      stop("The variable to be optimised must be specified")
    }
    if (length(object$stages[[v]]) %% n_init != 0) {
      stop("Starting number of stages must be a multiple of n_init")
    }
    best_score <- score(object)
    n_subtrees <- length(object$stages[[v]]) / n_init
    done <- FALSE
    iter <- 0
    temp <- object # clone the object
    while (!done && iter < max_iter) {
      done <- TRUE
      iter <- iter + 1
      stages <- unique(temp$stages[[v]][1:n_init])
      n_stages <- length(stages)
      if (length(stages) > 1) {
        for (i in 2:length(stages)) {
          # scan the neighbouring stages
          s1 <- as.numeric(stages[i - 1])
          s2 <- as.numeric(stages[i])
          try_object <- temp
          for (k in seq(1, n_subtrees))  { # join the 2 stages in every subtree
            s1a <- s1 + n_init * (k - 1)
            s2a <- s2 + n_init * (k - 1)
            try_object <- join_stages(try_object, v, as.character(s1a), as.character(s2a))
          }
          try_score <- score(try_object)
          if (try_score >= best_score) {
            best_object <- try_object
            best_score <- score(best_object)
            done <- FALSE
          }
        }
      } ## end if there are more than 1 stage
      temp <- best_object
      best_score <- score(best_object)
    } ## end while
    return(best_object)
}

#' Join stages of a staged event tree according to a partition
#'
#' Join multiple stages in a staged event tree object, updating
#' probabilities and log-likelihood accordingly.
#'
#' @param object an object of class \code{sevt}.
#' @param variable variable.
#' @param stages vector of stages considered for joining.
#' @param part vector of partitioning the stages into n bins.
#' @return the staged event tree where the stages are joined according to 
#' the partition vector
#' @details This function joins all in one partition defined by the `part` vector,
#'          updating probabilities and log-likelihood if 
#'          the object was fitted.
#' @export
partition_stages <- function(object, variable, stages, part) {
  v <- variable
  for (i in seq(1, length(part) - 1)) {
    join_stages <- stages[(part[i] + 1):part[i + 1]]
    object <- join_multiple_stages(object, v, join_stages)
  }
  return(object)
}

#' Join multiple stages
#'
#' Join multiple stages in a staged event tree object, updating
#' probabilities and log-likelihood accordingly.
#'
#' @param object an object of class \code{sevt}.
#' @param variable name of a variables that should be considered for the optimization
#' @param join_stages vector of stages to be merged.
#' @return the staged event tree where all stages in `join_stages` are joined
#' @details This function joins all in one partition defined by the `part` vector,
#'          updating probabilities and log-likelihood if 
#'          the object was fitted.
#' @export
join_multiple_stages <- function(object, variable, join_stages) {
  v <- variable
  check_sevt(object)
  if (length(join_stages) < 2) {  # no stages to be merged
    return(object)
  }
  join_stages <- as.character(join_stages)
  s1 <- join_stages[1]
  s2 <- join_stages[2:length(join_stages)]
  k <- length(object$tree[[v]])
  st <- object$stages[[v]]

  if (!is.null(object$prob)) {
    probs <- as.matrix(expand_prob(object)[[v]])
    counts <- as.matrix(object$ctables[[v]])
    rownames(counts) <- object$stages[[v]]
    rownames(probs) <- object$stages[[v]]
    if (is.null(object$lambda)) {
      object$lambda <- 0
    }
    old_prob <- probs[join_stages,]
    old_ct <- counts[join_stages,]
    if (is.null(object$lambda)) {
      object$lambda <- 0
    }
    dll <- sum(old_ct[old_ct > 0] * log(old_prob[old_ct > 0]))
    
    object$prob[[v]][[s1]] <- colSums(old_ct) + object$lambda
    attr(object$prob[[v]][[s1]], "n") <- sum(old_ct)
    object$prob[[v]][[s1]] <-
      object$prob[[v]][[s1]] / sum(object$prob[[v]][[s1]])
    for (s in s2) {
      object$prob[[v]][[s]] <- NULL ## delete remaining stages
    }
    object$stages[[v]][st %in% s2] <- s1
    if (!is.null(object$ll)) {
      ## update log likelihood
      new_ct <- colSums(old_ct)
      object$ll <-
        object$ll - dll + sum(new_ct[new_ct > 0] *
                                log(object$prob[[v]][[s1]][new_ct > 0]))
      attr(object$ll, "df") <-
        attr(object$ll, "df") - (k - 1) * (length(join_stages) - 1)
    }
  }
  return(object)
}


#' Plotting of staged event tree for every step of the BHC algorithm
#' 
#' Displaying the plot of a staged event tree together with the 
#' probabilities respective to the `plot_var` variable.
#' 
#' @param object an object of class \code{sevt} with fitted probabilities and 
#' data, as returned by \code{full} or \code{sevt_fit}.
#' @param plot_var variable for plotting the barplots.
#' @param score the score function to be maximized.
#' @param max_iter the maximum number of iterations per variable.
#' @param scope names of variables that should be considered for the optimization.
#' @param ignore vector of stages which will be ignored and left untouched,
#'               by default the name of the unobserved stages stored in
#'               `object$name_unobserved`.
#' @param trace if >0 increasingly amount of info
#' is printed (via \code{message}).
#' @details The function simultaniously performs Backward Hill Climbing as
#' implemented in `stages_bhc()`. Additionally, at every merge of two stages the
#' current staged tree together with the probabilities corresponding to the 
#' `plot_var` are plotted.
#' @return The final staged event tree obtained and print a series of plots.
#' @importFrom stats  BIC
#' @importFrom graphics par
#' @export
stages_bhc_plot <-
  function(object,
           plot_var = NULL,
           score = function(x) {return(-BIC(x))},
           max_iter = Inf,
           scope = NULL,
           ignore = object$name_unobserved,
           trace = 0) {
    check_sevt_fit(object)
    now_score <- score(object)
    if (is.null(scope)) {
      scope <- sevt_varnames(object)[-1]
    }
    stopifnot(all(scope %in% sevt_varnames(object)[-1]))
    for (v in scope) {
      r <- 1
      iter <- 0
      done <- FALSE
      while (!done && iter < max_iter) {
        iter <- iter + 1
        temp <- object # clone the object
        temp_score <- now_score
        done <- TRUE
        stages <- unique(object$stages[[v]])
        stages <- stages[!(stages %in% ignore)]
        if (length(stages) > 1) {
          for (i in 2:length(stages)) {
            ## try all stages pair
            s1 <- stages[i]
            for (j in 1:(i - 1)) {
              s2 <- stages[j]
              try <-
                join_stages(object, v, s1, s2) ## join the 2 stages
              try_score <- score(try)
              if (try_score >= temp_score) {
                temp <- try
                temp_score <- try_score
                s1a <- s1
                s2a <- s2
                done <- FALSE
              }
            }
          }
        } ## end if there are more than 1 stage
        object <- temp
        now_score <- temp_score
        if (v == plot_var) {
          par(mfrow = c(1,2))
          plot(object)
          barplot(object, plot_var)
        }
        if ((trace > 1) && !done) {
          message(v, " joined stages: ", s1a, " and ", s2a)
        }
      } ## end while
      if (trace > 0) {
        message("BHC over ", v, " done after ", iter, " iterations")
      }
    } ## end for over variables
    if (trace > 0) {
      message("BHC done")
    }
    object$call <- sys.call()
    object$score <- list(value = now_score, f = score)
    par(mfrow = c(1,1))
    return(object)
  }


#' #' Backward hill-climbing for ordered variables
#' #'
#' #' Greedy search on one level of a staged event tree with
#' #' iterative joining of "adjacent" stages for an ordererd variable.
#' #'
#' #' @param object an object of class \code{sevt} with fitted probabilities and
#' #' data, as returned by \code{full} or \code{sevt_fit}.
#' #' @param variable name of a variables that should be considered for the optimization.
#' #' @param n_init initial number of stages per subset considered for the optimization.
#' #' @param per_subset if TRUE performs the merging of the stages in each subtree separately,
#' #' if FALSE in first step merges all stages associated with the same outcome into the same stage.
#' #' @param score the score function to be maximized.
#' #' @param max_iter the maximum number of iterations per variable.
#' #' @param ignore vector of stages which will be ignored and left untouched,
#' #'               by default the name of the unobserved stages stored in
#' #'               `object$name_unobserved`.
#' #' @param trace if >0 increasingly amount of info
#' #' is printed (via \code{message}).
#' #' @details For the given variable the algorithm tries to join "adjacent" stages
#' #' and moves to the best model that increases the score. When no
#' #' increase is possible it terminates and returns the best scoring model.
#' #' @return The final staged event tree obtained.
#' #' @importFrom stats  BIC
#' #' @export
#' stages_ordered_bhc <-
#'   function(object,
#'            variable = NULL,
#'            n_init = NULL,
#'            per_subset = FALSE,
#'            score = function(x) {
#'              return(-BIC(x))
#'            },
#'            max_iter = Inf,
#'            ignore = object$name_unobserved,
#'            trace = 0) {
#'     check_sevt_fit(object)
#'     if (is.null(n_init)) {
#'       stop("Initial number of stages per subset must be specified")
#'     }
#'     now_score <- score(object)
#'     if (is.null(variable)) {
#'       stop("The variable to be optimised must be specified")
#'     }
#'     v <- variable
#' 
#'     # initial merge of stages associated with the same outcomes
#'     if (!per_subset) {
#'       for (i in seq(1, n_init)) {
#'         join_stages <- seq(i, length(object$stages[[v]]), n_init)
#'         object <- join_multiple_stages(object, v, join_stages)
#'       }
#'     }
#' 
#'     stages <- unique(object$stages[[v]])
#'     stages <- stages[!(stages %in% ignore)]
#'     subset_stages_index <- seq(1,length(stages), n_init)  #starting indices for subsets of stages
#'     n_subsets <- length(subset_stages_index)
#'     start_stage <- 1
#'     for (k in seq(1, n_subsets)) {
#'       iter <- 0
#'       done <- FALSE
#'       n_stages <- n_init
#'       while (!done && iter < max_iter) {
#'         iter <- iter + 1
#'         temp <- object # clone the object
#'         temp_score <- now_score
#'         done <- TRUE
#'         stages <- unique(object$stages[[v]])[seq(start_stage, start_stage + n_stages - 1)]
#'         stages <- stages[!(stages %in% ignore)]
#'         if (length(stages) > 1) {
#'           for (i in 2:length(stages)) {
#'             ## scan the neighbouring stages
#'             s1 <- stages[i - 1]
#'             s2 <- stages[i]
#'             try <- join_stages(object, v, s1, s2) ## join the 2 stages
#'             try_score <- score(try)
#'             if (try_score >= temp_score) {
#'               temp <- try
#'               temp_score <- try_score
#'               s1a <- s1
#'               s2a <- s2
#'               done <- FALSE
#'             }
#'           }
#'         } ## end if there are more than 1 stage
#'         if (!done) {
#'           n_stages <- n_stages - 1
#'         }
#'         object <- temp
#'         now_score <- temp_score
#'         if ((trace > 1) && !done) {
#'           message(v, " joined stages: ", s1a, " and ", s2a)
#'         }
#'       } ## end while
#'       start_stage <- start_stage + n_stages
#'     }
#' 
#'     if (trace > 0) {
#'       message("BHC over ", v, " done after ", iter, " iterations")
#'     }
#'     if (trace > 0) {
#'       message("BHC done")
#'     }
#'     object$call <- sys.call()
#'     object$score <- list(value = now_score, f = score)
#'     return(object)
#' }

#' #' Full search of an ordered variable.
#' #'
#' #' Full search on one level of a staged event tree with
#' #' iterative joining of "adjacent" stage for an ordered variable.
#' #' 
#' #' @param object an object of class \code{sevt} with fitted probabilities and 
#' #' data, as returned by \code{full} or \code{sevt_fit}.
#' #' @param variable name of a variables that should be considered for the optimization.
#' #' @param n_bins final number of stages.
#' #' @param n_init initial number of stages per subset considered for the optimization.
#' #' @param per_subset if TRUE performs the merging of the stages in each subtree separately,
#' #' @param score the score function to be maximized.
#' #' @param ignore vector of stages which will be ignored and left untouched,
#' #'               by default the name of the unobserved stages stored in
#' #'               `object$name_unobserved`.
#' #' @details For the given variable the algorithm separates the stages into bins 
#' #' and joins all stages in one bin. It searches all possible splittings of the 
#' #' stages into the bins and returns the best scoring model.
#' #' @return The final staged event tree obtained.
#' #' @importFrom stats  BIC
#' #' @export
#' full_ordered_search <- function(
#'                               object,
#'                               n_bins = 3,
#'                               n_init = NULL,
#'                               variable = NULL,
#'                               per_subset = FALSE, 
#'                               ignore = object$name_unobserved,
#'                               score = function(x) {return(logLik(x))}) {
#'   v <- variable
#'   
#'   if (is.null(v)) {
#'     stop("Variable to be optimised must be specified")
#'   }
#'   if (is.null(n_init)) {
#'     stop("Initial number of stages must be specified")
#'   }
#'   # initial merge of stages associated with the same outcomes
#'   if (!per_subset) {
#'     for (i in seq(1, n_init)) {
#'       join_stages <- seq(i, length(object$stages[[v]]), n_init)
#'       object <- join_multiple_stages(object, v, join_stages)
#'     }
#'   }
#'   stages <- unique(object$stages[[v]])
#'   subset_stages_index <- seq(1,length(stages), n_init)  #starting indices for subsets of stages
#'   n_subsets <- length(subset_stages_index)
#'   start_stage <- 1
#'   for (k in seq(1, n_subsets)) {
#'     subset_stages <- stages[seq(start_stage, start_stage + n_init - 1)]
#'     partitions <- as.matrix(gtools::combinations(n = n_init - 1, r = n_bins - 1, repeats.allowed = FALSE))
#'     maxScore <- -Inf
#'     best_part <- NULL
#'     best_object <- NULL
#'     for (i in seq(1, nrow(partitions))) {
#'       part <- partitions[i,]
#'       part <- c(0, part, n_init)
#'       try_object <- partition_stages(object, v, subset_stages, part)
#'       try_score <- score(try_object)
#'       if (try_score > maxScore) {
#'           maxScore <- try_score
#'           best_part <- part
#'           best_object <- try_object
#'       }
#'     }
#'     start_stage <- start_stage + n_init
#'     object <- best_object
#'     print(best_part)
#'   }
#'   return(best_object)
#' }

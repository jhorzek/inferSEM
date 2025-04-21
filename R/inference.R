#' get_implied
#'
#' Computes the model implied means and covariances. In contrast to mxGetExpected,
#' this also includes the means and covariances for the latent variables
#'
#' @param A matrix with directed effects
#' @param S matrix with undirected effects
#' @param M matrix with intercepts
#' @returns list with implied means and covariances
#' @noRd
#' @examples
#' library(mxsem)
#'
#' model <- '
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + a1*y2 + b*y3 + c1*y4
#'      dem65 =~ y5 + a2*y6 + b*y7 + c2*y8
#'
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#'
#'   # residual correlations
#'     y1 ~~ y5
#'     y2 ~~ y4 + y6
#'     y3 ~~ y7
#'     y4 ~~ y8
#'     y6 ~~ y8
#' '
#'
#' fit <- mxsem(model = model,
#'              data  = OpenMx::Bollen) |>
#'   mxTryHard()
#' inferSEM:::get_implied(A = fit$A$values,
#'                        S = fit$S$values,
#'                        M = fit$M$values)
get_implied <- function(A, S, M){
  I <- diag(nrow(A))
  implied_means <- t(solve(I - A) %*% t(M))
  implied_covariances <- solve(I - A) %*% S %*% t(solve(I - A))

  return(list(
    implied_means = implied_means,
    implied_covariances = implied_covariances
  ))
}

#' check_specification
#'
#' Checks that the setup was correct.
#'
#' @param model OpenMx RAM model
#' @param variable_names names of the variables in the model
#' @param input_vector interventional or conditional vector
#' @returns nothing. Thows error in case of misspecification
#' @noRd
check_specification <- function(model, variable_names, input_vector){
  if(is.null(names(input_vector)))
    stop("input_vector must be a vector with names.")
  if(any(!names(input_vector) %in% variable_names))
    stop("The following variables were not found in the model: ",
         paste0(names(input_vector)[!names(input_vector) %in% variable_names],
                collapse = ", "))
}

#' infer
#'
#' Computes the interventional and / or conditional distribution of a SEM based
#' on an intervention and / or an observation of variable values.
#'
#' For a thorough introduction to interventional and conditional distributions in
#' SEM, please see Gische et al. (2021).
#'
#'
#' References:
#'
#' Gische, C., West, S. G., & Voelkle, M. C. (2021). Forecasting causal effects
#' of interventions versus predicting future outcomes. Structural Equation Modeling:
#' A Multidisciplinary Journal, 28(3), 475-492.
#'
#' @param model OpenMX RAM model. The model must have matrices A, S, and M.
#' @param intervene named vector with variables on which is intervened and the
#' respective intervention levels
#' @param observe named vector with variables for which a specific value is observed
#' (conditioning)
#' @returns list with expected means and covariances
#' @import mxsem
#' @import OpenMx
#' @importFrom condMVNorm condMVN
#' @export
#' @examples
#' library(mxsem)
#'
#' model <- '
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + a1*y2 + b*y3 + c1*y4
#'      dem65 =~ y5 + a2*y6 + b*y7 + c2*y8
#'
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#'
#'   # residual correlations
#'     y1 ~~ y5
#'     y2 ~~ y4 + y6
#'     y3 ~~ y7
#'     y4 ~~ y8
#'     y6 ~~ y8
#' '
#'
#' fit <- mxsem(model = model,
#'              data  = OpenMx::Bollen) |>
#'   mxTryHard()
#'
#' # intervention on dem60
#' infer(model = fit,
#'       intervene = c("dem60" = 2))
#' # condition on dem60
#' infer(model = fit,
#'       observe = c("dem60" = 2))
#' # condition and intervene
#' infer(model = fit,
#'       intervene = c("ind60" = -1),
#'       observe = c("dem60" = 2))
infer <- function(model,
                  intervene = NULL,
                  observe = NULL){

  A <- model$A$values
  S <- model$S$values
  M <- model$M$values

  variable_names <- colnames(A)

  if(!is.null(intervene))
    check_specification(model = model,
                        variable_names = variable_names,
                        input_vector = intervene)
  if(!is.null(observe))
    check_specification(model = model,
                        variable_names = variable_names,
                        input_vector = observe)

  if(!is.null(intervene)){
    for(i in names(intervene)){
      # Set the expected mean to the intervention level
      M[, i] <- intervene[i]
      # Set all effects on the intervention variable to 0
      A[i, ] <- 0
      # Set all variances and covariances with the intervention variable to 0
      S[i, ] <- 0
      S[, i] <- 0
      # If we want to compute the condintional distrubtion, we have to make
      # the matrix positive definite
      if(!is.null(observe))
        S[i, i] <- 1e-6
    }
  }

  implied <- get_implied(A = A, S = S, M = M)

  if(!is.null(observe)){
    conditional <- condMVNorm::condMVN(mean = implied$implied_means,
                                       sigma = implied$implied_covariances,
                                       dependent.ind = which(!variable_names %in% names(observe)),
                                       given.ind = which(variable_names %in% names(observe)),
                                       X.given = observe)
    implied$implied_means <- matrix(conditional$condMean, nrow = 1)
    colnames(implied$implied_means) <- colnames(conditional$condVar)
    implied$implied_covariances <- conditional$condVar
  }

  return(
    list(
      means = implied$implied_means,
      covariances = implied$implied_covariances
    )
  )
}

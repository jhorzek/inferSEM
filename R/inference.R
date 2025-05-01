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
#' @param variable_names names of the variables in the model
#' @param input_vector interventional or conditional vector
#' @returns nothing. Thows error in case of misspecification
#' @noRd
check_specification <- function(variable_names, input_vector){
  if(is.null(names(input_vector)))
    stop("input_vector must be a vector with names.")
  if(any(!names(input_vector) %in% variable_names))
    stop("The following variables were not found in the model: ",
         paste0(names(input_vector)[!names(input_vector) %in% variable_names],
                collapse = ", "))
}

#' get_ram_matrices_openmx
#'
#' Extract RAM matrices from an OpenMx model
#' @param model OpenMx model
#' @returns list with A, S, and M matrix
#' @noRd
get_ram_matrices_openmx <- function(model){
  return(
    list(
      A = model[[model$expectation@A]]$values,
      S = model[[model$expectation@S]]$values,
      M = model[[model$expectation@M]]$values
    )
  )
}

#' get_ram_matrices_lavaan
#'
#' Extract RAM matrices from a lavaan model
#' @param model lavaan model
#' @returns list with A, S, and M matrix
#' @noRd
#' @importFrom lavaan lavMatrixRepresentation
#' @importFrom lavaan parTable
get_ram_matrices_lavaan <- function(model){
  ram_parameter_table <- lavaan::lavMatrixRepresentation(lavaan::parTable(model),
                                                         representation = "RAM")
  ram_parameter_table$mat <- toupper(ram_parameter_table$mat)
  if(!"M" %in% ram_parameter_table$mat)
    stop("Your lavaan model must be estimated with intercepts. Please use meanstructure = TRUE.")

  # Dimensions of the matrices:
  observed <- model@Model@dimNames[[1]][[1]]
  if(length(model@Model@dimNames[[1]]) > 1){
    latent <- model@Model@dimNames[[1]][[2]]
  }else{
    latent <- c()
  }
  variable_names <- c(observed, latent)
  n_variables <- length(variable_names)

  ram_matrices <- list(
    M = matrix(data = 0, nrow = n_variables, ncol = 1, dimnames = list(variable_names, NULL)),
    A = matrix(data = 0, nrow = n_variables, ncol = n_variables, dimnames = list(variable_names, variable_names)),
    S = matrix(data = 0, nrow = n_variables, ncol = n_variables, dimnames = list(variable_names, variable_names))
  )

  for(i in 1:nrow(ram_parameter_table)){
    if(ram_parameter_table$mat[i] == "")
      next
    ram_matrices[[ram_parameter_table$mat[i]]][ram_parameter_table$row[i], ram_parameter_table$col[i]] <- ram_parameter_table$est[i]
    if(ram_parameter_table$mat[i] == "S")
      ram_matrices[[ram_parameter_table$mat[i]]][ram_parameter_table$col[i], ram_parameter_table$row[i]] <- ram_parameter_table$est[i]
  }

  ram_matrices$M <- t(ram_matrices$M)
  return(ram_matrices)
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
#' @param model OpenMX RAM model or lavaan model
#' @param intervene named vector with variables on which is intervened and the
#' respective intervention levels
#' @param observe named vector with variables for which a specific value is observed
#' (conditioning)
#' @returns list with expected means and covariances
#' @import mxsem
#' @importFrom methods is
#' @importFrom condMVNorm condMVN
#' @export
#' @examples
#' library(mxsem)
#'
#' model <- '
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'      dem65 =~ y5 + a*y6 + b*y7 + c*y8
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

  if(is(model, "MxRAMModel")){
    ram <- get_ram_matrices_openmx(model)
    A <- ram$A
    S <- ram$S
    M <- ram$M
  }else if(is(model, "lavaan")){
    ram <- get_ram_matrices_lavaan(model)
    A <- ram$A
    S <- ram$S
    M <- ram$M
  }

  variable_names <- colnames(A)

  if(!is.null(intervene))
    check_specification(variable_names = variable_names,
                        input_vector = intervene)
  if(!is.null(observe))
    check_specification(variable_names = variable_names,
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
    }
  }

  implied <- get_implied(A = A, S = S, M = M)

  if(any(eigen(implied$implied_covariances, only.values = TRUE)$values <= 0)){
    warning("The implied covariance matrix is not positive definite. The matrix",
            " will be made positive definite by adding a small constant to the diagonal.")
    if(any(diag(implied$implied_covariances == 0))){
      diag(implied$implied_covariances)[diag(implied$implied_covariances) == 0] <- 1e-6
    }
    for(i in 1:100){
      if(!any(eigen(implied$implied_covariances, only.values = TRUE)$values <= 0))
        break
      diag(implied$implied_covariances) <- diag(implied$implied_covariances) + 1e-6
    }
  }

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

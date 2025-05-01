test_that("Comparing lavaan and OpenMx", {
  library(mxsem)
  library(lavaan)
  library(inferSEM)
  set.seed(2342)

  model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a1*y2 + b*y3 + c1*y4
     dem65 =~ y5 + a2*y6 + b*y7 + c2*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

  fit_mx <- mxsem(model = model,
                  data  = OpenMx::Bollen) |>
    mxTryHard(exhaustive = TRUE)

  fit_lavaan <- sem(model = model,
                    data = OpenMx::Bollen,
                    meanstructure = TRUE)

  infer_mx <- inferSEM::infer(model = fit_mx, intervene = c("dem60" = 2))
  infer_lavaan <- inferSEM::infer(model = fit_lavaan, intervene = c("dem60" = 2))

  testthat::expect_lt(max(abs(infer_mx$means -
                                  infer_lavaan$means[,colnames(infer_mx$means)])), .001)
  testthat::expect_lt(max(abs(infer_mx$covariances -
                                infer_lavaan$covariances[rownames(infer_mx$covariances),
                                                         colnames(infer_mx$covariances)])), .001)
})

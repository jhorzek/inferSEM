test_that("Test computation of distributions in Gische 2021", {
  # The following code is copy-pasted from:
  # Gische, C., West, S. G., & Voelkle, M. C. (2021). Forecasting causal
  # effects of interventions versus predicting future outcomes.
  # Structural Equation Modeling: A Multidisciplinary Journal, 28(3), 475-492.
  library(matrixcalc)
  library(condMVNorm)
  library(inferSEM)

  library(matrixcalc)
  #################
  # read in parameter values for the data generation
  # and obtain matrices needed for data generation
  theta_dgp<-c(0.05,0.4,-0.6,1.2,20,40)
  Sigma_epsilon<-matrix(nrow=2,c(theta_dgp[5],0,0,theta_dgp[6]))
  B<-matrix(nrow=2,theta_dgp[c(1,2,3,4)],byrow = T)
  # check if process satisfies regularity conditions
  abs(eigen(B)$values)
  # compute the value of the auto-covariance function at lag zero
  # assuming that the process has an infinite past
  I_4<-diag(1, nrow = (4))
  Gamma_0<-solve(I_4-B%x%B)%*%vec(Sigma_epsilon)
  Gamma_0<-matrix(Gamma_0, nrow=2)
  # add parameters for the initial values that correspond to the
  # long term stable values of the system
  theta_dgp<-round(c(B[1,1],B[1,2],B[2,1],B[2,2],Gamma_0[1,1],Gamma_0[2,2],Gamma_0[1,2],Sigma_epsilon[1,1],Sigma_epsilon[2,2]),2)
  theta_dgp
  # add variances and covariances of the person-specific
  # time invariant variables; since Ito et al. (1998) contains only a single time series
  # we postulate a hypothetical between-person variance:
  # exx=0.25*varepsilonxx, eyy=0.25*varepsilonyy
  Sigma_eta<-matrix(nrow=2, c(5,2.5,2.5,10))
  is.positive.definite(Sigma_eta)
  eigen(Sigma_eta)
  correxey<-Sigma_eta[1,2]/(sqrt(Sigma_eta[1,1]*Sigma_eta[2,2]))
  # add loadings of the person-specific time invariant variables
  # onto X and Y; the loading on the initial variable is set
  # equal to the long term effect
  I_2<-diag(1, nrow=2)
  tao<-solve(I_2-B)
  # the long run cumulative effect of a unit impulse on x
  # is -4; that of a unit impulse on y is 19;
  theta_dgp<-round(c(tao[1,1],tao[2,2],tao[1,2],tao[2,1],B[1,1],B[1,2],B[2,1],
                     B[2,2],Sigma_eta[1,1],Sigma_eta[2,2],Sigma_eta[1,2],
                     Gamma_0[1,1],Gamma_0[2,2],Gamma_0[1,2],Sigma_epsilon[1,1],Sigma_epsilon[2,2]),2)
  # create data frame for parameters
  thetas<-data.frame(t(1:16))
  names(thetas)<-c("cx0ex","cy0ey","cx0ey","cy0ex",
                   "cxx","cxy","cyx","cyy","psiexex","psieyey",
                   "psiexey","psixx","psiyy","psixy","psixr","psiyr")
  thetas[2,]<-theta_dgp

  #####################################################################
  # Part II: Results for a homogeneous population
  #####################################################################
  # Define the matrix of structural coefficients,
  # the error covariance matrix and
  # compute the model implied moments
  k<-2 # number of processes (insulin, glucose)
  T<-4 # number of time points / measurement occasions
  n<-k*T # number of variables in the system
  r<-2 # row number from which the theta values are exctracted
  Psi_hom<-matrix(nrow=n,ncol=n,0)
  Psi_hom[1,1]<-thetas$psixx[r]
  Psi_hom[3,3]<-thetas$psixr[r]
  Psi_hom[5,5]<-thetas$psixr[r]
  Psi_hom[7,7]<-thetas$psixr[r]
  Psi_hom[2,2]<-thetas$psiyy[r]
  Psi_hom[4,4]<-thetas$psiyr[r]
  Psi_hom[6,6]<-thetas$psiyr[r]
  Psi_hom[8,8]<-thetas$psiyr[r]
  Psi_hom[1,2]<-thetas$psixy[r]
  Psi_hom[2,1]<-thetas$psixy[r]
  Psi_hom
  C_hom<-matrix(nrow=n,ncol=n,0)
  C_hom[3,1]<-thetas$cxx[r]
  C_hom[3,2]<-thetas$cxy[r]
  C_hom[4,1]<-thetas$cyx[r]
  C_hom[4,2]<-thetas$cyy[r]
  C_hom[5,3]<-thetas$cxx[r]
  C_hom[5,4]<-thetas$cxy[r]
  C_hom[6,3]<-thetas$cyx[r]
  C_hom[6,4]<-thetas$cyy[r]
  C_hom[7,5]<-thetas$cxx[r]
  C_hom[7,6]<-thetas$cxy[r]
  C_hom[8,5]<-thetas$cyx[r]
  C_hom[8,6]<-thetas$cyy[r]
  C_hom
  # expected values of error terms are zero
  Eepsilon = rep(0,n)
  # calculating the model implied covariance matrix
  # and model implied mean vector of the joint distribution
  I_Cinv_hom<-solve(diag(1,nrow=n)-C_hom)
  Sigma_V_hom = I_Cinv_hom%*%Psi_hom%*%t(I_Cinv_hom)
  E_V_hom = I_Cinv_hom%*%Eepsilon
  # define the interventional values / values of the conditional
  # variables one population sd of V3
  x2<-round(sqrt(Sigma_V_hom[3,3]),2)
  # Calculate the interventional mean and covariance matrix
  # for do(x2)
  e83<-c(0,0,1,0,0,0,0,0)
  IN83<-diag(c(1,1,0,1,1,1,1,1))
  T1_hom = solve(diag(1,nrow=n)-IN83%*%C_hom)%*%IN83
  a1_hom = solve(diag(1,nrow=n)-IN83%*%C_hom)%*%e83
  Edox2_hom=a1_hom*x2

  Vdox2_hom = T1_hom%*%Psi_hom%*%t(T1_hom)

  # We want to reproduce this distribution with inferSEM. To this end,
  # we have to first set up the model in OpenMx:
  manifest_names <- paste0(rep(c("x", "y"), 4), rep(1:4, each = 2))
  dimnames(C_hom) <- list(manifest_names,
                          manifest_names)
  dimnames(Psi_hom) <- list(manifest_names,
                            manifest_names)
  F <- diag(8)
  dimnames(F) <- list(manifest_names,
                      manifest_names)
  M <- matrix(0, nrow = 1, ncol = 8)
  dimnames(M) <- list(NULL, manifest_names)
  model <- OpenMx::mxModel(
    type = "RAM",
    manifestVars = manifest_names,
    OpenMx::mxMatrix(type = "Full", values = C_hom, free = FALSE, name = "A"),
    OpenMx::mxMatrix(type = "Full", values = Psi_hom, free = FALSE, name = "S"),
    OpenMx::mxMatrix(type = "Full", values = F, free = FALSE, name = "F"),
    OpenMx::mxMatrix(type = "Full", values = M, free = FALSE, name = "M"),
    OpenMx::mxExpectationRAM(A = "A", S = "S", F = "F", M = "M"),
    OpenMx::mxFitFunctionML()
  )

  # Next, we compute the interventional distribution given x2:
  inferred <- inferSEM::infer(model = model,
                              intervene = c("x2" = x2))
  testthat::expect_true(all(abs(t(Edox2_hom) - inferred$means) < 1e-6))
  testthat::expect_true(all(abs(Vdox2_hom - inferred$covariances) < 1e-6))

  # Next, we compute the conditional distributions
  # Calculate the conditional mean and covariance matrix
  # for X2=x2
  # reorder mean vector and covariance matrix from
  # 1,2,3,4,5,6 to 1,2,4,5,6,3 and obtain partitions
  E_V_hom_reorder = E_V_hom[c(1,2,4,5,6,7,8,3)]
  E_V_hom_reorder1 = E_V_hom_reorder[c(1,2,3,4,5,7,8)]
  E_V_hom_reorder2 = E_V_hom_reorder[8]
  Sigma_V_hom_reorder = Sigma_V_hom[, c(1,2,4,5,6,7,8,3)]
  Sigma_V_hom_reorder = Sigma_V_hom_reorder[c(1,2,4,5,6,7,8,3), ]
  # define selection matrices to obtain
  # partitions of the reordered covariance matrix
  e88 = c(0,0,0,0,0,0,0,1)
  IS8o8 = matrix(c(1,0,0,0,0,0,0,0,
                   0,1,0,0,0,0,0,0,
                   0,0,1,0,0,0,0,0,
                   0,0,0,1,0,0,0,0,
                   0,0,0,0,1,0,0,0,
                   0,0,0,0,0,1,0,0,
                   0,0,0,0,0,0,1,0), ncol=7)
  Sigma_V_hom_reorder22 = e88%*%Sigma_V_hom_reorder%*%e88
  Sigma_V_hom_reorder11 = t(IS8o8)%*%Sigma_V_hom_reorder%*%IS8o8
  Sigma_V_hom_reorder21 = t(e88)%*%Sigma_V_hom_reorder%*%IS8o8
  Sigma_V_hom_reorder12 = t(IS8o8)%*%Sigma_V_hom_reorder%*%e88
  Ex2_hom = E_V_hom_reorder1 + Sigma_V_hom_reorder12/as.numeric(Sigma_V_hom_reorder22)*(x2 - E_V_hom_reorder2)
  Vx2_hom = Sigma_V_hom_reorder11 - Sigma_V_hom_reorder12%*%Sigma_V_hom_reorder21/as.numeric(Sigma_V_hom_reorder22)

  # Again, we compare this to inferSEM
  inferred <- inferSEM::infer(model = model,
                              observe = c("x2" = x2))
  testthat::expect_true(all(abs(t(Ex2_hom) - inferred$means) < 1e-6))
  testthat::expect_true(all(abs(Vx2_hom - inferred$covariances) < 1e-6))

  #####################################################################
  # Part III: Results for a heterogeneous population
  #####################################################################
  # Define the matrix of structural coefficients,
  # the error covariance matrix and compute the model implied moments
  l<-2 # number of random intercepts
  n_het<-k*T+l # number of variables in the system
  r<-2 # row number from which the theta values are exctracted
  Psi_het<-matrix(nrow=n_het,ncol=n_het,0)
  Psi_het[1,1]<-thetas$psiexex[r]
  Psi_het[2,2]<-thetas$psieyey[r]
  Psi_het[1,2]<-thetas$psiexey[r]
  Psi_het[2,1]<-thetas$psiexey[r]
  Psi_het[3,3]<-thetas$psixx[r]
  Psi_het[4,4]<-thetas$psiyy[r]
  Psi_het[3,4]<-thetas$psixy[r]
  Psi_het[4,3]<-thetas$psixy[r]
  Psi_het[5,5]<-thetas$psixr[r]
  Psi_het[7,7]<-thetas$psixr[r]
  Psi_het[9,9]<-thetas$psixr[r]
  Psi_het[6,6]<-thetas$psiyr[r]
  Psi_het[8,8]<-thetas$psiyr[r]
  Psi_het[10,10]<-thetas$psiyr[r]
  C_het<-matrix(nrow=n_het,ncol=n_het,0)
  C_het[3,1]<-thetas$cx0ex[r]
  C_het[3,2]<-thetas$cx0ey[r]
  C_het[4,1]<-thetas$cy0ex[r]
  C_het[4,2]<-thetas$cy0ey[r]
  C_het[5,1]<-1
  C_het[6,2]<-1
  C_het[5,3]<-thetas$cxx[r]
  C_het[5,4]<-thetas$cxy[r]
  C_het[6,3]<-thetas$cyx[r]
  C_het[6,4]<-thetas$cyy[r]
  C_het[7,1]<-1
  C_het[8,2]<-1
  C_het[7,5]<-thetas$cxx[r]
  C_het[7,6]<-thetas$cxy[r]
  C_het[8,5]<-thetas$cyx[r]
  C_het[8,6]<-thetas$cyy[r]
  C_het[9,1]<-1
  C_het[10,2]<-1
  C_het[9,7]<-thetas$cxx[r]
  C_het[9,8]<-thetas$cxy[r]
  C_het[10,7]<-thetas$cyx[r]
  C_het[10,8]<-thetas$cyy[r]
  # expected values of error terms are zero
  Eepsilon = rep(0,n_het)
  ######################################################
  # calculating the model implied covariance matrix
  # and model implied mean vector of the joint distribution
  I_Cinv_het<-solve(diag(1,nrow=n_het)-C_het)
  Sigma_V_het = I_Cinv_het%*%Psi_het%*%t(I_Cinv_het)
  E_V_het = I_Cinv_het%*%Eepsilon
  # define the interventional values / values of the conditional
  # variables one population sd of V5
  x2<-round(sqrt(Sigma_V_het[5,5]),2)
  # define the values of the person-specific characteristics
  # for threee prototypical individuals from the population
  z_Amy<-c(-round(sqrt(diag(Sigma_V_het))[1],2),-round(sqrt(diag(Sigma_V_het))[2],2))
  z_Sam<-c(round(sqrt(diag(Sigma_V_het))[1],2),round(sqrt(diag(Sigma_V_het))[2],2))
  z_Joe<-c(0,0)
  # in the following we calculate the person-specific values for Amy
  z<-z_Amy
  #####################################################################################################

  # Calculate the interventional mean and covariance matrix
  # for do(X2=x2)
  e105<-c(0,0,0,0,1,0,0,0,0,0)
  IN105<-diag(c(1,1,1,1,0,1,1,1,1,1))
  T1_het = solve(diag(1,nrow=n_het)-IN105%*%C_het)%*%IN105
  a1_het = solve(diag(1,nrow=n_het)-IN105%*%C_het)%*%e105
  Edox2_het=a1_het*x2
  Vdox2_het = T1_het%*%Psi_het%*%t(T1_het)
  # get the interventional distribution of all
  # non-interventional variables
  one10_o5<-matrix(ncol=(n_het-1), c( 1,0,0,0,0,0,0,0,0,0,
                                      0,1,0,0,0,0,0,0,0,0,
                                      0,0,1,0,0,0,0,0,0,0,
                                      0,0,0,1,0,0,0,0,0,0,
                                      0,0,0,0,0,1,0,0,0,0,
                                      0,0,0,0,0,0,1,0,0,0,
                                      0,0,0,0,0,0,0,1,0,0,
                                      0,0,0,0,0,0,0,0,1,0,
                                      0,0,0,0,0,0,0,0,0,1))
  # via selection from the interventional moments
  Edox2_non_het = t(one10_o5)%*%Edox2_het
  Vdox2_non_het = t(one10_o5)%*%Vdox2_het%*%one10_o5

  # Comparison to inferSEM
  latent_names <- c("eta_x", "eta_y")
  manifest_names <- c(paste0(rep(c("x", "y"), 4), rep(1:4, each = 2)))
  variable_names <- c(latent_names, manifest_names)
  dimnames(C_het) <- list(variable_names,
                          variable_names)
  dimnames(Psi_het) <- list(variable_names,
                            variable_names)
  F <- cbind(diag(length(manifest_names)), matrix(0, nrow = length(manifest_names), ncol = length(latent_names)))
  dimnames(F) <- list(manifest_names,
                      c(manifest_names, latent_names))
  M <- matrix(0, nrow = 1, ncol = nrow(C_het))
  dimnames(M) <- list(NULL, variable_names)
  model <- OpenMx::mxModel(
    type = "RAM",
    manifestVars = manifest_names,
    latentVars = latent_names,
    OpenMx::mxMatrix(type = "Full", values = C_het[colnames(F), colnames(F)], free = FALSE, name = "A"),
    OpenMx::mxMatrix(type = "Full", values = Psi_het[colnames(F), colnames(F)], free = FALSE, name = "S"),
    OpenMx::mxMatrix(type = "Full", values = F, free = FALSE, name = "F"),
    OpenMx::mxMatrix(type = "Full", values = M[1, colnames(F), drop = FALSE], free = FALSE, name = "M"),
    OpenMx::mxExpectationRAM(A = "A", S = "S", F = "F", M = "M"),
    OpenMx::mxFitFunctionML()
  )
  inferred <- inferSEM::infer(model = model,
                              intervene = c("x2" = x2))
  testthat::expect_true(all(abs(Edox2_non_het[,1] - inferred$means[,variable_names[!variable_names == "x2"]]) < 1e-6))
  testthat::expect_true(all(abs(Vdox2_non_het - inferred$covariances[variable_names[!variable_names == "x2"],
                                                                     variable_names[!variable_names == "x2"]]) < 1e-6))

  # Distribution for Amy
  # Calculate the conditional mean and covariance matrix
  # given etax=zx and etay=zy of the interventional distribution etax=zx and etay=zy
  Edox2c12_re = Edox2_non_het[c(3,4,5,6,7,8,9,1,2)]
  # obtain partitions of the reordered mean vector
  Edox2c12_1 = Edox2c12_re[1:7]
  Edox2c12_2 = Edox2c12_re[c(8,9)]
  # change the oder of variabels in the interventional covariance
  # matrix from 1,2,,4,5,6 to 1,2,5,6,4
  Vdox2c12_re = Vdox2_non_het[, c(3,4,5,6,7,8,9,1,2)]
  Vdox2c12_re = Vdox2c12_re[c(3,4,5,6,7,8,9,1,2), ]
  # define selection matrices to obtain
  # partitions of the reordered covariance matrix
  I9c12_1 = matrix(c(1, 0, 0, 0, 0,0,0, 0, 0,
                     0, 1, 0, 0, 0,0,0, 0, 0,
                     0, 0, 1, 0, 0,0,0, 0, 0,
                     0, 0, 0, 1, 0,0,0, 0, 0,
                     0, 0, 0, 0, 1,0,0, 0, 0,
                     0, 0, 0, 0, 0,1,0, 0, 0,
                     0, 0, 0, 0, 0,0,1, 0, 0), ncol=7)
  I9c12_2 = matrix(c(0,0,0,0,0,0,0,1,0,
                     0,0,0,0,0,0,0,0,1),ncol=2)
  Vdox2c12_22 = t(I9c12_2)%*%Vdox2c12_re%*%I9c12_2
  Vdox2c12_11 = t(I9c12_1)%*%Vdox2c12_re%*%I9c12_1
  Vdox2c12_21 = t(I9c12_2)%*%Vdox2c12_re%*%I9c12_1
  Vdox2c12_12 = t(I9c12_1)%*%Vdox2c12_re%*%I9c12_2
  Edox2c12 = Edox2c12_1 + Vdox2c12_12%*%solve(Vdox2c12_22)%*%(z - Edox2c12_2)
  Vdox2c12 = Vdox2c12_11 - Vdox2c12_12%*%solve(Vdox2c12_22)%*%Vdox2c12_21

  # Same thing with inferSEM
  inferred <- infer(model = model,
                    intervene = c("x2" = x2),
                    observe = c(eta_x = z[1], eta_y = z[2]))

  testthat::expect_true(all(abs(Edox2c12[,1] - inferred$means[,colnames(inferred$means) != "x2"]) < 1e-3))
  testthat::expect_true(all(abs(Vdox2c12 - inferred$covariances[rownames(inferred$covariances) != "x2",
                                                                colnames(inferred$covariances) != "x2"]) < 1e-3))
})

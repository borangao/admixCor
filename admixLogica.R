library(data.table)
library(snpStats)
library(mvtnorm)
library(dplyr)
library(peakRAM)
library(Matrix)
library(CompQuadForm)
library(Rcpp)
sourceCpp("./admixLogica.cpp")
run_admix_analysis <- function( X1ty, X2ty, X1tX1, X2tX2, X1tX2, N, max_iter = 100, tol = 1e-6) {
  
  
  ld_score1 <- apply(X1tX1 / N, 2, function(x) sum(x^2))
  ld_score2 <- apply(X2tX2 / N, 2, function(x) sum(x^2))
  ld_score3 <- apply(X1tX2 / N, 2, function(x) sum(x^2))
  N_SNP <- ncol(X1tX1)
  sigma_g1_initial_MoM <- mean((X1ty/sqrt(N))^2 - 1) / mean(ld_score1) / N_SNP
  sigma_g2_initial_MoM <- mean((X2ty/sqrt(N))^2 - 1) / mean(ld_score2) / N_SNP

  # --- Ancestry 1 ---
  # 1. eigen‐decompose X1ᵀX1
  eig1   <- eigen(X1tX1, symmetric=TRUE)
  eigval1<- eig1$values
  U1     <- eig1$vectors
  # 2. rotate X1ty
  ut1    <- as.vector(t(U1) %*% X1ty)
  
  # 3. optimize h2 in (0,1)
  res1 <- optimize(
    f        = function(h2) -llk_uni_h2(h2, eigval1, ut1, N),
    interval = c(0.001/N_SNP, 0.1/N_SNP),
    tol      = 1e-8
  )
  sigma_g1_initial <- res1$minimum
  
  # --- Ancestry 2 ---
  eig2    <- eigen(X2tX2, symmetric=TRUE)
  eigval2 <- eig2$values
  ut2     <- as.vector(t(eig2$vectors) %*% X2ty)
  
  res2 <- optimize(
    f        = function(h2) -llk_uni_h2(h2, eigval2, ut2, N),
    interval = c(0.001/N_SNP, 0.1/N_SNP),
    tol      = 1e-8
  )
  sigma_g2_initial <- res2$minimum
 # admixres <- admixLogica(X1ty, X2ty, X1tX1, X2tX2, X1tX2, N, max_iter=100,
#                       sigma_g1_initial, sigma_g2_initial,
 #                      rho_g_initial, sigma_e_sq = 1, tol = tol)
  rho_g_initial=0
  admixres <-admixLogica_PXEM(X1ty, X2ty, X1tX1, X2tX2, X1tX2, N, max_iter=300,
                   sigma_g1_initial, sigma_g2_initial,
                   rho_g_initial, sigma_e_sq = 1 -N_SNP*(sigma_g2_initial+sigma_g1_initial) , tol = tol)
  
  Q1 <- as.numeric(crossprod(X1ty))
  Q2 <- as.numeric(crossprod(X2ty))
  Q12 <- as.numeric(t(X1ty) %*% X2ty)
  
  lambda1 <- eigen(X1tX1)$values
  lambda1<-lambda1[lambda1 > mean(lambda1) * 1e-5]
  lambda2 <- eigen(X2tX2)$values
  lambda2<-lambda2[lambda2 > mean(lambda2) * 1e-5]
  pvalue1 <- CompQuadForm::davies(Q1, lambda = lambda1)$Qq
  pvalue2 <- CompQuadForm::davies(Q2, lambda = lambda2)$Qq

  null_est<-admixLogica_PXEM_null(X1ty, X2ty, X1tX1, X2tX2, X1tX2, N,
                         sigma_g1_sq = sigma_g1_initial,
                         sigma_g2_sq = sigma_g2_initial,
                         sigma_e_sq = 1 - (sigma_g1_initial + sigma_g2_initial) * N_SNP,
                         max_iter = 200)
  
  res <- score_test_full(EUR_c, AFR_c, y,
                         sigma_g1_sq = null_est$sigma_g1_sq,
                         sigma_g2_sq = null_est$sigma_g2_sq,
                         sigma_e_sq  = 1-(null_est$sigma_g1_sq+null_est$sigma_g2_sq)*N_SNP)
  pvalue12<-res$p_two
  cat("Rho P value is ", pvalue12, "\n")
  res_summary<-score_test_summary(X1ty, X2ty, X1tX1, X2tX2, X1tX2,
                                  sigma_g1_sq = null_est$sigma_g1_sq,
                                  sigma_g2_sq = null_est$sigma_g2_sq,
                                  sigma_e_sq  = 1-(null_est$sigma_g1_sq+null_est$sigma_g2_sq)*N_SNP)
  pvalue12_summary<-res_summary$p_summary
  cat("Rho P value Summary is ",  pvalue12_summary, "\n")

  
  result_vector <- c(
    h1_est = admixres$sigma_g1_sq * N_SNP,
    h1_pval = pvalue1,
    h2_est = admixres$sigma_g2_sq * N_SNP,
    h2_pval = pvalue2,
    rho_est = admixres$rho_g,
    rho_pval = pvalue12,
    rho_pval_summary = pvalue12_summary
  #  rho_pval_llk = p_llk
  )
  return(result_vector)
}









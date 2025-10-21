// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List admixLogica_h2_initial(const arma::vec& Xty,
                 const arma::mat& XtX,
                 int N, 
                 double sigma_g_sq, 
                 double sigma_e_sq,
                 int max_iter = 100,
                 double tol = 1e-6){
  int num_SNP = XtX.n_rows;
  arma::vec eigval;
  arma::mat eigvec;
  
  eig_sym(eigval,eigvec,XtX);
  arma::vec UtXty = eigvec.t()*Xty;
  arma::vec llk(max_iter+1);
  llk(0) = -datum::inf;
  
  
  arma::vec lambda(num_SNP), UtXty_sigma(num_SNP);
  double sigma_g_sq_update = sigma_g_sq, sigma_e_sq_update = sigma_e_sq;
  for(int i=0;i<max_iter;i++){
    lambda = 1.0/(eigval/sigma_e_sq_update + 1.0/sigma_g_sq_update);
    UtXty_sigma = UtXty/sigma_e_sq_update;
    llk(i+1) = -N/2*log(sigma_e_sq_update)-num_SNP/2*log(sigma_g_sq_update)+accu(log(lambda))/2-N/2/sigma_e_sq_update+accu(pow(UtXty_sigma,2)%lambda)/2;
    if(abs(llk(i+1)-llk(i))<pow(10,-5)){
      llk.resize(i+2);
      break;
    }
    sigma_e_sq_update = 1 + ((accu(pow(UtXty_sigma%lambda,2)%eigval)+accu(lambda%eigval))-2.0*sigma_e_sq_update*accu(pow(UtXty_sigma,2)%lambda))/N;
    sigma_g_sq_update = (accu(pow(UtXty_sigma%lambda,2)) + accu(lambda))/num_SNP;
  }
  return List::create(Named("sigma_g_sq")=sigma_g_sq_update,
                      Named("sigma_e_sq")=sigma_e_sq_update,
                      Named("llk")=llk);
  
}
// [[Rcpp::export]]
List admixLogica_h2_joint(const arma::vec& X1ty,
                          const arma::vec& X2ty,
                          const arma::mat& X1tX1,
                          const arma::mat& X2tX2,
                          const arma::mat& X1tX2,
                          int N,
                          double sigma_g1_sq,
                          double sigma_g2_sq,
                          double sigma_e_sq,
                          int max_iter = 100,
                          double tol = 1e-6) {
  
  int p = X1ty.n_elem;
  
  // Concatenate Xty
  arma::vec Xty(2*p);
  Xty.subvec(0, p - 1) = X1ty;
  Xty.subvec(p, 2*p - 1) = X2ty;
  
  // Build full XtX matrix
  arma::mat XtX(2*p, 2*p);
  XtX.submat(0, 0, p-1, p-1) = X1tX1;
  XtX.submat(p, p, 2*p-1, 2*p-1) = X2tX2;
  XtX.submat(0, p, p-1, 2*p-1) = X1tX2;
  XtX.submat(p, 0, 2*p-1, p-1) = X1tX2.t();
  
  // Eigen-decomposition
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, XtX);
  arma::vec UtXty = eigvec.t() * Xty;
  
  arma::vec llk(max_iter + 1);
  llk(0) = -datum::inf;
  
  arma::vec lambda(2*p), UtXty_sigma(2*p);
  double sigma_g1_sq_update = sigma_g1_sq;
  double sigma_g2_sq_update = sigma_g2_sq;
  double sigma_e_sq_update = sigma_e_sq;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    // lambda = 1 / (eigenvalue / sigma_e_sq + 1/sigma_g1 + 1/sigma_g2)
    for (int j = 0; j < 2*p; ++j) {
      if (j < p) {
        lambda(j) = 1.0 / (eigval(j) / sigma_e_sq_update + 1.0 / sigma_g1_sq_update);
      } else {
        lambda(j) = 1.0 / (eigval(j) / sigma_e_sq_update + 1.0 / sigma_g2_sq_update);
      }
    }
    
    UtXty_sigma = UtXty / sigma_e_sq_update;
    llk(iter+1) = -N/2.0 * log(sigma_e_sq_update)
      - p/2.0 * log(sigma_g1_sq_update)
      - p/2.0 * log(sigma_g2_sq_update)
      + 0.5 * accu(log(lambda))
      - N / (2.0 * sigma_e_sq_update)
      + 0.5 * accu(pow(UtXty_sigma, 2) % lambda);
      
      if (std::abs(llk(iter+1) - llk(iter)) < tol) {
        llk.resize(iter+2);
        break;
      }
      
      double tr_g1 = 0.0, tr_g2 = 0.0;
      double num_g1 = 0.0, num_g2 = 0.0;
      for (int j = 0; j < 2*p; ++j) {
        double eig = eigval(j);
        double lmb = lambda(j);
        double us = UtXty_sigma(j);
        double tmp = pow(us * lmb, 2);
        if (j < p) {
          tr_g1 += lmb;
          num_g1 += tmp;
        } else {
          tr_g2 += lmb;
          num_g2 += tmp;
        }
      }
      
      sigma_g1_sq_update = (num_g1 + tr_g1) / p;
      sigma_g2_sq_update = (num_g2 + tr_g2) / p;
      
      double num_e = accu(pow(UtXty_sigma % lambda, 2) % eigval);
      double tr_e = accu(lambda % eigval);
      sigma_e_sq_update = 1.0 + (num_e + tr_e - 2.0 * sigma_e_sq_update * accu(pow(UtXty_sigma, 2) % lambda)) / N;
  }
  
  return List::create(Named("sigma_g1_sq") = sigma_g1_sq_update,
                      Named("sigma_g2_sq") = sigma_g2_sq_update,
                      Named("sigma_e_sq") = sigma_e_sq_update,
                      Named("llk") = llk);
}
// [[Rcpp::export]]
List admixLogica_h2_joint_PXEM(const arma::vec& X1ty,
                               const arma::vec& X2ty,
                               const arma::mat& X1tX1,
                               const arma::mat& X2tX2,
                               const arma::mat& X1tX2,
                               int N,
                               double sigma_g1_sq,
                               double sigma_g2_sq,
                               double sigma_e_sq,
                               int max_iter = 100,
                               double tol = 1e-6) {
  int p = X1ty.n_elem;
  
  // 1) Build joint XtX and Xty
  arma::mat XtX(2*p,2*p);
  XtX.submat(0,0,p-1,p-1) = X1tX1;
  XtX.submat(p,p,2*p-1,2*p-1) = X2tX2;
  XtX.submat(0,p,p-1,2*p-1) = X1tX2;
  XtX.submat(p,0,2*p-1,p-1) = X1tX2.t();
  arma::vec Xty(2*p);
  Xty.subvec(0,p-1) = X1ty;
  Xty.subvec(p,2*p-1) = X2ty;
  
  // 2) Eigen-decompose XtX once
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, XtX);
  arma::vec Utx = eigvec.t() * Xty;
  
  // 3) PX parameters
  double gamma1 = 1.0, gamma2 = 1.0;
  arma::vec lambda(2*p), UtXty_sigma(2*p);
  arma::vec llk(max_iter+1);
  llk(0) = -datum::inf;
  
  // 4) PX-EM loop in eigen-space
  for(int iter=0; iter<max_iter; ++iter){
    // E-step: compute diagonal of posterior covariance
    for(int j=0; j<2*p; ++j){
      bool is1 = (j < p);
      double prior_prec = is1 
      ? 1.0/(gamma1 * sigma_g1_sq) 
        : 1.0/(gamma2 * sigma_g2_sq);
      lambda(j) = 1.0 / (eigval(j)/sigma_e_sq + prior_prec);
    }
    UtXty_sigma = Utx / sigma_e_sq;
    
    // llk (optional monitor)
    double logdet = accu(log(lambda));
    double quad   = accu(square(UtXty_sigma) % lambda);
    llk(iter+1) = -N/2.0*log(sigma_e_sq)
      -p/2.0*log(sigma_g1_sq)
      -p/2.0*log(sigma_g2_sq)
      +0.5*logdet
      -N/(2.0*sigma_e_sq)
      +0.5*quad;
      if(std::abs(llk(iter+1)-llk(iter)) < tol){
        llk.resize(iter+2);
        break;
      }
      
      // M-step: accumulate moments
      double sum1=0, tr1=0, sum2=0, tr2=0, se_num=0, se_tr=0;
      for(int j=0; j<2*p; ++j){
        double d = eigval(j), l = lambda(j), u = UtXty_sigma(j);
        double u2l2 = std::pow(u*l,2);
        if(j < p){ sum1 += u2l2; tr1 += l; }
        else      { sum2 += u2l2; tr2 += l; }
        se_num += u2l2 * d;
        se_tr  += l * d;
      }
      
      // PX expansions
      double g1_exp = (sum1 + tr1)/p;
      double g2_exp = (sum2 + tr2)/p;
      gamma1 = g1_exp / sigma_g1_sq;
      gamma2 = g2_exp / sigma_g2_sq;
      sigma_g1_sq = g1_exp / gamma1;
      sigma_g2_sq = g2_exp / gamma2;
      
      // residual update
      sigma_e_sq = 1.0 + (se_num + se_tr - 2.0*sigma_e_sq*accu(square(UtXty_sigma)%lambda))/N;
  }
  
  return List::create(
    _["sigma_g1_sq"] = sigma_g1_sq,
    _["sigma_g2_sq"] = sigma_g2_sq,
    _["sigma_e_sq"]  = sigma_e_sq,
    _["llk"]         = llk
  );
}

// [[Rcpp::export]]
double llk_rho(const arma::vec& X1ty,
               const arma::vec& X2ty,
               const arma::mat& X1tX1,
               const arma::mat& X2tX2,
               const arma::mat& X1tX2,
               int N,
               double sigma_g1_sq,
               double sigma_g2_sq,
               double sigma_e_sq,
               double rho_g) {
  int p = X1ty.n_elem;
  // 构造 XtX 和 Xty
  arma::mat XtX(2*p, 2*p);
  XtX.submat(0,0,    p-1,p-1)     = X1tX1;
  XtX.submat(p,  p,  2*p-1,2*p-1) = X2tX2;
  XtX.submat(0,  p,  p-1,2*p-1)   = X1tX2;
  XtX.submat(p,  0,  2*p-1,p-1)   = X1tX2.t();
  arma::vec Xty(2*p);
  Xty.subvec(0,   p-1)  = X1ty;
  Xty.subvec(p, 2*p-1)  = X2ty;
  
  // 构造 Σ_g^{-1}
  double denom = sigma_g1_sq * sigma_g2_sq * (1 - rho_g*rho_g);
  double s11   = sigma_g2_sq / denom;
  double s22   = sigma_g1_sq / denom;
  double s12   = -rho_g * std::sqrt(sigma_g1_sq * sigma_g2_sq) / denom;
  arma::mat I_p = arma::eye(p,p);
  arma::mat Sg_inv(2*p,2*p, fill::zeros);
  Sg_inv.submat(0,0,    p-1,p-1)     = s11 * I_p;
  Sg_inv.submat(p,  p,  2*p-1,2*p-1) = s22 * I_p;
  Sg_inv.submat(0,  p,  p-1,2*p-1)   = s12 * I_p;
  Sg_inv.submat(p,  0,  2*p-1,p-1)   = s12 * I_p;
  
  arma::mat V = inv_sympd(Sg_inv + XtX / sigma_e_sq);
  arma::vec m = V * (Xty / sigma_e_sq);
  
  // 计算对数似然： 
  //  L = -N/2*log(σ_e^2) - p/2*log(σ_g1^2) - p/2*log(σ_g2^2)
  //      +1/2 log det(V) -N/(2σ_e^2) +1/2 mᵀ Sg_inv m
  double ll  = - (double)N/2 * std::log(sigma_e_sq)
    - p/2.0 * std::log(sigma_g1_sq)
    - p/2.0 * std::log(sigma_g2_sq)
    + 0.5 * std::log(det(V))
    - (double)N/(2.0*sigma_e_sq)
    + 0.5 * as_scalar(m.t() * Sg_inv * m);
    return ll;
}

// [[Rcpp::export]]
List admixLogica(const arma::vec& X1ty,
              const arma::vec& X2ty,
              const arma::mat& X1tX1,
              const arma::mat& X2tX2,
              const arma::mat& X1tX2,
              int N, int max_iter,
              double sigma_g1_sq, double sigma_g2_sq,
              double rho_g, double sigma_e_sq,
              double tol = 1e-6){
  
  int p = X1ty.n_elem;
  
  arma::mat I_p = arma::eye(p, p);
  
  double diff = 1.0;
  int iter = 0;
  
  arma::mat Sigma_g(2*p, 2*p);
  arma::mat Sigma_g_inv(2*p, 2*p);
  
  arma::mat XtX(2*p, 2*p);
  XtX.submat(0,0,p-1,p-1) = X1tX1;
  XtX.submat(p,p,2*p-1,2*p-1) = X2tX2;
  XtX.submat(0,p,p-1,2*p-1) = X1tX2;
  XtX.submat(p,0,2*p-1,p-1) = X1tX2.t();
  
  arma::vec Xty(2*p);
  Xty.subvec(0,p-1) = X1ty;
  Xty.subvec(p,2*p-1) = X2ty;
  
  arma::mat Sigma_beta_y(2*p,2*p);
  arma::vec mu_beta_y(2*p);
  
  while(diff > tol && iter < max_iter){
    
    // E-step (optimized inverse):
    double denom = sigma_g1_sq * sigma_g2_sq * (1 - rho_g * rho_g);
    double s11 = sigma_g2_sq / denom;
    double s22 = sigma_g1_sq / denom;
    double s12 = -rho_g * sqrt(sigma_g1_sq * sigma_g2_sq) / denom;
    
    Sigma_g_inv.zeros();
    Sigma_g_inv.submat(0,0,p-1,p-1) = s11 * I_p;
    Sigma_g_inv.submat(p,p,2*p-1,2*p-1) = s22 * I_p;
    Sigma_g_inv.submat(0,p,p-1,2*p-1) = s12 * I_p;
    Sigma_g_inv.submat(p,0,2*p-1,p-1) = s12 * I_p;
    
    Sigma_beta_y = inv_sympd(Sigma_g_inv + XtX/sigma_e_sq);
    mu_beta_y = Sigma_beta_y * (Xty/sigma_e_sq);
    
    arma::vec mu_b1 = mu_beta_y.subvec(0,p-1);
    arma::vec mu_b2 = mu_beta_y.subvec(p,2*p-1);
    arma::mat Sigma_b1 = Sigma_beta_y.submat(0,0,p-1,p-1);
    arma::mat Sigma_b2 = Sigma_beta_y.submat(p,p,2*p-1,2*p-1);
    arma::mat Sigma_b12 = Sigma_beta_y.submat(0,p,p-1,2*p-1);
    
    // M-step updates:
    double sigma_g1_sq_new = trace(Sigma_b1 + mu_b1*mu_b1.t())/p;
    double sigma_g2_sq_new = trace(Sigma_b2 + mu_b2*mu_b2.t())/p;
    double rho_g_new = trace(Sigma_b12 + mu_b1*mu_b2.t()) / 
      sqrt(trace(Sigma_b1 + mu_b1*mu_b1.t()) * trace(Sigma_b2 + mu_b2*mu_b2.t()));
    
    double sigma_e_sq_new = (1.0 -
                             2.0 * dot(mu_beta_y, Xty)/N +
                             trace(XtX * (Sigma_beta_y + mu_beta_y * mu_beta_y.t()))/N);
    
    diff = std::abs(sigma_g1_sq_new-sigma_g1_sq) +
      std::abs(sigma_g2_sq_new-sigma_g2_sq) +
      std::abs(rho_g_new-rho_g) +
      std::abs(sigma_e_sq_new-sigma_e_sq);
    
    sigma_g1_sq = sigma_g1_sq_new;
    sigma_g2_sq = sigma_g2_sq_new;
    rho_g = rho_g_new;
    sigma_e_sq = sigma_e_sq_new;
    
    iter++;
  }
  
  return List::create(Named("sigma_g1_sq")=sigma_g1_sq,
                      Named("sigma_g2_sq")=sigma_g2_sq,
                      Named("rho_g")=rho_g,
                      Named("sigma_e_sq")=sigma_e_sq,
                      Named("iter")=iter);
}

// [[Rcpp::export]]
List admixLogica_PXEM(const arma::vec& X1ty,
                      const arma::vec& X2ty,
                      const arma::mat& X1tX1,
                      const arma::mat& X2tX2,
                      const arma::mat& X1tX2,
                      int N, int max_iter,
                      double sigma_g1_sq, 
                      double sigma_g2_sq,
                      double rho_g, 
                      double sigma_e_sq,
                      double tol = 1e-6) {
  
  int p = X1ty.n_elem;
  arma::mat I_p = arma::eye(p, p);
  
  // Build full XtX and Xty
  arma::mat XtX(2*p,2*p);
  XtX.submat(0,0,    p-1,p-1)    = X1tX1;
  XtX.submat(p,  p,  2*p-1,2*p-1)= X2tX2;
  XtX.submat(0,  p,  p-1,2*p-1)  = X1tX2;
  XtX.submat(p,  0,  2*p-1,p-1)  = X1tX2.t();
  
  arma::vec Xty(2*p);
  Xty.subvec(0,   p-1)   = X1ty;
  Xty.subvec(p, 2*p-1)   = X2ty;
  
  // PX‐EM expansions
  double gamma1 = 1.0, gamma2 = 1.0;
  
  arma::mat Sigma_g_inv(2*p,2*p);
  arma::mat Sigma_beta_y(2*p,2*p);
  arma::vec mu_beta_y(2*p);
  
  double llk_old = -datum::inf, llk_new;

  for( int iter = 0; iter < max_iter; ++iter){
    //
    // E-step with expanded prior variances
    //
    double vg1 = gamma1 * sigma_g1_sq;
    double vg2 = gamma2 * sigma_g2_sq;
    double denom = vg1 * vg2 * (1 - rho_g*rho_g);
    double s11 = vg2 / denom;
    double s22 = vg1 / denom;
    double s12 = -rho_g * std::sqrt(vg1*vg2) / denom;
    
    Sigma_g_inv.zeros();
    Sigma_g_inv.submat(0,0,    p-1,p-1)    = s11 * I_p;
    Sigma_g_inv.submat(p,  p,  2*p-1,2*p-1)= s22 * I_p;
    Sigma_g_inv.submat(0,  p,  p-1,2*p-1)  = s12 * I_p;
    Sigma_g_inv.submat(p,  0,  2*p-1,p-1)  = s12 * I_p;
    
    Sigma_beta_y = inv_sympd(Sigma_g_inv + XtX / sigma_e_sq);
    mu_beta_y    = Sigma_beta_y * (Xty / sigma_e_sq);
    
    //
    // M-step: compute posterior 2nd moments
    //
    arma::vec mu_b1 = mu_beta_y.subvec(0, p-1);
    arma::vec mu_b2 = mu_beta_y.subvec(p, 2*p-1);
    arma::mat S_b1  = Sigma_beta_y.submat(0,0,    p-1,p-1);
    arma::mat S_b2  = Sigma_beta_y.submat(p,p,    2*p-1,2*p-1);
    arma::mat S_b12 = Sigma_beta_y.submat(0,p,    p-1,2*p-1);
    
    double exp_vg1 = (trace(S_b1  + mu_b1*mu_b1.t())) / p;
    double exp_vg2 = (trace(S_b2  + mu_b2*mu_b2.t())) / p;
    double exp_cov=  trace(S_b12 + mu_b1*mu_b2.t()) / p;
    
    // Project back to original scale
    double gamma1_new = exp_vg1 / sigma_g1_sq;
    double gamma2_new = exp_vg2 / sigma_g2_sq;
    double rho_g_new  = exp_cov / std::sqrt(exp_vg1*exp_vg2);
    
    sigma_g1_sq = exp_vg1 / gamma1_new;
    sigma_g2_sq = exp_vg2 / gamma2_new;
    rho_g       = rho_g_new;
    gamma1      = gamma1_new;
    gamma2      = gamma2_new;
    
    // Residual variance update
    double term = trace(XtX * (Sigma_beta_y + mu_beta_y*mu_beta_y.t())) / N;
    sigma_e_sq = 1.0
    - 2.0 * dot(mu_beta_y, Xty) / N
    + term;
    
    //
    // compute current llk
    //
    double logdetV, signV;
    log_det(logdetV, signV, Sigma_beta_y);
    double quad = as_scalar(mu_beta_y.t() * Sigma_g_inv * mu_beta_y);
    llk_new = 
      - double(N)/2.0 * std::log(sigma_e_sq)
      - double(p)/2.0 * std::log(sigma_g1_sq)
      - double(p)/2.0 * std::log(sigma_g2_sq)
      + 0.5 * logdetV
      - double(N)/(2.0 * sigma_e_sq)
      + 0.5 * quad;
      
      double diff_llk = std::abs(llk_new - llk_old);
      if(diff_llk < tol) break;
      llk_old = llk_new;
  }
  
  return List::create(
    Named("sigma_g1_sq") = sigma_g1_sq,
    Named("sigma_g2_sq") = sigma_g2_sq,
    Named("rho_g")       = rho_g,
    Named("sigma_e_sq")  = sigma_e_sq,
    Named("gamma1")      = gamma1,
    Named("gamma2")      = gamma2,
    Named("llk")         = llk_new
  );
}

// [[Rcpp::export]]
double llk_uni_h2(double h2,
                  const arma::vec& eigval,
                  const arma::vec& utXty,
                  int N) {
  int p = eigval.n_elem;
  // parameterize σ_g² = h2, σ_e² = 1−h2
  double sg2 = h2;
  double se2 = 1.0 - h2;
  
  // constant terms
  double ll = 
    - double(N)/2.0 * std::log(se2)
    - double(p)/2.0 * std::log(sg2)
    - double(N)/(2.0 * se2);
    
    // diagonal posterior variances + quadratic
    double sum_log = 0.0, sum_quad = 0.0;
    for(int j=0; j<p; ++j){
      // prior precision = 1/σ_g²
      double prior_prec = 1.0/sg2;
      // posterior variance in eigen‐space
      double lam = 1.0 / (eigval[j]/se2 + prior_prec);
      sum_log  += std::log(lam);
      double u = utXty[j] / se2;
      sum_quad += u*u * lam;
    }
    return ll + 0.5*(sum_log + sum_quad);
}




// [[Rcpp::export]]
List score_test_full(const arma::mat& X1,
                     const arma::mat& X2,
                     const arma::vec&  y,
                     double sigma_g1_sq,
                     double sigma_g2_sq,
                     double sigma_e_sq) {
  int n = y.n_elem;
  // 1) Null covariance V0 (assume centered)
  arma::mat V0 = sigma_g1_sq*(X1*X1.t())
    + sigma_g2_sq*(X2*X2.t())
    + sigma_e_sq*arma::eye(n,n);
    
    // 2) Build C = ∂V/∂ρ at ρ=0
    arma::mat C = std::sqrt(sigma_g1_sq*sigma_g2_sq) 
      * (X1*X2.t() + X2*X1.t());
    
    // 3) Invert V0
    arma::mat iV0;
    if(!inv_sympd(iV0, V0))
      stop("V0 inversion failed");
    
    // 4) Build score‐matrix A = 0.5 * iV0 * C * iV0
    arma::mat A = 0.5 * (iV0 * C * iV0);
    
    // 5) Compute T = yᵀ A y
    double T = as_scalar(y.t() * (A * y));
    
    // 6) Compute V0^{-1/2} via eigen‐decomposition
    arma::vec d; arma::mat U;
    eig_sym(d, U, V0);    // V0 = U diag(d) Uᵀ
    arma::vec inv_sqrt_d = 1.0/sqrt(d);
    arma::mat iVsq = U * diagmat(inv_sqrt_d) * U.t();  // V0^{-1/2}
    
    // 7) Build M = 0.5 * iVsq * C * iVsq
    arma::mat M = 0.5 * (iVsq * C * iVsq);
    
    // 8) Get real eigenvalues of M
    arma::vec lam = eig_sym(M);
    // drop exact zeros if any
    lam = lam.elem(find(lam != 0.0));
    
    // 9) Call Davies
    Function davies = Environment::namespace_env("CompQuadForm")["davies"];
    List out = davies(T, wrap(lam), Named("acc")=1e-6);
    double p_up  = as<double>(out["Qq"]);
    double p_two = 2.0 * std::min(p_up, 1.0 - p_up);
    
    return List::create(
      _["T"]      = T,
      _["lambda"] = lam,
      _["p_one"]  = p_up,
      _["p_two"]  = p_two
    );
}


// [[Rcpp::export]]
SEXP score_test_summary(arma::vec& Xty_1,
                        arma::vec& Xty_2,
                        const arma::mat& XtX_1,
                        const arma::mat& XtX_2,
                        const arma::mat& X1tX2,
                        double sigma_g1_sq,
                        double sigma_g2_sq,
                        double sigma_e_sq){
  
  arma::vec eigval_1,eigval_2;
  arma::mat eigvec_1,eigvec_2;
  
  eig_sym(eigval_1,eigvec_1,XtX_1);
  eig_sym(eigval_2,eigvec_2,XtX_2);
  
  arma::mat eigvec_1_trans = trans(eigvec_1);
  arma::mat eigvec_2_trans = trans(eigvec_2);
  // Define Xty, UtXty, and mu_beta, lambda
  
  arma::vec UtXty_1 = eigvec_1.t()*Xty_1;
  arma::vec UtXty_2 = eigvec_2.t()*Xty_2;
  
  Xty_1 = Xty_1/sigma_e_sq;
  Xty_2 = Xty_2/sigma_e_sq;
  
  int num_SNP = XtX_1.n_rows;
  arma::mat XtX_1_quad_1_XtX_1(num_SNP,num_SNP),XtX_2_quad_2_XtX_2(num_SNP,num_SNP),B_mat(num_SNP,num_SNP);
  
  
  XtX_1_quad_1_XtX_1 = eigvec_1*diagmat(eigval_1/sigma_e_sq - pow(eigval_1/sigma_e_sq,2)/(eigval_1/sigma_e_sq + 1.0/sigma_g1_sq))*eigvec_1.t();
  XtX_2_quad_2_XtX_2 = eigvec_2*diagmat(eigval_2/sigma_e_sq - pow(eigval_2/sigma_e_sq,2)/(eigval_2/sigma_e_sq + 1.0/sigma_g2_sq))*eigvec_2.t();
  double rho_score = dot(Xty_1.t()-UtXty_1.t()/sigma_e_sq*(eigvec_1_trans.each_col()%(eigval_1/sigma_e_sq/(eigval_1/sigma_e_sq + 1.0/sigma_g1_sq))),Xty_2.t()-UtXty_2.t()/sigma_e_sq*(eigvec_2_trans.each_col()%(eigval_2/sigma_e_sq/(eigval_2/sigma_e_sq + 1.0/sigma_g2_sq))));
  
  B_mat = XtX_1_quad_1_XtX_1 * XtX_2_quad_2_XtX_2;
  
  cx_vec eigval =eig_gen(B_mat); 
  vec evals = real(eigval);
  vec pos_evals = sqrt(evals(find(evals>0)))/2.0;
  pos_evals = join_cols( pos_evals,-pos_evals);

  // 8) Davies for two‐sided p‐value
  double T = std::abs(rho_score);
  Function davies = Environment::namespace_env("CompQuadForm")["davies"];
  List dres = davies(T, wrap(pos_evals), Named("acc") = 1e-6);
  double p_one = as<double>(dres["Qq"]);
  double p_two = 2.0 * p_one;
  
  
  List output = List::create(
    _["rho_score"] = rho_score,
    _["rho_eigvals"] = pos_evals,
    _["p_summary"]  = p_two
  );
  return(output);
}


// [[Rcpp::export]]
List admixLogica_PXEM_null(const arma::vec& X1ty,
                           const arma::vec& X2ty,
                           const arma::mat& X1tX1,
                           const arma::mat& X2tX2,
                           const arma::mat& X1tX2,
                           int N,
                           double sigma_g1_sq,
                           double sigma_g2_sq,
                           double sigma_e_sq,
                           int max_iter=100,
                           double tol = 1e-6) {
  int p = X1ty.n_elem;
  arma::mat I_p = arma::eye(p,p);
  
  // Build big XtX and Xty
  arma::mat XtX(2*p,2*p,fill::zeros);
  XtX.submat(0,0,    p-1,p-1)     = X1tX1;
  XtX.submat(p,p,  2*p-1,2*p-1)   = X2tX2;
  XtX.submat(0,p,  p-1,2*p-1)     = X1tX2;
  XtX.submat(p,0,  2*p-1,p-1)     = X1tX2.t();
  arma::vec Xty(2*p);
  Xty.subvec(0,   p-1) = X1ty;
  Xty.subvec(p, 2*p-1) = X2ty;
  
  // PX‐EM “expansion” parameters
  double gamma1 = 1.0, gamma2 = 1.0;
  
  arma::mat Sigma_g_inv(2*p,2*p);
  arma::mat Sigma_beta(2*p,2*p);
  arma::vec mu_beta(2*p);
  
  double llk_old = -datum::inf, llk_new = 0.0;
  
  int iter=0;
  for(; iter<max_iter; ++iter) {
    // E‐step: rho=0 => Sigma_g is block‐diag
    double vg1 = gamma1 * sigma_g1_sq;
    double vg2 = gamma2 * sigma_g2_sq;
    // inverse of block‐diag [vg1 I, vg2 I]
    Sigma_g_inv.zeros();
    Sigma_g_inv.submat(0,0,    p-1,p-1)     = (1.0/vg1) * I_p;
    Sigma_g_inv.submat(p,p,    2*p-1,2*p-1) = (1.0/vg2) * I_p;
    
    // posterior cov and mean
    Sigma_beta = inv_sympd(Sigma_g_inv + XtX/sigma_e_sq);
    mu_beta    = Sigma_beta * (Xty/sigma_e_sq);
    
    // M‐step: second moments
    arma::vec mu1 = mu_beta.subvec(0,   p-1);
    arma::vec mu2 = mu_beta.subvec(p, 2*p-1);
    arma::mat S1  = Sigma_beta.submat(0,0,    p-1,p-1);
    arma::mat S2  = Sigma_beta.submat(p,p,  2*p-1,2*p-1);
    
    double exp_vg1 = trace(S1 + mu1*mu1.t())/p;
    double exp_vg2 = trace(S2 + mu2*mu2.t())/p;
    
    // update gammas (PX step)
    double gamma1_new = exp_vg1 / sigma_g1_sq;
    double gamma2_new = exp_vg2 / sigma_g2_sq;
    
    // project back onto original scale
    sigma_g1_sq = exp_vg1 / gamma1_new;
    sigma_g2_sq = exp_vg2 / gamma2_new;
    gamma1      = gamma1_new;
    gamma2      = gamma2_new;
    
    // residual variance
    double term = trace(XtX * (Sigma_beta + mu_beta*mu_beta.t()))/N;
    sigma_e_sq = 1.0 
    - 2.0 * dot(mu_beta,Xty)/N
    + term;
    
    // log‐likelihood (optional, for convergence)
    double logdet, sign;
    log_det(logdet,sign, Sigma_beta);
    double quad = as_scalar(mu_beta.t() * Sigma_g_inv * mu_beta);
    llk_new = 
      - double(N)/2.0 * log(sigma_e_sq)
      - double(p)/2.0 * log(sigma_g1_sq)
      - double(p)/2.0 * log(sigma_g2_sq)
      + 0.5*logdet
      - double(N)/(2.0*sigma_e_sq)
      + 0.5*quad;
      
      if(std::abs(llk_new - llk_old) < tol) break;
      llk_old = llk_new;
  }
  
  return List::create(
    _["sigma_g1_sq"] = sigma_g1_sq,
    _["sigma_g2_sq"] = sigma_g2_sq,
    _["sigma_e_sq"]  = sigma_e_sq,
    _["gamma1"]      = gamma1,
    _["gamma2"]      = gamma2,
    _["llk"]         = llk_new,
    _["iter"]        = iter
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Matrix Inverse
//' 
//' Calcualtes \eqn{A^{-1}}.
//'
//' @param A Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP matInv(const arma::mat A){
  const arma::mat Ai = arma::pinv(A);
  return Rcpp::wrap(Ai);
}

//' Matrix Inner Product
//'
//' Calculates the product \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP matIP(const arma::mat A, const arma::mat B){
  const arma::mat AtB = A.t()*B;
  return Rcpp::wrap(AtB);
}

//' Matrix Matrix Product
//'
//' Calculates the product \eqn{AB}. 
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP MMP(const arma::mat A, const arma::mat B){
  const arma::mat C = A*B;
  return Rcpp::wrap(C);
}

//' Matrix Outer Product
//' 
//' Calculates the outer product \eqn{AB'}.
//' 
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP matOP(const arma::mat A, const arma::mat B){
  const arma::mat ABt = A*B.t();
  return Rcpp::wrap(ABt);
}

//' Quadratic Form
//' 
//' Calculates the quadratic form \eqn{X'AX}.
//' 
//' @param X Numeric matrix.
//' @param A Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP matQF(const arma::mat X, const arma::mat A){
  const arma::mat xAx = X.t()*A*X;
  return Rcpp::wrap(xAx);
}
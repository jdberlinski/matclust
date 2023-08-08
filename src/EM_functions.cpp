#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::uvec cl);
double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::mat z);
double log_f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p);
double f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p);
double f(arma::mat xi, arma::vec pr, arma::mat mu, arma::cube sig, int R, int p, int K);
arma::mat make_mask(arma::urowvec inds, int ncols);
List em_step(arma::cube x, arma::mat mu, arma::cube Sigma,  arma::mat z, arma::vec pr, arma::vec cl,
    arma::cube A, int n, int K, int R, int p, int iter);

// function for obtaining the log-liklihood for the matrix-variate normal
// distribution given a set of data, using a mixture model
//
// double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::uvec cl) {
// [[Rcpp::export]]
double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::mat z) {

  // TODO: this should probably include the class probabilities?
  arma::uword n = x.n_slices;
  double ll = 0.0;

  // for (arma::uword i = 0; i < n; i++)
  //   ll += log_f_k(x.slice(i), mu.row(cl(i)), sig.slice(cl(i)), R, p);
  //
  // this is a bandaid :)
  for (arma::uword i = 0; i < n; i++)
    ll += log(f(x.slice(i), arma::conv_to<arma::colvec>::from(z.row(i)), mu, sig, R, p, sig.n_slices));

  return ll;
}

double f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p) {
  // double result = exp(log_f_k(xi, mu, sig, R, p));
  double result, quad_trace;

  arma::mat M(R, p, arma::fill::zeros);
  M.each_row() += mu;

  // arma::mat sig_inv = arma::inv_sympd(arma::symmatu(sig));
  arma::mat sig_inv = arma::inv(sig);
  quad_trace = arma::trace((xi - M) * sig_inv * (xi - M).t());

  result = pow(2.0*M_PI, -2*R*p) * pow(exp(arma::log_det_sympd(arma::symmatu(sig))), -0.5*p) * exp(-0.5*quad_trace);

  return result;
}

double f(arma::mat xi, arma::vec pr, arma::mat mu, arma::cube sig, int R, int p, int K) {
  double acc = 0.0;

  for (arma::uword k = 0; k < K; k++)
    acc += pr(k) * f_k(xi, mu.row(k), sig.slice(k), R, p);

  return acc;
}

// [[Rcpp::export]]
double log_f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p) {
  // mean matrix
  // TODO: is there a better way to create a matrix that has identical rows?
  //       perhaps one way is to create wiht identical columns and transpose it.
  /* arma::mat M = arma::vec(R, arma::fill::ones) * mu.t(); */
  arma::mat M(R, p, arma::fill::zeros);
  M.each_row() += mu;

  arma::mat sig_inv = arma::inv_sympd(arma::symmatu(sig));

  double quad_trace = arma::trace((xi - M) * sig_inv * (xi - M).t());

  double result = (-R * p * 0.5) * std::log(2.0 * M_PI) +
    (-p * 0.5) * arma::log_det_sympd(arma::symmatu(sig)) +
    (-0.5 * quad_trace);

  return result;
}

arma::mat make_mask(arma::urowvec inds, int ncols) {
  arma::mat e(inds.n_elem, ncols, arma::fill::zeros);
  for (arma::uword i = 0; i < inds.n_elem; i++)
    e(i, inds(i)) = 1;

  return e;
}

// main function for EM
//
// [[Rcpp::export]]
List em_step(
  arma::cube x,      // data
  arma::mat mu,      // cluster means
  arma::cube Sigma,  // cluster covariances
  arma::mat z,       // cluster probs. (for each observation)
  arma::vec pr,      // cluster probs. (overall)
  arma::uvec cl,     // cluster membership
  arma::cube A,      // missingness indicator for x
  int n,             // number of observations
  int K,             // number of clusters
  int R,             // number of rows in data
  int p,             // number of columns in data
  int iter           // iteration number
) {

  // loop variables
  arma::uword i, j, k;

  // denominator and accumulators for updating mean and covariance
  double den;
  arma::rowvec acc(p, arma::fill::zeros);
  arma::mat sacc(p, p, arma::fill::zeros);

  // matrices to be used in EM updates
  arma::mat M(R, p);     // current mean matrix
  arma::mat m_k(R, p);   // work mean matrix, used for updating condtitional expectations
  arma::mat ec;          // indicator matrix to extract relevant values for E(X'X)
  arma::mat diff(R, p);
  arma::uvec miss;       // index of missing values (column-major)
  arma::uvec nmiss;      // index of nonmissing values (column-major)
  arma::umat miss_ind;   // matrix indices of missing values

  arma::mat Phi;         // conditional covariance

  arma::cube S(p*R, p*R, K); // covariance matrix of vec(X)


  // constant identity matrix for row covariance
  arma::mat I(R, R, arma::fill::eye);

  double ll, bic;

  // (E step)
  //   update class memebership probabilities,
  //   to be used for caculating conditional expectations
  for (i = 0; i < n; i++) {
    for (k = 0; k < K; k++) {
      z(i, k) = pr(k) * f_k(x.slice(i), mu.row(k), Sigma.slice(k), R, p);
      S.slice(k) = arma::kron(Sigma.slice(k), I);
    }
  }
  // NOTE: is there a precision problem here?

  z = arma::normalise(z, 1, 1);

  // M step
  for (k = 0; k < K; k++) {
    acc.zeros();
    sacc.zeros();
    M.zeros();
    M.each_row() += mu.row(k);

    den = arma::sum(z.col(k));

    for (i = 0; i < n; i++) {
      // get missing and nonmissing indices in various forms
      miss = arma::find(A.slice(i));
      nmiss = arma::find(A.slice(i) < 0.5);
      miss_ind = arma::ind2sub(size(A.slice(i)), miss);

      // impute the missing elements of x_i with their conditional expectations
      // given the observed elements
      x.slice(i).elem(miss).zeros();
      for (j = 0; j < K; j++) {
        m_k.zeros();
        m_k.each_row() += mu.row(j);
        x.slice(i).elem(miss) += z(i, j) * (m_k.elem(miss) + S.slice(j).submat(miss, nmiss) * arma::inv_sympd(arma::symmatu(S.slice(j).submat(nmiss, nmiss))) * (x.slice(i).elem(nmiss) - m_k.elem(nmiss)));
      }

      // conditional variance of the missing portion of x_i given the observed
      // portion. used to calculate E(X'X) in covariance estimation
      Phi = S.slice(k).submat(miss, miss) - S.slice(k).submat(miss, nmiss) * S.slice(k).submat(nmiss, nmiss).i() * S.slice(k).submat(nmiss, miss);

      acc += z(i, k) * arma::mean(x.slice(i), 0);
      diff = x.slice(i) - M;

      // covariance update requires the conditional variance to be added in
      // places corresponding to two missing values (if there are any missing
      // values)
      if (!miss.is_empty()) {
        ec = make_mask(miss_ind.row(1), p);
        sacc += z(i, k) * (diff.t() * diff + ec.t() * Phi * ec);
      } else {
        sacc += z(i, k) * diff.t() * diff;
      }
    }

    mu.row(k) = acc / den;
    Sigma.slice(k) = sacc / (R * den);
    pr(k) = den / n;
  }

  // assign "hard" clusters
  cl = arma::index_max(z, 1);

  ll = get_ll(x, mu, Sigma, R, p, z);
  bic = -2 * ll + log(n) * (K + K*p + K*p*(p + 1)/2);

  return List::create(
      Named("x") = x,
      Named("mu") = mu,
      Named("Sigma") = Sigma,
      Named("z") = z,
      Named("pr") = pr,
      Named("cl") = cl,
      Named("ll") = ll,
      Named("bic") = bic
  );
}

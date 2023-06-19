#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::uvec cl);
double log_f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p);
double f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p);
double f(arma::mat xi, arma::vec pr, arma::mat mu, arma::cube sig, int R, int p, int K);
arma::mat sweep(arma::mat A, arma::vec k);
arma::mat make_mask(arma::urowvec inds, int ncols);
List em_step(arma::cube x, arma::mat mu, arma::cube Sigma,  arma::mat z, arma::vec pr, arma::vec cl,
    arma::cube A, int n, int K, int R, int p, int iter);

// function for obtaining the log-liklihood for the matrix-variate normal
// distribution given a set of data, using a mixture model
//
// [[Rcpp::export]]
double get_ll(arma::cube x, arma::mat mu, arma::cube sig, int R, int p, arma::uvec cl) {

  // TODO: this should probably include the class probabilities?
  arma::uword n = x.n_slices;
  double ll = 0.0;

  for (arma::uword i = 0; i < n; i++)
    ll += log_f_k(x.slice(i), mu.row(cl(i)), sig.slice(cl(i)), R, p);

  return ll;
}

double f_k(arma::mat xi, arma::rowvec mu, arma::mat sig, int R, int p) {
  double result = exp(log_f_k(xi, mu, sig, R, p));
  return result;
}

double f(arma::mat xi, arma::vec pr, arma::mat mu, arma::cube sig, int R, int p, int K) {
  double acc = 0.0;

  for (arma::uword k = 0; k < K; k++)
    acc += pr(k) * f_k(xi, mu.row(k), sig.slice(k), R, p);

  return acc;
}

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

arma::mat sweep(arma::mat A, arma::uvec k) {

  double d, B;
  int i, j;

  for (int iter1 = 0; iter1 < k.n_elem; iter1++) {
    i = k(iter1);
    d = A(i, i);
    A.row(i) /= d;

    for (int iter2 = 0; iter2 < k.n_elem; iter2++) {
      if (iter1 == iter2) continue;
      j = k(iter2);
      B = A(j, i);

      A.row(j) -= B * A.row(i);
      A(j, i) = -B / d;
    }

    A(i, i) = 1.0 / d;
  }

  return A;
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
List em_step(arma::cube x,      // data
             arma::mat mu,      // cluster means
             arma::cube Sigma,  // cluster covariances
             arma::mat z,       // cluster probs. (for each observation)
             arma::vec pr,      // cluster probs. (overall)
             arma::uvec cl,      // cluster membership
             arma::cube A,      // missingness indicator for x
             int n,             // number of observations
             int K,             // number of clusters
             int R,             // number of rows in data
             int p,             // number of columns in data
             int iter
    ) {

  // loop variables
  arma::uword i, k;

  // accumulators and work variables
  double den;
  arma::rowvec acc(p, arma::fill::zeros);
  arma::mat sacc(p, p, arma::fill::zeros);

  // matrices to be used in EM updates
  arma::mat Sigma_i(p, p);
  arma::mat M(R, p);
  arma::mat S(p*R, p*R);
  arma::mat Rmat(p*R, p*R);
  arma::mat Rm;
  arma::mat Scm;
  arma::mat ec;
  arma::mat diff(R, p);
  arma::uvec miss;
  arma::uvec nmiss;
  arma::umat miss_ind;

  // constant matrix
  arma::mat I(R, R, arma::fill::eye);

  double ll, bic;

  // E step: update class membership probabilities
  for (k = 0; k < K; k++)
    for (i = 0; i < n; i++)
      z(i, k) = pr(k) * f_k(x.slice(i), mu.row(k), Sigma.slice(k), R, p) / f(x.slice(i), pr, mu, Sigma, R, p, K);

  if (iter == 1)
    cl = arma::index_max(z, 1);

  // M step
  for (k = 0; k < K; k++) {
    acc.zeros();
    sacc.zeros();
    M.zeros();
    M.each_row() += mu.row(k);

    Sigma_i = arma::inv_sympd(arma::symmatu(Sigma.slice(k)));

    S = arma::kron(Sigma_i, I);

    den = arma::sum(z.col(k));

    for (i = 0; i < n; i++) {
      // get missing and nonmissing indices in various forms
      miss = arma::find(A.slice(i));
      nmiss = arma::find(A.slice(i) < 0.5);
      miss_ind = arma::ind2sub(size(A.slice(i)), miss);

      // sweep the missing rows of S
      Rmat = sweep(S, miss);

      // update the missing values of x[,,i] to be the class mean
      if (cl(i) == k)
        x.slice(i).elem(miss) = M.elem(miss) + Rmat.submat(miss, nmiss) * (x.slice(i).elem(nmiss) - M.elem(nmiss));

      Rm = Rmat.submat(miss, miss);
      Scm = Sigma_i.submat(miss_ind.row(1), miss_ind.row(1));

      acc += z(i, k) * arma::mean(x.slice(i), 0);
      diff = x.slice(i) - M;

      if (!miss.is_empty()) {
        ec = make_mask(miss_ind.row(1), p);
        sacc += z(i, k) * (diff.t() * diff + ec.t() * (Rm * Scm) * ec);
      } else {
        sacc += z(i, k) * diff.t() * diff;
      }
    }

    mu.row(k) = acc / den;
    Sigma.slice(k) = sacc / (R * den);
    pr(k) = den / n;
  }

  cl = arma::index_max(z, 1);

  ll = get_ll(x, mu, Sigma, R, p, cl);
  bic = -2 * ll + log(n) * (K + K*p + K*p*(p + 1)/2);

  return List::create(Named("mu") = mu,
                      Named("Sigma") = Sigma,
                      Named("z") = z,
                      Named("pr") = pr,
                      Named("cl") = cl,
                      Named("ll") = ll,
                      Named("bic") = bic);
}

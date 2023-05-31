#include "arma.h"
using namespace Rcpp;

//' C++ function for computing AR coefficients from PACF.
//' See Section III in Barndorff-Nielsen and Schou (1973)
//' @references O. Barndorff-Nielsen and G. Schou (1973)
//' \emph{On the Parametrization of Autoregressive Models by Partial Autocorrelations}
//' Journal of Multivariate Analysis (3),408-419
//' <doi:10.1016/0047-259X(73)90030-4>
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix pacf2AR(NumericVector pacf) {
  //#include <cassert>
  unsigned p = pacf.size();
  NumericMatrix arCoef(p, p);
  if (p==0) return arCoef;
  arCoef(p-1,p-1) = pacf[p-1];
  if (p==1) return arCoef;
  NumericVector pacfTmp(p-1);
  for(unsigned j=0; j<p-1; ++j) pacfTmp[j] = pacf[j];
  NumericMatrix coefTmp = pacf2AR(pacfTmp);
  //assert(coefTmp.nrow() == p-1 && coefTmp.ncol == p-1);
  for(unsigned i=0; i < p-1; ++i) {
    for(unsigned j=0; j < p-1; ++j) {
      arCoef(i,j)=coefTmp(i,j);
    }
  }
  if (p==2) {
    arCoef(p-1,0) = pacf[0] * (1 - pacf[1]);
  } 
  if (p > 2) {
    for (int j=p-1; j >= 1; --j) {
      arCoef(p-1,j-1) = pacf[j-1];
      for (int r=1; r <= p-j; ++r) {
        arCoef(p-1,j-1) = arCoef(p-1,j-1) - pacf[j+r-1]*arCoef(j+r-2,r-1);
      }
    }
  }
  return arCoef;
}

//' Get epsilon process (i.e. model residuals) for ARMA(p,q)
//' @keywords internal
// [[Rcpp::export]]
NumericVector genEpsARMAC(NumericVector zt, NumericVector ar, NumericVector ma) {
  int m = zt.size(), p = ar.size(), q = ma.size();
  NumericVector epsilon_t(m + q - p), arsum(m - p), masum(m - q);
  for (int ii = 0; ii < q - p; ++ii) {
    epsilon_t[ii] = 0;
  }
  for (int tt = 0; tt < m - p; ++tt) {
    for (int jj = 0; jj < p; ++jj) {
      arsum[tt] += ar[jj] * zt[tt + p - jj - 1];
    }
    for (int kk = 0; kk < q; ++kk) {
      masum[tt] += ma[kk] * epsilon_t[tt + q - kk - 1];
    }
    epsilon_t[tt + q] = zt[tt + p] - arsum[tt] - masum[tt];
  }            
  NumericVector::const_iterator first = epsilon_t.begin() + q;
  NumericVector::const_iterator last = epsilon_t.begin() + m + q - p;
  NumericVector epsilon_s(first, last);         
  return epsilon_s;
}

//' ARMA(p,q) spectral density function
//' 
//' Evaluate the ARMA(p,q) spectral density at some frequencies freq in [0,pi),
//' Note that no test for model stationarity is performed.
//' @details See section 4.4 in the referenced book
//' @param freq numeric vector of frequencies to evaluate the psd, 0 <= freq < pi
//' @param ar autoregressive coefficients of ARMA model (use numeric(0) for empty AR part)
//' @param ma moving average coefficients of ARMA model (use numeric(0) for empty MA part)
//' @param sigma2 the model innovation variance
//' @return numeric vector of the (real-valued) spectral density values
//' @references P. J. Brockwell and R. Davis (1996)
//' \emph{Time Series: Theory and Methods (Second Edition)}
//' @export
// [[Rcpp::export]]
NumericVector psd_arma(NumericVector freq, NumericVector ar, NumericVector ma, double sigma2 = 1.0) {
  const unsigned n = freq.size();
  const unsigned p = ar.size();
  const unsigned q = ma.size();
  const double constant = sigma2 / (2*M_PI);
  NumericVector psd(n);
  for (unsigned j = 0; j < n; ++j) {
    const double lambda = freq[j];
    // compute numerator (MA part)
    std::complex<double> numerator_c(1.0, 0.0);
    for (unsigned i = 0; i < q; ++i) {
      numerator_c += ma[i]*std::polar<double>(1.0, -lambda*(double)(i+1));
    }
    // compute denominator (AR part)
    std::complex<double> denominator_c(1.0, 0.0);
    for (unsigned i = 0; i < p; ++i) {
      denominator_c -= ar[i]*std::polar<double>(1.0, -lambda*(double)(i+1));
    }
    psd[j] = constant * std::norm(numerator_c) / std::norm(denominator_c);
  }
  return(psd);
}

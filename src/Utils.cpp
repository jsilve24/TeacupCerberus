#include <CoDA.h>

using namespace Rcpp;

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// [[Rcpp::export]]
Eigen::MatrixXd tc2iqlrvar_internal(Eigen::MatrixXd Sigma,
                        int D1, int D2,
                        double qLow1=0.25, double qHigh1=0.75,
                        double qLow2=0.25, double qHigh2=0.75){
  int P = Sigma.rows();
  int N = Sigma.cols();
  if ( (N % P) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
  if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
  N = N/P; // safe after above validation

  MatrixXd O(P, P*N);
  MatrixXd V1, V2;
  VectorXd s;
  for (int n=0; n<N; n++){
    Map<MatrixXd> S(Sigma.middleCols(n*P, P).data(), P,P);
    s = S.diagonal().head(D1);
    V1 = coda::var2iqlrContrast(s, qLow1, qHigh1);
    s = S.diagonal().tail(D2);
    V2 = coda::var2iqlrContrast(s, qLow2, qHigh2);
    O.middleCols(n*P, P) = coda::quadForm(S, V1, V2);
  }
  return O;
}


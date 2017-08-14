// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
#include <Rcpp.h>
#include <cmath>

double integrand_o4 (double x, double nu, double t1, double t2,
            double delta1, double delta2) {
  double f = (R::pnorm(t2*x /sqrt(nu) - delta2, 0.0, 1.0, 1, 0) -
              R::pnorm(t1*x /sqrt(nu) - delta1, 0.0, 1.0, 1, 0))*
              exp((nu-1)*log(x) - x*x/2 - (nu/2 - 1) * log(2) - lgamma(nu/2));
  return f;
}

class Integrand_o4: public Func
{
private:
  double nu;
  double t1;
  double delta1;
  double t2;
  double delta2;
public:
  Integrand_o4 (double nu_, double t1_, double delta1_, double t2_, double delta2_) :
  nu(nu_), t1(t1_), delta1(delta1_), t2(t2_), delta2(delta2_) {}

  double operator()(const double& x) const
  {
    return integrand_o4(x, nu, t1, t2, delta1, delta2);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector ipowen4(double nu, double t1, double t2, double delta1, double delta2,
                           int subdiv=100, double eps_abs=1e-14, double eps_rel=1e-14){
  Integrand_o4 f(nu, t1, delta1, t2, delta2);
  const double R = (delta1-delta2)/(t1-t2)*sqrt(nu);
  double err_est;
  int err_code;
  const double res = integrate(f, 0, R, err_est, err_code, subdiv, eps_abs, eps_rel,
                               Integrator<double>::GaussKronrod201);
  Rcpp::NumericVector out =
    Rcpp::NumericVector::create(res);
  out.attr("err_est") = err_est;
  out.attr("err_code") = err_code;
  return out;
}


double integrand_Q (double x, double nu, double t, double delta) {
  double f = R::pnorm(t*x /sqrt(nu) - delta, 0.0, 1.0, 1, 0) *
              exp((nu-1)*log(x) - x*x/2 - (nu/2 - 1) * log(2) - lgamma(nu/2));
  return f;
}

class Integrand_Q1: public Func
{
private:
  double nu;
  double t;
  double delta;
public:
  Integrand_Q1 (double nu_, double t_, double delta_) :
  nu(nu_), t(t_), delta(delta_) {}

  double operator()(const double& x) const
  {
    return integrand_Q(x, nu, t, delta);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector iOwenQ1(double nu, double t, double delta, double R,
                            int subdiv=100, double eps_abs=1e-14, double eps_rel=1e-14){
  Integrand_Q1 f(nu, t, delta);
  double err_est;
  int err_code;
  const double res = integrate(f, 0, R, err_est, err_code, subdiv, eps_abs, eps_rel,
                               Integrator<double>::GaussKronrod201);
  Rcpp::NumericVector out =
    Rcpp::NumericVector::create(res);
  out.attr("err_est") = err_est;
  out.attr("err_code") = err_code;
  return out;
}

class Integrand_Q2: public Func
{
private:
  double nu;
  double t;
  double delta;
  double R;
public:
  Integrand_Q2 (double nu_, double t_, double delta_, double R_) :
  nu(nu_), t(t_), delta(delta_), R(R_) {}

  double operator()(const double& x) const
  {
    return integrand_Q(R + x/(1-x), nu, t, delta)/(1-x)/(1-x);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector iOwenQ2(double nu, double t, double delta, double R,
                            int subdiv=100, double eps_abs=1e-14, double eps_rel=1e-14){
  Integrand_Q2 f(nu, t, delta, R);
  double err_est;
  int err_code;
  const double res = integrate(f, 0, 1, err_est, err_code, subdiv, eps_abs, eps_rel,
                               Integrator<double>::GaussKronrod201);
  Rcpp::NumericVector out =
    Rcpp::NumericVector::create(res);
  out.attr("err_est") = err_est;
  out.attr("err_code") = err_code;
  return out;
}

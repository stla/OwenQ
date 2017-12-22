// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
#include <Rcpp.h>
#include <cmath>
#include <functional>


//---------- powen4 ---------//
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


//--------- OwenQ1 ------------//
double integrand_Q (double x, double nu, double t, double delta) {
  double f = exp(R::pnorm(t*x /sqrt(nu) - delta, 0.0, 1.0, 1, 1) +
              (nu-1)*log(x) - x*x/2 - (nu/2 - 1) * log(2) - lgamma(nu/2));
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

//------ OwenQ2 -------//
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

//------ powen 2 -------//
double integrand_o2 (double x, double nu, double t1, double t2,
            double delta1, double delta2) {
  double f = (R::pnorm(t2*x /sqrt(nu) - delta2, 0.0, 1.0, 1, 0) -
              R::pnorm(t1*x /sqrt(nu) - delta1, 0.0, 1.0, 1, 0))*
              exp((nu-1)*log(x) - x*x/2 - (nu/2 - 1) * log(2) - lgamma(nu/2));
  return f;
}

class Integrand_o2: public Func
{
private:
  double nu;
  double t1;
  double delta1;
  double t2;
  double delta2;
  double R;
public:
  Integrand_o2 (double nu_, double t1_, double delta1_, double t2_, double delta2_, double R_) :
  nu(nu_), t1(t1_), delta1(delta1_), t2(t2_), delta2(delta2_), R(R_) {}

  double operator()(const double& x) const
  {
    return integrand_o2(R + x/(1-x), nu, t1, t2, delta1, delta2)/(1-x)/(1-x);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector ipowen2(double nu, double t1, double t2, double delta1, double delta2,
                           int subdiv=100, double eps_abs=1e-14, double eps_rel=1e-14){
  const double R = (delta1-delta2)/(t1-t2)*sqrt(nu);
  Integrand_o2 f(nu, t1, delta1, t2, delta2, R);
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

// // test
// class Integrand: public Func
// {
// private:
//   std::function<double(double)> f;
// public:
//   Integrand (std::function<double(double)> f_) :
//   f(f_) {}
//
//   double operator()(const double& x) const
//   {
//     return f(x);
//   }
// };
//
// double SEXPtoReal (SEXP vectorR){
//   double* vector = REAL(vectorR);
//   return vector[0];
// }
//
// Rcpp::NumericVector integral(Rcpp::Function f,
//                             int subdiv=100, double eps_abs=1e-14, double eps_rel=1e-14){
//   std::function<double(double)> ff = [&](double x){ return REAL(f(x))[0]; };
//   Integrand I(ff);
//   double err_est;
//   int err_code;
//   const double res = integrate(I, 0, 1, err_est, err_code, subdiv, eps_abs, eps_rel,
//                                Integrator<double>::GaussKronrod201);
//   Rcpp::NumericVector out =
//     Rcpp::NumericVector::create(res);
//   out.attr("err_est") = err_est;
//   out.attr("err_code") = err_code;
//   return out;
// }

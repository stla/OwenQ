#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <cfloat>

// [[Rcpp::export]]
double RcppOwenT01(double h, double a, int jmax, double cutpoint){
  const double twopi = 2*3.14159265358979323846;
  if(h > cutpoint){
    return atan(a) * exp(-0.5 * (h*h) * a / atan(a)) *
      (1 + 0.00868 * pow(h*a,4)) / twopi;
  }
  double cumsum = 1.0;
  const double K = exp(-h*h/2);
  double crossprod = (1-K)*a;
  int i;
  for(i=1; i<=jmax; i++){
    cumsum += exp(2*i*log(h) - i*log(2.0) - lgamma(i+1));
    double jk = 1.0/(2.0*i+1);
    if(i%2==1){
      jk = -jk;
    }
    crossprod += (1.0 - K*cumsum)*jk*pow(a,2*i+1);
  }
  return (atan(a)-crossprod)/twopi;
}
double RcppOwenT01(double h, double a, int jmax, double cutpoint);

//****************************************************************************80
double pNorm(double q){
  return R::pnorm(q,0.0,1.0,1,0);
}
double pNorm(double q);

//****************************************************************************80
// [[Rcpp::export]]
double RcppOwenT(double h, double a, int jmax, double cutpoint){
  if(a == 0.0){
    return 0.0;
  }
  double absa = fabs(a);
  double absh = fabs(h);
  if(absa > DBL_MAX){
    return R::sign(a) * pNorm(-absh) / 2;
  }
  if(absa <= 1.0){
    return R::sign(a) * RcppOwenT01(absh, absa, jmax, cutpoint);
  }
  return R::sign(a) * (pNorm(absh)/2 + pNorm(absa*absh) * (0.5 - pNorm(absh)) -
               RcppOwenT01(absa*absh, 1/absa, jmax, cutpoint));
}
double RcppOwenT(double h, double a, int jmax, double cutpoint);

//****************************************************************************80
NumericVector isPositive(NumericVector x){
  int n = x.size();
  NumericVector out(n);
  int i;
  for(i=0; i<n; i++){
    out[i] = x[i] >= 0;
  }
  return out;
}
NumericVector isPositive(NumericVector x);

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenQ1(int nu, double t, NumericVector delta, NumericVector R,
    int jmax, double cutpoint){
  const double a = R::sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sB = sqrt(b);
  double ab;
  double asB;
  if(fabs(t)>DBL_MAX){
    ab = 0;
    asB = sign(t);
  }else{
    ab = a*b;
    asB = R::sign(t)*sqrt(t*t/(nu+t*t));
  }
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(R) - isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      double C1 =
        RcppOwenT(delta[i]*sB, a, jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M(n,J);
  H(0,_) = -dnorm(R) * pnorm(a*R-delta);
  M(0,_) = asB*dnorm(delta*sB)*(pnorm(delta*asB)-pnorm((delta*ab-R)/sB));
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M(1,_) = delta*ab*M(0,_) + ab*dnorm(delta*sB)*(dnorm(delta*asB)-dnorm((delta*ab-R)/sB));
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L(0,_) = ab * R * dnorm(R) * dnorm(a*R-delta) / 2;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L(k,_) = A[k+2] * R * L(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        H(k,_) = A[k] * R * H(k-1,_);
        M(k,_) = (k-1.0)/k * (A[k-2] * delta * ab * M(k-1,_) + b*M(k-2,_)) - L(k-2,_);
      }
    }
  }
  if(nu % 2 == 0){
    double sqrt2pi = 2.506628274631000502415765284811;
    NumericVector sum(J);
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M(i,_)+H(i,_);
    }
    return pnorm(-delta) + sqrt2pi * sum;
  }else{
    NumericVector sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_)+H(i,_);
    }
    NumericVector C = pnorm(R) - isPositive(delta);
    for(i=0; i<J; i++){
      double C1 =
        RcppOwenT(delta[i]*sB, a, jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C+2*sum;
  }
}

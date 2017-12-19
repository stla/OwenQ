// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <cfloat>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/constants/constants.hpp>

const double sqrt2pi = boost::math::constants::root_two_pi<double>();
const double onedivsqrt2pi = 
             boost::math::constants::one_div_root_two_pi<double>();
const double logsqrt2pi = boost::math::constants::log_root_two_pi<double>();
const double logtwo = 0.693147180559945309417232121458176568;

NumericVector xdnormx(NumericVector x){
  return exp(log(x) - 0.5*x*x - logsqrt2pi);
}

//NumericVector xndnormx(NumericVector x, int n){
//  return exp(n*log(x) - 0.5*x*x - logsqrt2pi);
//}
//
//double logAA(int k){
//  if(k == 0){
//    return 0.0;
//  }
//  if(k % 2 == 1){
//    return (k-1.0)*logtwo + 2*lgamma((k+1.0)/2.0) - lgamma(k+1.0);
//  }else{
//    return lgamma(k) - (k-2)*logtwo - 2.0*lgamma(k/2.0) - log(k);
//  }
//}
//
//double AA(int k){
//  return exp(logAA(k));
//}


// [[Rcpp::export]]
double RcppOwenT(double h, double a){
  return boost::math::owens_t(h, a);
}
double RcppOwenT(double h, double a);

// //****************************************************************************80
// double pNorm(double q){
//   return R::pnorm(q,0.0,1.0,1,0);
// }
// double pNorm(double q);

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenStudent(double q, int nu, NumericVector delta){
  const double a = R::sign(q)*sqrt(q*q/nu);
  const double b = nu/(nu+q*q);
  NumericVector dsB = delta*sqrt(b);
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(-dsB);
    int i;
    for(i=0; i<J; i++){
      C[i] += 2*RcppOwenT(dsB[i], a);
    }
    return C;
  }
  NumericMatrix M(nu-1,J);
  M(0,_) = a * sqrt(b) * dnorm(dsB) * pnorm(a*dsB);
  if(nu>2){
    M(1,_) = b * (delta * a * M(0,_) + a * dnorm(delta) * onedivsqrt2pi);
    if(nu>3){
      NumericVector A(nu-3);
      A[0] = 1;
      int k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/k/A[k-1];
        }
      }
      for(k=2; k<nu-1; k++){
        M(k,_) = (k-1) * b * (A[k-2] * delta * a * M(k-1,_) + M(k-2L,_)) / k;
      }
    }
  }
  int i;
  if(nu%2==1){
    NumericVector C = pnorm(-dsB);
    for(i=0; i<J; i++){
      C[i] += 2*RcppOwenT(dsB[i], a);
    }
    NumericVector sum(J);
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_);
    }
    return C + 2*sum;
  }
  NumericVector sum(J);
  for(i=0; i<nu-1; i+=2){
    sum += M(i,_);
  }
  return pnorm(-delta) + sqrt2pi*sum;
}

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
NumericVector RcppOwenQ1
    (int nu, double t, NumericVector delta, NumericVector R, int algo=1){
  const double a = R::sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sB = sqrt(b);
  const double ab = sqrt(nu)/(nu/t + t);
  const double asB = R::sign(t)*1/sqrt(nu/(t*t)+1);
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(R) - isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta[i]*sB, a);
      double C2 = RcppOwenT(R[i], a-delta[i]/R[i]);
      double C3 = RcppOwenT(delta[i]*sB, (ab-R[i]/delta[i])/b);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M(n,J);
  NumericVector Hfactor = -pnorm(a*R-delta); 
  NumericVector Lfactor = ab * dnorm(a*R-delta);
  H(0,_) = dnorm(R);
  M(0,_) = asB*dnorm(delta*sB)*(pnorm(delta*asB)-pnorm((delta*ab-R)/sB));
  if(nu >= 3){
    H(1,_) = xdnormx(R);
    M(1,_) = delta*ab*M(0,_) + 
                ab*dnorm(delta*sB)*(dnorm(delta*asB)-dnorm((delta*ab-R)/sB));
    if(nu >= 4){
      int k;
      if(algo==1){
        NumericVector A(n);
 //       NumericVector prodA(n);
        NumericMatrix L(n-2,J);
        A[0] = 1.0;
        A[1] = 1.0;
 //       prodA[0] = 1.0;
 //       prodA[1] = 1.0;
        L(0,_) = 0.5 * H(1,_);
        for(k=2; k<n; k++){
          A[k] = 1.0/k/A[k-1];
 //         prodA[k] = A[k] * prodA[k-1];
        }
        if(nu >= 5){
          for(k=1; k<n-2; k++){
            L(k,_) = A[k+2] * R * L(k-1,_); 
//              L(k,_) = prodA[k+2] * xndnormx(R,k+1);
          }
        }
        for(k=2; k<n; k++){
//            H(k,_) = prodA[k] * xndnormx(R,k);
            H(k,_) = R * H(k-1,_);
            M(k,_) = (k-1.0)/k * (A[k-2] * delta * ab * M(k-1,_) + b*M(k-2,_)) - 
                        Lfactor*L(k-2,_);
        }
      }else{ // algo 2
        NumericVector A(n-1);
        A[0] = 1.0;
        NumericVector halfRR = 0.5*R*R;
        NumericVector logR = log(R);
        for(k=0; k<n-2; k++){
          A[k+1] = 1.0/(k+1.0)/A[k]; // un de trop
          //double Ak = AA(k);
          double ldf;
          if(k % 2 == 0){
            ldf = (0.5*k+1.0)*logtwo + lgamma(2.0+0.5*k);
          }else{
            ldf = lgamma(k+3.0) - 0.5*(k+1.0)*logtwo - lgamma(0.5*(k+3.0));
          }
          double r = (k+1.0)/(k+2.0);
          NumericVector K = 
                  exp(-ldf + (k+1.0)*logR - halfRR - logsqrt2pi);
          H(k+2,_) = K*R;
          M(k+2,_) = r*(A[k]*delta*ab*M(k+1,_)+b*M(k,_)) - K*Lfactor;
        }          
      }
    }
  }
  NumericVector sumM(J);
  NumericVector sumH(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sumM += M(i,_);
      sumH += H(i,_);
    }
    return pnorm(-delta) + sqrt2pi * (sumM + Hfactor*sumH);
  }else{
    for(i=1; i<nu-1; i+=2){
      sumM += M(i,_);
      sumH += H(i,_);
    }
    NumericVector C = pnorm(R) - isPositive(delta);
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta[i]*sB, a);
      double C2 = RcppOwenT(R[i], a-delta[i]/R[i]);
      double C3 = RcppOwenT(delta[i]*sB, (ab-R[i]/delta[i])/b);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C + 2*(sumM + Hfactor*sumH);
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenQ2(int nu, double t, NumericVector delta, NumericVector R){
  const double a = R::sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sB = sqrt(b);
  double ab;
  double asB;
  if(fabs(t)>DBL_MAX){
    ab = 0;
    asB = R::sign(t);
  }else{
    ab = a*b;
    asB = R::sign(t)*sqrt(t*t/(nu+t*t));
  }
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(-delta*sB) - pnorm(R) + isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      double C2 = RcppOwenT(R[i], (a*R[i]-delta[i])/R[i]);
      double C3 = RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i]);
      C[i] += 2*(C2 + C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M(n,J);
  H(0,_) = -dnorm(R) * pnorm(a*R-delta);
  M(0,_) = asB*dnorm(delta*sB)*pnorm((delta*ab-R)/sB);
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M(1,_) = delta*ab*M(0,_) + ab*dnorm(delta*sB)*dnorm((delta*ab-R)/sB);
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
        M(k,_) = (k-1.0)/k * (A[k-2] * delta * ab * M(k-1,_) + b*M(k-2,_)) + L(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    return sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    NumericVector C = pnorm(-delta*sB) - pnorm(R) + isPositive(delta);
    for(i=0; i<J; i++){
      double C2 = RcppOwenT(R[i], (a*R[i]-delta[i])/R[i]);
      double C3 = RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i]);
      C[i] += 2*(C2 + C3);
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF4(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2){
  const double a1 = R::sign(t1)*sqrt(t1*t1/nu);
  const double b1 = nu/(nu+t1*t1);
  const double sB1 = sqrt(b1);
  const double ab1 = a1*b1;
  const double asB1 = R::sign(t1)*sqrt(t1*t1/(nu+t1*t1));
  const double a2 = R::sign(t2)*sqrt(t2*t2/nu);
  const double b2 = nu/(nu+t2*t2);
  const double sB2 = sqrt(b2);
  const double ab2 = a2*b2;
  const double asB2 = R::sign(t2)*sqrt(t2*t2/(nu+t2*t2));
  const int J = delta1.size();
  const NumericVector R = sqrt(nu)*(delta1 - delta2)/(t1-t2);
  if(nu==1){
    NumericVector C = isPositive(delta1) - isPositive(delta2);
    int i;
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta2[i]*sB2, a2) - RcppOwenT(delta1[i]*sB1, a1);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
                    RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M1(n,J);
  NumericMatrix M2(n,J);
  H(0,_) = -dnorm(R) * (pnorm(a2*R-delta2) - pnorm(a1*R-delta1));
  M1(0,_) = asB1*dnorm(delta1*sB1)*(pnorm(delta1*asB1)-pnorm((delta1*ab1-R)/sB1));
  M2(0,_) = asB2*dnorm(delta2*sB2)*(pnorm(delta2*asB2)-pnorm((delta2*ab2-R)/sB2));
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M1(1,_) = delta1*ab1*M1(0,_) +
              ab1*dnorm(delta1*sB1)*(dnorm(delta1*asB1)-dnorm((delta1*ab1-R)/sB1));
    M2(1,_) = delta2*ab2*M2(0,_) +
              ab2*dnorm(delta2*sB2)*(dnorm(delta2*asB2)-dnorm((delta2*ab2-R)/sB2));
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L1(n-2,J);
      NumericMatrix L2(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L1(0,_) = ab1 * R * dnorm(R) * dnorm(a1*R-delta1) / 2;
      L2(0,_) = ab2 * R * dnorm(R) * dnorm(a2*R-delta2) / 2;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L1(k,_) = A[k+2] * R * L1(k-1,_);
          L2(k,_) = A[k+2] * R * L2(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        H(k,_) = A[k] * R * H(k-1,_);
        M1(k,_) = (k-1.0)/k *
                  (A[k-2] * delta1 * ab1 * M1(k-1,_) + b1*M1(k-2,_)) - L1(k-2,_);
        M2(k,_) = (k-1.0)/k *
                  (A[k-2] * delta2 * ab2 * M2(k-1,_) + b2*M2(k-2,_)) - L2(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += M2(i,_)-M1(i,_)+H(i,_);
    }
    return pnorm(-delta2) - pnorm(-delta1) + sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += M2(i,_)-M1(i,_)+H(i,_);
    }
    NumericVector C = isPositive(delta1) - isPositive(delta2);
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta2[i]*sB2, a2) - RcppOwenT(delta1[i]*sB1, a1);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
                    RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF3(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2){
  const double a1 = R::sign(t1)*sqrt(t1*t1/nu);
  const double b1 = nu/(nu+t1*t1);
  const double sB1 = sqrt(b1);
  const double ab1 = a1*b1;
  const double asB1 = R::sign(t1)*sqrt(t1*t1/(nu+t1*t1));
  const double a2 = R::sign(t2)*sqrt(t2*t2/nu);
  const double b2 = nu/(nu+t2*t2);
  const double sB2 = sqrt(b2);
  const double ab2 = a2*b2;
  const double asB2 = R::sign(t2)*sqrt(t2*t2/(nu+t2*t2));
  const int J = delta1.size();
  const NumericVector R = sqrt(nu)*(delta1 - delta2)/(t1-t2);
  if(nu==1){
    NumericVector C = - isPositive(delta1) + isPositive(delta2) - pnorm(-delta1*sB1);
    int i;
    for(i=0; i<J; i++){
      double C1 = - RcppOwenT(delta2[i]*sB2, a2);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
                    RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 + C2 + C3) + 1.0;
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M1(n,J);
  NumericMatrix M2(n,J);
  H(0,_) = -dnorm(R) * (pnorm(a2*R-delta2) - pnorm(a1*R-delta1));
  M1(0,_) = asB1*dnorm(delta1*sB1)*pnorm((delta1*ab1-R)/sB1);
  M2(0,_) = asB2*dnorm(delta2*sB2)*(pnorm(delta2*asB2)-pnorm((delta2*ab2-R)/sB2));
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M1(1,_) = delta1*ab1*M1(0,_) + ab1*dnorm(delta1*sB1)*dnorm((delta1*ab1-R)/sB1);
    M2(1,_) = delta2*ab2*M2(0,_) +
              ab2*dnorm(delta2*sB2)*(dnorm(delta2*asB2)-dnorm((delta2*ab2-R)/sB2));
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L1(n-2,J);
      NumericMatrix L2(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L1(0,_) = 0.5 * ab1 * R * dnorm(R) * dnorm(a1*R-delta1);
      L2(0,_) = 0.5 * ab2 * R * dnorm(R) * dnorm(a2*R-delta2);
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L1(k,_) = A[k+2] * R * L1(k-1,_);
          L2(k,_) = A[k+2] * R * L2(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        H(k,_) = A[k] * R * H(k-1,_);
        M1(k,_) = (k-1.0)/k *
                  (A[k-2] * delta1 * ab1 * M1(k-1,_) + b1*M1(k-2,_)) + L1(k-2,_);
        M2(k,_) = (k-1.0)/k *
                  (A[k-2] * delta2 * ab2 * M2(k-1,_) + b2*M2(k-2,_)) - L2(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += -M2(i,_)-M1(i,_)+H(i,_);
    }
    return 1.0 - pnorm(-delta2) + sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += -M2(i,_)-M1(i,_)+H(i,_);
    }
    NumericVector C = - isPositive(delta1) + isPositive(delta2) - pnorm(-delta1*sB1);
    for(i=0; i<J; i++){
      double C1 = - RcppOwenT(delta2[i]*sB2, a2);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
                    RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 + C2 + C3) + 1.0;
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF2(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2){
  const double a1 = R::sign(t1)*sqrt(t1*t1/nu);
  const double b1 = nu/(nu+t1*t1);
  const double sB1 = sqrt(b1);
  const double ab1 = a1*b1;
  const double asB1 = R::sign(t1)*sqrt(t1*t1/(nu+t1*t1));
  const double a2 = R::sign(t2)*sqrt(t2*t2/nu);
  const double b2 = nu/(nu+t2*t2);
  const double sB2 = sqrt(b2);
  const double ab2 = a2*b2;
  const double asB2 = R::sign(t2)*sqrt(t2*t2/(nu+t2*t2));
  const int J = delta1.size();
  const NumericVector R = sqrt(nu)*(delta1 - delta2)/(t1-t2);
  if(nu==1){
    NumericVector C = isPositive(delta1) - isPositive(delta2) +
                        pnorm(-delta1*sB1) - pnorm(-delta2*sB2);
    int i;
    for(i=0; i<J; i++){
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
                    RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += -2*(C2 + C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M1(n,J);
  NumericMatrix M2(n,J);
  H(0,_) = -dnorm(R) * (pnorm(a2*R-delta2) - pnorm(a1*R-delta1));
  M1(0,_) = asB1*dnorm(delta1*sB1)*pnorm((delta1*ab1-R)/sB1);
  M2(0,_) = asB2*dnorm(delta2*sB2)*pnorm((delta2*ab2-R)/sB2);
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M1(1,_) = delta1*ab1*M1(0,_) + ab1*dnorm(delta1*sB1)*dnorm((delta1*ab1-R)/sB1);
    M2(1,_) = delta2*ab2*M2(0,_) + ab2*dnorm(delta2*sB2)*dnorm((delta2*ab2-R)/sB2);
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L1(n-2,J);
      NumericMatrix L2(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L1(0,_) = 0.5 * ab1 * R * dnorm(R) * dnorm(a1*R-delta1);
      L2(0,_) = 0.5 * ab2 * R * dnorm(R) * dnorm(a2*R-delta2);
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L1(k,_) = A[k+2] * R * L1(k-1,_);
          L2(k,_) = A[k+2] * R * L2(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        H(k,_) = A[k] * R * H(k-1,_);
        M1(k,_) = (k-1.0)/k *
                  (A[k-2] * delta1 * ab1 * M1(k-1,_) + b1*M1(k-2,_)) + L1(k-2,_);
        M2(k,_) = (k-1.0)/k *
                  (A[k-2] * delta2 * ab2 * M2(k-1,_) + b2*M2(k-2,_)) + L2(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += -M2(i,_)+M1(i,_)-H(i,_);
    }
    return sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += -M2(i,_)+M1(i,_)-H(i,_);
    }
    NumericVector C = isPositive(delta1) - isPositive(delta2) +
                          pnorm(-delta1*sB1) - pnorm(-delta2*sB2);
    for(i=0; i<J; i++){
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += -2*(C2 + C3);
    }
    return C+2*sum;
  }
}


//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF1(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2){
  const double a1 = R::sign(t1)*sqrt(t1*t1/nu);
  const double b1 = nu/(nu+t1*t1);
  const double sB1 = sqrt(b1);
  const double ab1 = a1*b1;
  const double asB1 = R::sign(t1)*sqrt(t1*t1/(nu+t1*t1));
  const double a2 = R::sign(t2)*sqrt(t2*t2/nu);
  const double b2 = nu/(nu+t2*t2);
  const double sB2 = sqrt(b2);
  const double ab2 = a2*b2;
  const double asB2 = R::sign(t2)*sqrt(t2*t2/(nu+t2*t2));
  const int J = delta1.size();
  const NumericVector R = sqrt(nu)*(delta1 - delta2)/(t1-t2);
  if(nu==1){
    NumericVector C = -isPositive(delta1) + isPositive(delta2) + 
                        pnorm(-delta2*sB2);
    int i;
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta1[i]*sB1, a1);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 + C2 + C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix M1(n,J);
  NumericMatrix M2(n,J);
  M1(0,_) = asB1*dnorm(delta1*sB1) * 
                (pnorm(delta1*asB1)-pnorm((delta1*ab1-R)/sB1));
  M2(0,_) = asB2*dnorm(delta2*sB2)*pnorm((delta2*ab2-R)/sB2);
  if(nu >= 3){
    M1(1,_) = delta1*ab1*M1(0,_) +
              ab1*dnorm(delta1*sB1)*(dnorm(delta1*asB1)-dnorm((delta1*ab1-R)/sB1));
    M2(1,_) = delta2*ab2*M2(0,_) +
              ab2*dnorm(delta2*sB2)*dnorm((delta2*ab2-R)/sB2);
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L1(n-2,J);
      NumericMatrix L2(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L1(0,_) = 0.5 * ab1 * R * dnorm(R) * dnorm(a1*R-delta1);
      L2(0,_) = 0.5 * ab2 * R * dnorm(R) * dnorm(a2*R-delta2);
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L1(k,_) = A[k+2] * R * L1(k-1,_);
          L2(k,_) = A[k+2] * R * L2(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        M1(k,_) = (k-1.0)/k *
                  (A[k-2] * delta1 * ab1 * M1(k-1,_) + b1*M1(k-2,_)) - L1(k-2,_);
        M2(k,_) = (k-1.0)/k *
                  (A[k-2] * delta2 * ab2 * M2(k-1,_) + b2*M2(k-2,_)) + L2(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += M2(i,_)+M1(i,_);
    }
    return pnorm(-delta1) + sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += M2(i,_)+M1(i,_);
    }
    NumericVector C = -isPositive(delta1) + isPositive(delta2) + pnorm(-delta2*sB2);
    for(i=0; i<J; i++){
      double C1 = RcppOwenT(delta1[i]*sB1, a1);
      double C2 = RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i]) -
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i]);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i]) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i]);
      C[i] += 2*(C1 + C2 + C3);
    }
    return C+2*sum;
  }
}


//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppSpecialOwenCDF2(int nu, double t, NumericVector delta){
  const double a = sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sB = sqrt(b);
  const double ab = a*b;
  const double asB = sqrt(t*t/(nu+t*t));
  const int J = delta.size();
  const NumericVector R = sqrt(nu)*delta/t;
  NumericVector C(J);
  if(nu % 2 == 1){
    C = 2*pnorm(-delta*sB);
    for(int i=0; i<J; i++){
      double C2 = 2*RcppOwenT(R[i], (a*R[i]-delta[i])/R[i]);
      double C3 =
        2*RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i]);
      C[i] += 2*(C2 + C3);
    }
    if(nu == 1){
      return C;
    }
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M(n,J);
  H(0,_) = -dnorm(R) * (1-2*pnorm(a*R-delta));
  M(0,_) = 2*asB*dnorm(delta*sB)*pnorm((delta*ab-R)/sB);
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M(1,_) = delta*ab*M(0,_) + 2*ab*dnorm(delta*sB)*dnorm((delta*ab-R)/sB);
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L(0,_) = ab * R * dnorm(R) * dnorm(a*R-delta);
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
        M(k,_) = (k-1.0)/k *
                  (A[k-2] * delta * ab * M(k-1,_) + b*M(k-2,_)) + L(k-2,_);
      }
    }
  }
  NumericVector sum(J);
  int i;
  if(nu % 2 == 0){
    for(i=0; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    return sqrt2pi * sum;
  }else{
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    return C+2*sum;
  }
}

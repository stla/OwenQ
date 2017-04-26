#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <cfloat>

// [[Rcpp::export]]
double RcppOwenT01(double h, double a, int jmax=50, double cutpoint=8){
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
double RcppOwenT(double h, double a, int jmax=50, double cutpoint=8){
  if(a == 0.0){
    return 0.0;
  }
  double absh = fabs(h);
  if(absh > DBL_MAX){
    return 0.0;
  }
  double absa = fabs(a);
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
// [[Rcpp::export]]
NumericVector RcppOwenStudent(double q, int nu, NumericVector delta,
    int jmax=50, double cutpoint=8){
  const double a = R::sign(q)*sqrt(q*q/nu);
  const double b = nu/(nu+q*q);
  NumericVector dsB = delta*sqrt(b);
  const int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(-dsB);
    int i;
    for(i=0; i<J; i++){
      C[i] += 2*RcppOwenT(dsB[i], a, jmax, cutpoint);
    }
    return C;
  }
  NumericMatrix M(nu-1,J);
  M(0,_) = a * sqrt(b) * dnorm(dsB) * pnorm(a*dsB);
  const double sqrt2pi = 2.506628274631000502415765284811;
  if(nu>2){
    M(1,_) = b * (delta * a * M(0,_) + a * dnorm(delta) / sqrt2pi);
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
      C[i] += 2*RcppOwenT(dsB[i], a, jmax, cutpoint);
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
NumericVector RcppOwenQ1(int nu, double t, NumericVector delta, NumericVector R,
    int jmax=50, double cutpoint=8){
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
    NumericVector C = pnorm(R) - isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      double C1 =
        RcppOwenT(delta[i]*sB, a, jmax, cutpoint);
      // double C2 =
      //   RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], a-delta[i]/R[i], jmax, cutpoint);
      double C3 =
        // RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
        RcppOwenT(delta[i]*sB, (ab-R[i]/delta[i])/b, jmax, cutpoint);
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
    const double sqrt2pi = 2.506628274631000502415765284811;
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
      // double C2 =
      //   RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], a-delta[i]/R[i], jmax, cutpoint);
      double C3 =
        // RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
        RcppOwenT(delta[i]*sB, (ab-R[i]/delta[i])/b, jmax, cutpoint);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenQ2(int nu, double t, NumericVector delta, NumericVector R,
    int jmax=50, double cutpoint=8){
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
      double C2 =
        RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
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
  if(nu % 2 == 0){
    const double sqrt2pi = 2.506628274631000502415765284811;
    NumericVector sum(J);
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    return sqrt2pi * sum;
  }else{
    NumericVector sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_)-H(i,_);
    }
    NumericVector C = pnorm(-delta*sB) - pnorm(R) + isPositive(delta);
    for(i=0; i<J; i++){
      double C2 =
        RcppOwenT(R[i], (a*R[i]-delta[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax, cutpoint);
      C[i] += 2*(C2 + C3);
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF4(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2, int jmax=50, double cutpoint=8){
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
      double C1 =
        RcppOwenT(delta2[i]*sB2, a2, jmax, cutpoint) -
          RcppOwenT(delta1[i]*sB1, a1, jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i], jmax, cutpoint) -
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i], jmax, cutpoint) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i], jmax, cutpoint);
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
  if(nu % 2 == 0){
    const double sqrt2pi = 2.506628274631000502415765284811;
    NumericVector sum(J);
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M2(i,_)-M1(i,_)+H(i,_);
    }
    return pnorm(-delta2) - pnorm(-delta1) + sqrt2pi * sum;
  }else{
    NumericVector sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M2(i,_)-M1(i,_)+H(i,_);
    }
    NumericVector C = isPositive(delta1) - isPositive(delta2);
    for(i=0; i<J; i++){
      double C1 =
        RcppOwenT(delta2[i]*sB2, a2, jmax, cutpoint) -
          RcppOwenT(delta1[i]*sB1, a1, jmax, cutpoint);
      double C2 =
        RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i], jmax, cutpoint) -
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i], jmax, cutpoint);
      double C3 =
        RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i], jmax, cutpoint) -
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i], jmax, cutpoint);
      C[i] += 2*(C1 - C2 - C3);
    }
    return C+2*sum;
  }
}

//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenCDF3(int nu, double t1, double t2, NumericVector delta1,
    NumericVector delta2, int jmax=50, double cutpoint=8){
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
    NumericVector C = 1 - isPositive(delta1) - isPositive(delta2) - pnorm(-delta*sB);
    int i;
    for(i=0; i<J; i++){
      double C1 =
        - RcppOwenT(delta2[i]*sB2, a2, jmax, cutpoint);
      double C2 =
        - RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i], jmax, cutpoint) +
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i], jmax, cutpoint);
      double C3 =
        - RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i], jmax, cutpoint) +
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i], jmax, cutpoint);
      C[i] += 2*(C1 + C2 + C3);
    }
    return C;
  }
  const int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M1(n,J);
  NumericMatrix M2(n,J);
  H(0,_) = -dnorm(R) * (pnorm(a2*R-delta2) + pnorm(a1*R-delta1));
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
  if(nu % 2 == 0){
    const double sqrt2pi = 2.506628274631000502415765284811;
    NumericVector sum(J);
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += -M2(i,_)-M1(i,_)+H(i,_);
    }
    return 1.0 - pnorm(-delta2) + sqrt2pi * sum;
  }else{
    NumericVector sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += -M2(i,_)-M1(i,_)+H(i,_);
    }
    NumericVector C = 1 - isPositive(delta1) - isPositive(delta2) - pnorm(-delta*sB);
    int i;
    for(i=0; i<J; i++){
      double C1 =
        - RcppOwenT(delta2[i]*sB2, a2, jmax, cutpoint);
      double C2 =
        - RcppOwenT(R[i], (a2*R[i]-delta2[i])/R[i], jmax, cutpoint) +
          RcppOwenT(R[i], (a1*R[i]-delta1[i])/R[i], jmax, cutpoint);
      double C3 =
        - RcppOwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i], jmax, cutpoint) +
          RcppOwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i], jmax, cutpoint);
      C[i] += 2*(C1 + C2 + C3);
    }
    return C+2*sum;
  }
}


// //****************************************************************************80
// //
// NumericVector Cconstant(double t, double delta, double R,
//     int jmax=50, double cutpoint=8){
//   int nu = 1;
//   const double a = R::sign(t)*sqrt(t*t/nu);
//   const double b = nu/(nu+t*t);
//   const double sB = sqrt(b);
//   double ab = a*b;
//   double C = pNorm(R) - (delta > 0);
//   double C1 =
//     RcppOwenT(delta*sB, a, jmax, cutpoint);
//   double C2 =
//     RcppOwenT(R, (a*R-delta)/R, jmax, cutpoint);
//   double C3 =
//     RcppOwenT(delta*sB, (delta*ab-R)/b/delta, jmax, cutpoint);
//   C += 2*(C1 - C2 - C3);
//   return NumericVector::create(C, C1, C2, C3, R, (a*R-delta)/R);
// }

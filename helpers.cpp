// Copyright (c) 2026 Dennis Dobler, Paavo Sattler, Nils Hichert, 
//                    Jörg-Tobias Kuhn, Lubna Amro
// Licensed under the MIT License. See LICENSE file for details.

//##############################################################################
//-----------------------------------------------------------------------------#
//--------------- This file contains all C++ helper functions -----------------#
//-----------------------------------------------------------------------------#
//##############################################################################

#include <Rcpp.h>
using namespace Rcpp;

inline double cfun(double x){
  if(x>0) return 1.0;
  if(x<0) return 0.0;
  return 0.5;
}

// [[Rcpp::export]]
NumericVector est_p_cpp(NumericMatrix X, NumericMatrix L) {
  
  int d = X.nrow(), n = X.ncol();
  double sum_k, sum_s, sum_l, diff;
  
  // pmax(rowSums(L), 1)
  NumericVector l_(d);
  for (int i = 0; i < d; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += L(i, j);
    }
    l_[i] = sum;
  }
  
  // calculate p 
  NumericVector p(d);
  for (int i = 0; i < d; i++) {
    if(l_[i] == 0) continue; // because then the summand is 0
    sum_k = 0.0;
    for (int k = 0; k < n; k++) {
      sum_s = 0.0;
      for (int s = 0; s < d; s++) {
        if(l_[s] == 0) continue; // because then the summand is 0
        sum_l = 0.0;
        for (int l = 0; l < n; l++) {
          diff = X(i, k) - X(s, l);
          sum_l += L(s, l) * cfun(diff);
        }
        sum_s += sum_l / l_[s];
      }
      sum_k += sum_s * L(i, k) / d;
    }
    p[i] = sum_k / l_[i];
  }
  return p;
}


// [[Rcpp::export]]
NumericMatrix est_V_cpp(NumericMatrix X, NumericMatrix L) {
  
  int d = X.nrow(), n = X.ncol();
  
  // lambda_i.
  NumericVector l_(d);
  for(int i=0;i<d;i++){
    double s=0;
    for(int j=0;j<n;j++) s+=L(i,j);
    l_[i]=s;
  }
  
  // pre-compute all weights W_1^ils = lambda_il/(lambda_s.*lambda_i.)
  // and store in matrix W_1
  NumericVector W1(d*d*n);
  for(int i=0;i<d;i++){
    for(int s=0;s<d;s++){
      for(int l=0;l<n;l++){
        W1[i + d*s + d*d*l] = L(i,l)/(l_(s)*l_(i));
      }
    }
  }
  
  
  
  // pre-compute all weights W_2^ij = n*lambda_ij/lambda_i.
  // and store in matrix W_2
  NumericMatrix W2(d,n);
  for(int s=0;s<d;s++){
    for(int j=0;j<n;j++){
      W2(s,j) = n*L(s,j)/l_(s);
    }
  }
  
  
  // pre-compute inner inner sum that does not depend on k
  NumericVector R(d*d*n);
  for(int i=0;i<d;i++){
    for(int s=0;s<d;s++){
      for(int l=0;l<n;l++){
        double sum=0;
        for(int j=0;j<n;j++){
          sum += cfun(X(i,l)-X(s,j))*W2(s,j);
        }
        R[i + d*s + d*d*l] = sum;
      }
    }
  }
  
  
  
  // pre-compute inner sum that does not depend on l
  NumericVector H(d*d*n);
  for(int i=0;i<d;i++){
    for(int s=0;s<d;s++){
      for(int k=0;k<n;k++){
        double sum=0;
        for(int l=0;l<n;l++){
          int idx = i + d*s + d*d*l;
          sum += W1[idx] * (n*cfun(X(i,l)-X(s,k)) - R[idx]);
        }
        H[i + d*s + d*d*k] = sum;
      }
    }
  }
  
  // Psi
  NumericMatrix Psi(d,n);
  for(int k=0;k<n;k++){
    for(int i=0;i<d;i++){
      double l_ik = L(i,k);
      double sum = 0;
      for(int s=0;s<d;s++){
        sum += L(s,k) * H[i + d*s + d*d*k] - l_ik * H[s + d*i + d*d*k];
      }
      Psi(i,k) = sum/d;
    }
  }
  
  // V
  NumericMatrix V(d,d);
  for(int i=0;i<d;i++){
    for(int j=i;j<d;j++){
      double sum=0;
      for(int k=0;k<n;k++)
        sum += Psi(i,k)*Psi(j,k);
      V(i,j)=sum/n;
      V(j,i)=V(i,j);
    }
  }
  
  return V;
}

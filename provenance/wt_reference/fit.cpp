#include <Rcpp.h>
using namespace Rcpp;
// Exact weighted least-squares condition coefficient after eliminating source intercepts.
// [[Rcpp::export]]
NumericVector weighted_condition(NumericMatrix y,NumericMatrix w,IntegerVector group,IntegerVector cond) {
 int n=y.nrow(),m=y.ncol(),G=max(group);NumericVector out(n);
 for(int i=0;i<n;i++) {
  double sw[6]={0},sx[6]={0},sy[6]={0},sxy[6]={0};
  for(int j=0;j<m;j++) {int g=group[j]-1;double a=w(i,j);sw[g]+=a;sx[g]+=a*cond[j];sy[g]+=a*y(i,j);sxy[g]+=a*y(i,j)*cond[j];}
  double num=0,den=0;
  for(int g=0;g<G;g++){num+=sxy[g]-sx[g]*sy[g]/sw[g];den+=sx[g]-sx[g]*sx[g]/sw[g];}
  out[i]=num/den;
 }
 return out;
}

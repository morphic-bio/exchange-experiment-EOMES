#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// Each column is a ranking encoded as gene-index -> one-based rank.
// [[Rcpp::export]]
NumericMatrix score_rank_matrix(IntegerMatrix ranks,List members) {
 int N=ranks.nrow(),B=ranks.ncol(),P=members.size();
 NumericMatrix out(P,B);
 for(int p=0;p<P;p++) {
  IntegerVector ids=members[p];int m=ids.size();std::vector<int> pos(m);
  for(int b=0;b<B;b++) {
   for(int k=0;k<m;k++)pos[k]=ranks(ids[k]-1,b);
   std::sort(pos.begin(),pos.end());double hi=0,lo=0;
   for(int k=0;k<m;k++) {
    double miss=double(pos[k]-(k+1))/(N-m);
    hi=std::max(hi,double(k+1)/m-miss);
    lo=std::min(lo,double(k)/m-miss);
   }
   out(p,b)=hi>-lo?hi:(hi<-lo?lo:0);
  }
  if(p%250==0)checkUserInterrupt();
 }
 return out;
}
// [[Rcpp::export]]
List calibrate_wt(NumericMatrix obs,NumericMatrix nulls) {
 int P=obs.nrow(),B=obs.ncol(),K=nulls.ncol();
 NumericMatrix pv(P,B),nes(P,B);IntegerMatrix exceed(P,B),nsign(P,B);
 for(int p=0;p<P;p++) {
  std::vector<double> v(K);int np=0,nn=0;double sp=0,sn=0;
  for(int k=0;k<K;k++){double x=nulls(p,k);v[k]=x;if(x>=0){np++;sp+=x;}if(x<=0){nn++;sn-=x;}}
  std::sort(v.begin(),v.end());
  for(int b=0;b<B;b++) {
   double e=obs(p,b);bool positive=e>=0;int n=positive?np:nn;
   int ne=positive ? v.end()-std::lower_bound(v.begin(),v.end(),e-1e-12) : std::upper_bound(v.begin(),v.end(),e+1e-12)-v.begin();
   nsign(p,b)=n;exceed(p,b)=ne;
   pv(p,b)=e==0?1:(n>0?std::min(1.0,double(1+ne)/(1+n)):NA_REAL);
   nes(p,b)=e==0?0:(n>0?e/(positive?sp/n:sn/n):NA_REAL);
  }
 }
 return List::create(_["pvalue"]=pv,_["NES"]=nes,_["n_exceedances"]=exceed,_["n_null_same_sign"]=nsign);
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = "distance")]]
double cal(NumericMatrix data_resample,NumericMatrix data_base, NumericMatrix data_experiment, int ancestry){
  int k = data_base.nrow();
  double RMSa=0;
  double RMSb=0;
  double a;
  double b;
  for (int i=0; i<k; i++){
    for (int j=0; j<ancestry; j++){
      a=pow(data_base(i,j)-data_experiment(i,j),2)+a;
      b=pow(data_resample(i,j)-data_experiment(i,j),2)+b;
    }
    RMSa = sqrt(a/ancestry)+RMSa;
    RMSb = sqrt(b/ancestry)+RMSb;
  }
  return (RMSb-RMSa)/RMSb;
}

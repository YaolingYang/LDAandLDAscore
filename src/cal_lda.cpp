#include <Rcpp.h>
using namespace Rcpp;

//' @title LDA of a pair of SNPs
//' @description Computation of the pairwise Linkage Disequilibrium of Ancestry (LDA) between a pair of single nucleotide polymorphisms (SNPs).
//' @param data_resample a data frame of the first SNP's ancestry probabilities after resampling.
//' Different ancestry probabilities are in different columns.
//' @param data_base a data frame of the first SNP's ancestry probabilities.
//' Different ancestry probabilities are in different columns.
//' @param data_experiment a data frame of the second SNP's ancestry probabilities.
//' Different ancestry probabilities are in different columns.
//' @param n_ancestry a positive integer representing the number of different ancestries.
//' @return a numeric number representing the pairwise LDA of the two SNPs.
//' @details This function computes the LDA between two SNPs.
//' Resampling of one of the SNPs' painting data is required prior to implementing this function.
//' To compute pairwise LDA between multiple pairs of SNPs, please use \code{\link{LDA}}.
//' @references Barrie W, Yang Y, Attfield K E, et al. Genetic risk for Multiple Sclerosis originated in Pastoralist Steppe populations. bioRxiv (2022).
//' @examples
//' \donttest{
//' # compute the LDA between the 50th SNP and the 55th SNP
//' # painting data for the 50th SNP (2 ancestries)
//' data_base <- cbind(LDAandLDAS::example_painting_p1[,50],
//'                    LDAandLDAS::example_painting_p2[,50])
//'
//' # painting data for the 55th SNP (2 ancestries)
//' data_experiment <- cbind(LDAandLDAS::example_painting_p1[,55],
//'                          LDAandLDAS::example_painting_p2[,55])
//'
//' # resample painting data for the 50th SNP
//' data_resample <- data_base[sample(1:nrow(data_base)),]
//'
//' #compute their pairwise LDA
//' LDA_value <- cal_lda(data_resample,data_base,data_experiment,2)
//' }
//' @export
// [[Rcpp::export(name = "cal_lda")]]
double cal(NumericMatrix data_resample,NumericMatrix data_base, NumericMatrix data_experiment, int n_ancestry){
  int k = data_base.nrow();
  double RMSa=0;
  double RMSb=0;
  double a=0;
  double b=0;
  for (int i=0; i<k; i++){
    for (int j=0; j<n_ancestry; j++){
      a=pow(data_base(i,j)-data_experiment(i,j),2)+a;
      b=pow(data_resample(i,j)-data_experiment(i,j),2)+b;
    }
    RMSa = sqrt(a/n_ancestry)+RMSa;
    RMSb = sqrt(b/n_ancestry)+RMSb;
  }
  return (RMSb-RMSa)/RMSb;
}

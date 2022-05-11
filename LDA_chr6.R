jobid <- as.integer(commandArgs(trailingOnly = TRUE))

library(Rcpp)

load('Steppe_chr6.Rdata')
load('African_chr6.Rdata')
load('CHG_chr6.Rdata')
load('EHG_chr6.Rdata')
load('WHG_chr6.Rdata')
load('Farmer_chr6.Rdata')

sourceCpp(code='#include <Rcpp.h>
 using namespace Rcpp;

 // [[Rcpp::export(name = "cal_d_originalrcpp")]]
            double cal(NumericMatrix data_resample,NumericMatrix data_base, NumericMatrix data_experiment){
              int k = data_base.nrow();
              double RMSa=0;
              double RMSb=0;
              double a;
              double b;
              for (int i=0; i<k; i++){
                a=pow(data_base(i,0)-data_experiment(i,0),2);
                a=pow(data_base(i,1)-data_experiment(i,1),2)+a;
                a=pow(data_base(i,2)-data_experiment(i,2),2)+a;
                a=pow(data_base(i,3)-data_experiment(i,3),2)+a;
                a=pow(data_base(i,4)-data_experiment(i,4),2)+a;
                a=pow(data_base(i,5)-data_experiment(i,5),2)+a;
                RMSa = sqrt(a/6)+RMS;
                b=pow(data_resample(i,0)-data_experiment(i,0),2);
                b=pow(data_resample(i,1)-data_experiment(i,1),2)+a;
                b=pow(data_resample(i,2)-data_experiment(i,2),2)+a;
                b=pow(data_resample(i,3)-data_experiment(i,3),2)+a;
                b=pow(data_resample(i,4)-data_experiment(i,4),2)+a;
                b=pow(data_resample(i,5)-data_experiment(i,5),2)+a;
                RMSb = sqrt(b/6)+RMS;
              }
              return (RMSb-RMSa)/RMSb;
              }')



nsnp<-c(0,575,575,585,585,590,590,700,700,700,700,725,725,750,800,850,850,950,1000,
        1050,1100,1150,1200,1250,1250,1350,1400,1450,1550,1600,1700,1800,1900,2500,3776)

cal_lda_mhc <- function(bc4id){

    lda_mhc<- as.data.frame(matrix(NA,nrow=38976,ncol=nsnp[bc4id+1]))
    k=1
    sum_before_snp <- sum(nsnp[1:(bc4id)])+1
    sum_snp <- sum(nsnp[1:(bc4id+1)])
    for (i in sum_before_snp:sum_snp){
      data_snp1 <- cbind(African[,i+1],CHG[,i+1],
                         EHG[,i+1],Farmer[,i+1],
                         WHG[,i+1],Steppe[,i+1])
      data_resample <- data_snp1[sample(1:nrow(data_snp1)),]
      for(j in i:38976){
        if(j==i){
          lda_mhc[j,k] <- 0
        }else{
          data_snp2 <- cbind(African[,j+1],CHG[,j+1],
                             EHG[,j+1],Farmer[,j+1],
                             WHG[,j+1],Steppe[,j+1])
          lda[j,k] <- cal_d_originalrcpp(data_resample,data_snp1,data_snp2)
        }
      }
      k=k+1
      print(i)
    }

  filename1 <- paste("LDA_chr6",bc4id,".csv",sep='')
  #filename2 <- paste("d_chr6",bc4id,".csv",sep='')
  write.table(lda_mhc,file=filename1,sep=',',row.names=FALSE)
  #write.table(d,file=filename2,sep=',',row.names=FALSE)
}

set.seed(123)
for(q in 1:34)
{
  if(jobid==q){
    cal_lda_mhc(q)
  }
}


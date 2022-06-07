library(Rcpp)

findrate <- function(x){
        rate=chrmhcmap[nrow(chrmhcmap[chrmhcmap[,2]<=x,]),3]
        return(rate)
}

findmorgan <- function(x){
        rownumber <- nrow(chrmhcmap[chrmhcmap[,2]<=x,])
        morgan <- chrmhcmap[rownumber,1]+chrmhcmap[rownumber,3]*(x-chrmhcmap[rownumber,2])/1000000
        return(morgan)
}


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
                RMSa = sqrt(a/2)+RMSa;
                b=pow(data_resample(i,0)-data_experiment(i,0),2);
                b=pow(data_resample(i,1)-data_experiment(i,1),2)+b;
                RMSb = sqrt(b/2)+RMSb;
              }
              return (RMSb-RMSa)/RMSb;
              }')

#load painting data for ancestry 1 and 2
data1 <- read.table('painting_p1.csv',sep=',',header=TRUE)
data2 <- read.table('painting_p2.csv',sep=',',header=TRUE)

n_snp <- ncol(data1)-1

#calculate pairwise lda
cal_ldaold <- function(n_snp){
        lda<- as.data.frame(matrix(NA,nrow=n_snp,ncol=n_snp))
        k=1
        for (i in 1:(n_snp-1)){
                data_snp1 <- cbind(data1[,i+1],data2[,i+1])
                data_resample <- data_snp1[sample(1:nrow(data_snp1)),]
                for(j in (i+1):n_snp){
                        data_snp2 <- cbind(data1[,j+1],data2[,j+1])
                        lda[j,k] <- cal_d_originalrcpp(data_resample,data_snp1,data_snp2)

                }
                k=k+1
                print(i)
        }
        return(lda)
}

cal_lda <- function(n_snp){
        lda<- as.data.frame(matrix(NA,nrow=n_snp,ncol=n_snp))
        k=1
        for (i in 1:(n_snp-1)){
                data_snp1 <- cbind(data1[,i+1],data2[,i+1])
                data_resample <- data_snp1[sample(1:nrow(data_snp1)),]
                if(i<n_snp-1400){
                        q=i+1400
                }else{
                        q=n_snp
                }
                for(j in (i+1):q){
                        data_snp2 <- cbind(data1[,j+1],data2[,j+1])
                        lda[j,k] <- cal_d_originalrcpp(data_resample,data_snp1,data_snp2)
                        
                }
                k=k+1
                print(i)
        }
        return(lda)
}

LDA=cal_lda(n_snp)
LDA[is.na(LDA)] <- 0
LDA=LDA+t(LDA)
diag(LDA)=1
colnames(LDA)=row.names(LDA)=gsub('X','',colnames(data1)[-1])

#mapping data contains pd,gd,rate
colnames(data1) <- gsub('X','',colnames(data1))
colnames(data2) <- gsub('X','',colnames(data2))
map <- data.frame(pd=as.numeric(colnames(data1)[-1]),
                  gd=as.numeric(colnames(data1)[-1])/1000000,
                  rate=rep(1,n_snp))

#calculate LDA score
LDAS_gd <- function(n_snp,window){
        cal_average_lda <- function(n){
                ave <- (LDA_use[n]+LDA_use[(n+1)])/2
                return(ave)
        }
        LDA_score <- vector()
        for (j in 1:n_snp){
                #print(j)
                
                # the number of SNPs within 3Mb window left and right to the SNP
                n_snps1<-length(which(map[1:j,2]<=(map[j,2]+window)))-1
                n_snps2<-length(which(map[j:n_snp,2]>=(map[j,2]-window)))-1
                
                # the genetic distance gap between every two SNPs
                if(j==1){
                        snp_gd_gap <- abs(map[(j+1):(j+n_snps2),2]-
                                                  map[j:(j+n_snps2-1),2])
                        LDA_use <- as.numeric(c(1,LDA[(j+1):(j+n_snps2),j]))
                }else{
                        if(j==n_snp){
                                snp_gd_gap <- abs(map[(j-n_snps1):(j-1),2]-
                                                          map[(j-n_snps1+1):j,2])
                                LDA_use <- as.numeric(c(LDA[j,(j-n_snps1):(j-1)],1))
                        }else{
                                snp_gd_gap <- c(abs(map[(j-n_snps1):(j-1),2]-
                                                            map[(j-n_snps1+1):j,2]),
                                                abs(map[(j+1):(j+n_snps2),2]-
                                                            map[j:(j+n_snps2-1),2]))
                                LDA_use <- as.numeric(c(LDA[j,(j-n_snps1):(j-1)],1,LDA[(j+1):(j+n_snps2),j]))
                        }
                }
                
                #LDA_use <- LDA[(j-n_snps1):(j+n_snps2),j] #use these LDA data
                
                # the average LDA of two SNPs with respect to the jth SNP
                lda_score_j <- sapply(1:(n_snps1+n_snps2),cal_average_lda)
                
                #left end
                if(map[j,2]<window+map[n_snp,2]){
                        
                        LDA_right=as.numeric(LDA[j,j-n_snps1-1])
                        gd_gap=map[j-n_snps1-1,2]-map[j-n_snps1,2]
                        gd_to_end=window-map[j-n_snps1,2]+map[j,2]
                        LDA_right_ave = (LDA_use[1]+gd_to_end/gd_gap*(LDA_right-LDA_use[1]))/2
                        
                        
                        gd_left <- map[n_snp,2]
                        n_snps3<-n_snps1-length(which(map[1:j,2]<(map[j,2]*2-gd_left)))+1
                        gd_gap_add <- map[n_snp-n_snps1-n_snps2+n_snps3-1,2]-(map[j,2]*2-gd_left)
                        lda_add <- (LDA_use[n_snps1+n_snps2+1]+LDA_use[n_snps3])/2
                        if(n_snps3==1){
                                snp_gd_gap <- c(gd_to_end,snp_gd_gap,gd_gap_add,gd_to_end)
                                lda_score_j <- c(LDA_right_ave,lda_score_j,lda_add,LDA_right_ave)
                        }else{
                                if(n_snps3!=0){
                                        snp_gd_gap <- c(gd_to_end,snp_gd_gap,gd_gap_add,
                                                        snp_gd_gap[1:(n_snps3-1)],gd_to_end)
                                        lda_score_j <- c(LDA_right_ave,lda_score_j,lda_add,
                                                         lda_score_j[1:(n_snps3-1)],LDA_right_ave)
                                }else
                                {
                                        snp_gd_gap <- c(gd_to_end,snp_gd_gap)
                                        lda_score_j <- c(LDA_right_ave,lda_score_j)
                                        
                                        gd_to_end=gd_left+window-map[j,2]
                                        gd_gap=map[j-n_snps1-1,2]-map[j,2]-window+gd_to_end
                                        LDA_right_ave = (LDA_use[n_snps1+n_snps2+1]+
                                                                 gd_to_end/gd_gap*(LDA_right-LDA_use[n_snps1+n_snps2+1]))/2
                                        
                                        snp_gd_gap <- c(snp_gd_gap,gd_to_end)
                                        lda_score_j <- c(lda_score_j,LDA_right_ave)
                                }
                        }
                        
                }else{
                        if(map[j,2]>map[1,2]-window){
                                
                                #left in real, but right in our data (right is big pd)
                                LDA_left=as.numeric(LDA[j+n_snps2+1,j])
                                gd_gap=map[j+n_snps2,2]-map[j+n_snps2+1,2]
                                gd_to_end=map[j+n_snps2,2]-map[j,2]+window
                                LDA_left_ave = (LDA_use[n_snps1+n_snps2+1]+
                                                        gd_to_end/gd_gap*(LDA_left-LDA_use[n_snps1+n_snps2+1]))/2
                                
                                gd_right <- map[1,2]
                                #Assume A is the gd of the jth SNP, and B is the closest SNP to right end
                                #Then n_snps3 is the number of SNPs within (A-4,A-(B-A))
                                n_snps3<-n_snps2-length(which(map[j:n_snp,2]>(map[j,2]*2-gd_right)))+1
                                gd_gap_add <- (map[j,2]*2-gd_right)-map[n_snps1+n_snps2-n_snps3+2,2]
                                lda_add <- (LDA_use[1]+LDA_use[n_snps1+n_snps2-n_snps3+2])/2
                                
                                if(n_snps3==1){
                                        snp_gd_gap <- c(gd_to_end,gd_gap_add,snp_gd_gap,gd_to_end)
                                        lda_score_j <- c(LDA_left_ave,lda_add,lda_score_j,LDA_left_ave)
                                }else{
                                        if(n_snps3!=0){
                                                snp_gd_gap <- c(gd_to_end,snp_gd_gap[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                                                                gd_gap_add,snp_gd_gap,gd_to_end)
                                                lda_score_j <- c(LDA_left_ave,
                                                                 lda_score_j[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                                                                 lda_add,lda_score_j,LDA_left_ave)
                                        }else
                                        {
                                                snp_gd_gap <- c(snp_gd_gap,gd_to_end)
                                                lda_score_j <- c(lda_score_j,LDA_left_ave)
                                                
                                                gd_to_end=window+map[j,2]-gd_right
                                                gd_gap=map[j,2]-window-map[j+n_snps1+1,2]+gd_to_end
                                                LDA_left_ave = (LDA_use[1]+gd_to_end/gd_gap*(LDA_left-LDA_use[1]))/2
                                                
                                                snp_gd_gap <- c(gd_to_end,snp_gd_gap)
                                                lda_score_j <- c(LDA_left_ave,lda_score_j)
                                        }
                                }
                                
                        }else{
                                LDA_right=as.numeric(LDA[j,j-n_snps1-1])
                                gd_gap_right=map[j-n_snps1-1,2]-map[j-n_snps1,2]
                                gd_to_end_right=window-map[j-n_snps1,2]+map[j,2]
                                LDA_right_ave = (LDA_use[1]+gd_to_end_right/gd_gap_right*(LDA_right-LDA_use[1]))/2
                                
                                
                                LDA_left=as.numeric(LDA[j+n_snps2+1,j])
                                gd_gap_left=map[j+n_snps2,2]-map[j+n_snps2+1,2]
                                gd_to_end_left=map[j+n_snps2,2]-map[j,2]+window
                                LDA_left_ave = (LDA_use[n_snps1+n_snps2+1]+
                                                        gd_to_end_left/gd_gap_left*(LDA_left-LDA_use[n_snps1+n_snps2+1]))/2
                                
                                snp_gd_gap <- c(gd_to_end_right,snp_gd_gap,gd_to_end_left)
                                lda_score_j <- c(LDA_right_ave,lda_score_j,LDA_left_ave)
                        }
                }
                
                #right end
                
                if(length(lda_score_j)!=length(snp_gd_gap)){print(j)}
                # compute the LDA score
                LDA_score[j]=sum(lda_score_j*snp_gd_gap)
        }
        
        LDA_score <- cbind(map[,1],LDA_score)
        colnames(LDA_score)[2]='LDAS'
        return(LDA_score)
}

LDAscore <- LDAS_gd(n_snp,window=5)

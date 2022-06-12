#' Pairwise Linkage Disequilibrium of Ancestry (LDA)
#'
#' @param data data frame of the painting for all the ancestries (combined in a single file).
#' @param ancestry number of ancestries.
#' @param printProcess logical. Print the process of calculating the pairwise LDA for the i-th SNP.
#' @param SNPlimit maximum number of SNPs at each side of the SNP used to calculate the pairwise LDA. The value shouldn't be larger than the total number of SNPs.
#'
#' @return a dataframe of the correlation matrix.
#' @export
#'
#' @examples
#' library(LDAandLDAS)
#' #combine the painting data for two ancestries
#' data=cbind(painting_p1[,-1],painting_p2[,-1])
#' #calculate the pairwise LDA of SNPs
#' LDA_result <- LDA(data,ancestry=2)
#'
LDA <- function(data,ancestry,printProcess=TRUE,SNPlimit=NA){
  Rcpp::sourceCpp(file='src/distance.cpp')

  n_snp <- ncol(data)/ancestry

  if(is.na(SNPlimit)){
    SNPlimit=n_snp
  }

  if(SNPlimit>n_snp){
    SNPlimit=n_snp
    print('SNPlimit must not be bigger than the number of SNPs!')
  }

  #calculate pairwise lda
    lda<- as.data.frame(matrix(NA,nrow=n_snp,ncol=n_snp))
    k=1
    for (i in 1:(n_snp-1)){
      data_snp1 <- cbind(data[,i],data[,i+n_snp])
      if(ancestry>2){
        for(m in 2:ancestry){
          data_snp1 <- cbind(data_snp1,data[,i+(m-1)*n_snp])
        }
      }

      data_resample <- data_snp1[sample(1:nrow(data_snp1)),]
      if(i<n_snp-SNPlimit){
        SNPnumberlimit=i+SNPlimit
      }else{
        SNPnumberlimit=n_snp
      }
      for(j in (i+1):SNPnumberlimit){
        data_snp2 <- cbind(data[,j],data[,j+n_snp])
        if(ancestry>2){
          for(m in 2:ancestry){
            data_snp2 <- cbind(data_snp2,data[,j+(m-1)*n_snp])
          }
        }
        lda[j,k] <- distance(data_resample,data_snp1,data_snp2,ancestry)

      }
      k=k+1
      if(printProcess) cat("Calculating pairwise LDA of SNP",i,'\n');
    }

  lda[is.na(lda)] <- 0
  lda=lda+t(lda)
  diag(lda)=1
  colnames(lda)=row.names(lda)=gsub('X','',colnames(data)[1:n_snp])
  return(lda)
}


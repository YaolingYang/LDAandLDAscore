#' Pairwise Linkage Disequilibrium of Ancestry (LDA)
#'
#' @param data a list of data frames of the N*k painting data for ancestries (N is the number of genomes, k is the number of SNPs). SNPs should be in the decreasing order.
#' @param printProcess logical. Print the process of calculating the pairwise LDA for the i-th SNP.
#' @param SNPlimit maximum number of SNPs at each side of the SNP used to calculate the pairwise LDA. The value shouldn't be larger than the total number of SNPs.
#'
#' @return a dataframe of the pairwise LDA, with SNPs in the decreasing order.
#' @export
#'
#' @examples
#' library(LDAandLDAS)
#' #combine the painting data for two ancestries
#' data=list(painting_p1[,-1],painting_p2[,-1])
#' #calculate the pairwise LDA of SNPs
#' LDA_result <- LDA(data,SNPlimit=1000)
#'
LDA <- function(data,printProcess=TRUE,SNPlimit=NA){

  n_snp <- ncol(data[[1]])

  ancestry <- length(data)

  #combine the painting data
  data_combine <- Reduce(cbind,data)

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
      data_snp1 <- cbind(data_combine[,i],data_combine[,i+n_snp])
      if(ancestry>2){
        for(m in 2:ancestry){
          data_snp1 <- cbind(data_snp1,data_combine[,i+(m-1)*n_snp])
        }
      }

      data_resample <- data_snp1[sample(1:nrow(data_snp1)),]
      if(i<n_snp-SNPlimit){
        SNPnumberlimit=i+SNPlimit
      }else{
        SNPnumberlimit=n_snp
      }
      for(j in (i+1):SNPnumberlimit){
        data_snp2 <- cbind(data_combine[,j],data_combine[,j+n_snp])
        if(ancestry>2){
          for(m in 2:ancestry){
            data_snp2 <- cbind(data_snp2,data_combine[,j+(m-1)*n_snp])
          }
        }
        lda[j,k] <- cal_lda(data_resample,data_snp1,data_snp2,ancestry)

      }
      k=k+1
      if(printProcess) cat("Calculating pairwise LDA of SNP",i,'\n');
    }

  lda[is.na(lda)] <- 0
  lda=lda+t(lda)
  diag(lda)=1
  colnames(lda)=row.names(lda)=colnames(data_combine)[1:n_snp]
  return(lda)
}


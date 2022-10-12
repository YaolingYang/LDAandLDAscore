#' @title LDA of all pairs of SNPs
#' @description Computation of the pairwise Linkage Disequilibrium of Ancestry (LDA) between all pairs of single nucleotide polymorphisms (SNPs).
#' @param paintings a list of data frames of the N*k painting data
#' (average ancestry probabilities) from different populations
#' (N is the number of genomes, k is the number of SNPs).
#' @param SNPlimit a positive integer representing the maximum number of SNPs at each side of
#' a SNP that is used to calculate the pairwise LDA for the SNP.
#' The value shouldn't be larger than the total number of SNPs.
#' We may set a limit if the LDAs between SNPs far in distance are not to be investigated.
#' @param verbose logical. If verbose=TRUE, print the process of calculating the pairwise LDA for the i-th SNP.
#' By default, verbose=FALSE
#'
#' @return a data frame of the pairwise LDA, with SNPs in the decreasing
#' order of physical position on a chromosome.
#'
#' @details Linkage Disequilibrium of Ancestry (LDA) quantifies the correlations
#' between the ancestry of two SNPs, measuring the proportion of individuals who
#' have experienced a recombination leading to a change in ancestry,
#' relative to the genome-wide baseline.
#'
#' @references Barrie W, Yang Y, Attfield K E, et al. Genetic risk for Multiple Sclerosis originated in Pastoralist Steppe populations. bioRxiv (2022).
#'
#' @examples
#' \donttest{
#' # visualize the painting data
#' # Painting data are the average probabilities of different populations
#' head(LDAandLDAS::example_painting_p1[1:5,],10)
#'
#' # combine the painting data for two populations as a list
#' # to make to input data for function 'LDA'.
#' paintings=list(LDAandLDAS::example_painting_p1,
#'           LDAandLDAS::example_painting_p2)
#'
#' # calculate the pairwise LDA of SNPs
#' LDA_result <- LDA(paintings)
#' }
#'
#' @export
#'
LDA <- function(paintings,SNPlimit=NULL,verbose=FALSE){

  n_snp <- ncol(paintings[[1]])

  n_ancestry <- length(paintings)

  #combine the painting data
  data_combine <- Reduce(cbind,paintings)

  if(is.null(SNPlimit)){
    SNPlimit=n_snp
  }

  if(SNPlimit>n_snp){
    SNPlimit=n_snp
    warning('SNPlimit must not be bigger than the number of SNPs!')
  }

  #calculate pairwise lda
    lda<- as.data.frame(matrix(NA,nrow=n_snp,ncol=n_snp))
    k=1
    for (i in 1:(n_snp-1)){
      data_hap1 <- cbind(data_combine[,i],data_combine[,i+n_snp])
      if(n_ancestry>2){
        for(m in 2:n_ancestry){
          data_hap1 <- cbind(data_hap1,data_combine[,i+(m-1)*n_snp])
        }
      }

      data_resample <- data_hap1[sample(1:nrow(data_hap1)),]
      if(i<n_snp-SNPlimit){
        SNPnumberlimit=i+SNPlimit
      }else{
        SNPnumberlimit=n_snp
      }
      for(j in (i+1):SNPnumberlimit){
        data_hap2 <- cbind(data_combine[,j],data_combine[,j+n_snp])
        if(n_ancestry>2){
          for(m in 2:n_ancestry){
            data_hap2 <- cbind(data_hap2,data_combine[,j+(m-1)*n_snp])
          }
        }
        lda[j,k] <- cal_lda(data_resample,data_hap1,data_hap2,n_ancestry)

      }
      k=k+1
      if(verbose) cat("Calculating pairwise LDA of SNP",i,'\n');
    }

  lda[is.na(lda)] <- 0
  lda=lda+t(lda)
  diag(lda)=1
  colnames(lda)=row.names(lda)=colnames(data_combine)[1:n_snp]
  return(lda)
}


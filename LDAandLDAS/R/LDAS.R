#'Linkage Disequilibrium of Ancestry Score (LDAS)
#'
#' @param LDA data frame of the pairwise LDA of all the SNPs.
#' @param map data frame of the physical position and genetic distance of all the SNPs.
#' @param window number of the window of genetic distance. The default window is 5.
#' @param printProcess logical. Print the process of calculating the LDA score for the i-th SNP.
#'
#' @return a data frame of the LDA score and its upper and lower bound at the physical position of each SNP.
#' @export
#'
#' @examples
#' library(LDAandLDAS)
#' #combine the painting data for two ancestries
#' data=cbind(painting_p1[,-1],painting_p2[,-1])
#' #calculate the pairwise LDA of SNPs
#' LDA_result <- LDA(data,ancestry=2)
#' #map is the data containing the physical position and recombination distance of the SNPs
#' #calculate the LDA score for the SNPs
#' LDA_score <- LDAS(LDA_result,map,window=4)
#'

LDAS <- function(LDA,map,window,printProcess=TRUE){
  n_snp <- nrow(map)
  cal_average_lda <- function(n){
    ave <- (LDA_use[n]+LDA_use[n+1])/2
    return(ave)
  }
  max_lda <- function(n){
    max_lda <- max(LDA_use[n],LDA_use[n+1])
    return(max_lda)
  }
  min_lda <- function(n){
    min_lda <- min(LDA_use[n],LDA_use[n+1])
    return(min_lda)
  }

  LDA_score <- vector()
  LDA_score_max <- vector()
  LDA_score_min <- vector()

  for (j in 1:n_snp){
    if(printProcess) cat("Calculating LDA score of SNP",j,'\n');

    # the number of SNPs within 5cM window left and right to the SNP
    # Note: n_snps1 is left in our data but right in reality
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
    lda_score_max_j <- sapply(1:(n_snps1+n_snps2),max_lda)
    lda_score_min_j <- sapply(1:(n_snps1+n_snps2),min_lda)

    #left end in reality but right end in the data
    if(map[j,2]<window+map[n_snp,2]){

      LDA_right=as.numeric(LDA[j,j-n_snps1-1])
      gd_gap=map[j-n_snps1-1,2]-map[j-n_snps1,2]
      gd_to_end=window-map[j-n_snps1,2]+map[j,2]
      LDA_right_ave = (LDA_use[1]+gd_to_end/gd_gap*(LDA_right-LDA_use[1]))/2

      LDA_right_ave_max = max(LDA_right,LDA_use[1])
      LDA_right_ave_min = min(LDA_right,LDA_use[1])


      gd_left <- map[n_snp,2]
      n_snps3<-n_snps1-length(which(map[1:j,2]<(map[j,2]*2-gd_left)))+1
      gd_gap_add <- map[n_snp-n_snps1-n_snps2+n_snps3-1,2]-(map[j,2]*2-gd_left)
      lda_add <- (LDA_use[n_snps1+n_snps2+1]+LDA_use[n_snps3])/2
      lda_add_max <- max(LDA_use[n_snps1+n_snps2+1],LDA_use[n_snps3])
      lda_add_min <- 0

      if(n_snps3==1){
        snp_gd_gap <- c(gd_to_end,snp_gd_gap,gd_gap_add,gd_to_end)
        lda_score_j <- c(LDA_right_ave,lda_score_j,lda_add,LDA_right_ave)
        lda_score_max_j <- c(LDA_right_ave_max,lda_score_max_j,lda_add_max,LDA_right_ave_max)
        lda_score_min_j <- c(LDA_right_ave_min,lda_score_min_j,lda_add_min,LDA_right_ave_min)

      }else{
        if(n_snps3!=0){
          snp_gd_gap <- c(gd_to_end,snp_gd_gap,gd_gap_add,
                          snp_gd_gap[1:(n_snps3-1)],gd_to_end)
          lda_score_j <- c(LDA_right_ave,lda_score_j,lda_add,
                           lda_score_j[1:(n_snps3-1)],LDA_right_ave)
          lda_score_max_j <- c(LDA_right_ave_max,lda_score_max_j,lda_add_max,
                               lda_score_max_j[1:(n_snps3-1)],LDA_right_ave_max)
          lda_score_min_j <- c(LDA_right_ave_min,lda_score_min_j,lda_add_min,
                               lda_score_min_j[1:(n_snps3-1)],LDA_right_ave_min)
        }else
        {
          snp_gd_gap <- c(gd_to_end,snp_gd_gap)
          lda_score_j <- c(LDA_right_ave,lda_score_j)
          lda_score_max_j <- c(LDA_right_ave_max,lda_score_max_j)
          lda_score_min_j <- c(LDA_right_ave_min,lda_score_min_j)

          gd_to_end=gd_left+window-map[j,2]
          gd_gap=map[j-n_snps1-1,2]-map[j,2]-window+gd_to_end
          LDA_right_ave = (LDA_use[n_snps1+n_snps2+1]+
                             gd_to_end/gd_gap*(LDA_right-LDA_use[n_snps1+n_snps2+1]))/2

          LDA_right_ave_max = max(LDA_right,LDA_use[n_snps1+n_snps2+1])

          LDA_right_ave_min = min(LDA_right,LDA_use[n_snps1+n_snps2+1])

          snp_gd_gap <- c(snp_gd_gap,gd_to_end)
          lda_score_j <- c(lda_score_j,LDA_right_ave)
          lda_score_max_j <- c(lda_score_max_j,LDA_right_ave_max)
          lda_score_min_j <- c(lda_score_min_j,LDA_right_ave_min)
        }
      }

    }else{
      #right end
      if(map[j,2]>map[1,2]-window){

        #left in real, but right in our data (right is small pd)
        LDA_left=as.numeric(LDA[j+n_snps2+1,j])
        gd_gap=map[j+n_snps2,2]-map[j+n_snps2+1,2]
        gd_to_end=map[j+n_snps2,2]-map[j,2]+window
        LDA_left_ave = (LDA_use[n_snps1+n_snps2+1]+
                          gd_to_end/gd_gap*(LDA_left-LDA_use[n_snps1+n_snps2+1]))/2
        LDA_left_ave_max = max(LDA_left,LDA_use[n_snps1+n_snps2+1])
        LDA_left_ave_min = min(LDA_left,LDA_use[n_snps1+n_snps2+1])

        gd_right <- map[1,2]
        #Assume A is the gd of the jth SNP, and B is the gd of the closest SNP to right end
        #Then n_snps3 is the number of SNPs within (A-5,A-(B-A))
        n_snps3<-n_snps2-length(which(map[j:n_snp,2]>(map[j,2]*2-gd_right)))+1
        gd_gap_add <- (map[j,2]*2-gd_right)-map[n_snps1+n_snps2-n_snps3+2,2]
        lda_add <- (LDA_use[1]+LDA_use[n_snps1+n_snps2-n_snps3+2])/2
        lda_add_max <- max(LDA_use[1],LDA_use[n_snps1+n_snps2-n_snps3+2])
        lda_add_min <- 0

        if(n_snps3==1){
          snp_gd_gap <- c(gd_to_end,gd_gap_add,snp_gd_gap,gd_to_end)
          lda_score_j <- c(LDA_left_ave,lda_add,lda_score_j,LDA_left_ave)
          lda_score_max_j <- c(LDA_left_ave_max,lda_add_max,lda_score_max_j,LDA_left_ave_max)
          lda_score_min_j <- c(LDA_left_ave_min,lda_add_min,lda_score_min_j,LDA_left_ave_min)
        }else{
          if(n_snps3!=0){
            snp_gd_gap <- c(gd_to_end,snp_gd_gap[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                            gd_gap_add,snp_gd_gap,gd_to_end)
            lda_score_j <- c(LDA_left_ave,
                             lda_score_j[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                             lda_add,lda_score_j,LDA_left_ave)
            lda_score_max_j <- c(LDA_left_ave_max,
                                 lda_score_max_j[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                                 lda_add_max,lda_score_max_j,LDA_left_ave_max)
            lda_score_min_j <- c(LDA_left_ave_min,
                                 lda_score_min_j[(n_snps1+n_snps2-n_snps3+2):(n_snps1+n_snps2)],
                                 lda_add_min,lda_score_min_j,LDA_left_ave_min)
          }else
          {
            snp_gd_gap <- c(snp_gd_gap,gd_to_end)
            lda_score_j <- c(lda_score_j,LDA_left_ave)
            lda_score_max_j <- c(lda_score_max_j,LDA_left_ave_max)
            lda_score_min_j <- c(lda_score_min_j,LDA_left_ave_min)

            gd_to_end=window+map[j,2]-gd_right
            gd_gap=map[j,2]-window-map[j+n_snps1+1,2]+gd_to_end
            LDA_left_ave = (LDA_use[1]+gd_to_end/gd_gap*(LDA_left-LDA_use[1]))/2
            LDA_left_ave_max = max(LDA_use[1],LDA_left)
            LDA_left_ave_min = min(LDA_use[1],LDA_left)

            snp_gd_gap <- c(gd_to_end,snp_gd_gap)
            lda_score_j <- c(LDA_left_ave,lda_score_j)
            lda_score_max_j <- c(LDA_left_ave_max,lda_score_max_j)
            lda_score_min_j <- c(LDA_left_ave_min,lda_score_min_j)
          }
        }

      }else{
        LDA_right=as.numeric(LDA[j,j-n_snps1-1])
        gd_gap_right=map[j-n_snps1-1,2]-map[j-n_snps1,2]
        gd_to_end_right=window-map[j-n_snps1,2]+map[j,2]
        LDA_right_ave = (LDA_use[1]+gd_to_end_right/gd_gap_right*(LDA_right-LDA_use[1]))/2
        LDA_right_ave_max = max(LDA_use[1],LDA_right)
        LDA_right_ave_min = min(LDA_use[1],LDA_right)


        LDA_left=as.numeric(LDA[j+n_snps2+1,j])

        gd_gap_left=map[j+n_snps2,2]-map[j+n_snps2+1,2]
        gd_to_end_left=map[j+n_snps2,2]-map[j,2]+window
        LDA_left_ave = (LDA_use[n_snps1+n_snps2+1]+
                          gd_to_end_left/gd_gap_left*(LDA_left-LDA_use[n_snps1+n_snps2+1]))/2
        LDA_left_ave_max = max(LDA_use[n_snps1+n_snps2+1],LDA_left)
        LDA_left_ave_min = min(LDA_use[n_snps1+n_snps2+1],LDA_left)

        snp_gd_gap <- c(gd_to_end_right,snp_gd_gap,gd_to_end_left)
        lda_score_j <- c(LDA_right_ave,lda_score_j,LDA_left_ave)
        lda_score_max_j <- c(LDA_right_ave_max,lda_score_max_j,LDA_left_ave_max)
        lda_score_min_j <- c(LDA_right_ave_min,lda_score_min_j,LDA_left_ave_min)
      }
    }

    if(length(lda_score_j)!=length(snp_gd_gap)){print(j)}
    # compute the LDA score
    LDA_score[j]=sum(lda_score_j*snp_gd_gap)
    LDA_score_max[j]=sum(lda_score_max_j*snp_gd_gap)
    LDA_score_min[j]=sum(lda_score_min_j*snp_gd_gap)
  }

  LDA_score <- cbind(map,LDA_score,LDA_score_max,LDA_score_min)
  colnames(LDA_score)[3]='LDAS'
  colnames(LDA_score)[4]='LDAS_max'
  colnames(LDA_score)[5]='LDAS_min'
  return(LDA_score)
}

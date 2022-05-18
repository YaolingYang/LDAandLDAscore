library(reshape2)
library(openxlsx)

getmap <- function(snp){
  ancestrymap <- data.frame(snp)
  for (i in 1:nrow(ancestrymap)){
    ancestrymap[i,2]<-findmorgan(ancestrymap[i,1])
    ancestrymap[i,3]<-findrate(ancestrymap[i,1])
    print(i)
  }
  
  colnames(ancestrymap)[1] <- 'pd'
  colnames(ancestrymap)[2] <- 'gd'
  colnames(ancestrymap)[3] <- 'rate'
  
  
  breakpoint=c(27067780,30467778,32067778,33467778,36567778,40467740,46593783)
  region <- vector()
  
  for (i in 1:nrow(ancestrymap)){
    if(ancestrymap[i,1]<breakpoint[1]){
      region[i]<-'others'
    }else{
      if(ancestrymap[i,1]<breakpoint[2]){
        region[i]<-'p22.1'
      }else{
        if(ancestrymap[i,1]<breakpoint[3]){
          region[i]<-'p21.33'
        }else{
          if(ancestrymap[i,1]<breakpoint[4]){
            region[i]<-'p21.32'
          }else{
            if(ancestrymap[i,1]<breakpoint[5]){
              region[i]<-'p21.31'
            }else{
              if(ancestrymap[i,1]<breakpoint[6]){
                region[i]<-'p21.2'
              }else{
                if(ancestrymap[i,1]<breakpoint[7]){
                  region[i]<-'p21.1'
                }else{
                  region[i]<-'others'
                }
              }
            }
          }
        }
      }
    }
    print(i)
  }
  ancestrymap <- cbind(ancestrymap,region)
  
  return(ancestrymap)
}

#Function for calculating recombination rate given physical position
findrate <- function(x){
  rownumber <- nrow(map[map$pd>=x,])
  if(rownumber==0){
    rate=0
  }else{
    rate=map[rownumber,3]
  }
  return(rate)
}

#Function for calculating genetic distance given physical position
findmorgan <- function(x){
  rownumber <- nrow(map[map$pd>=x,])
  if(rownumber==0){
    gd=0
  }else{
    gd <- map$gd[rownumber]+map$rate[rownumber]*(x-map$pd[rownumber])/1000000
  }
  return(gd)
}


LDAS_gd <- function(n_snp,window){
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
    print(j)
    
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
  colnames(LDA_score)[4]='LDAS'
  colnames(LDA_score)[5]='LDAS_max'
  colnames(LDA_score)[6]='LDAS_min'
  return(LDA_score)
}

setwd("C:/Ubuntu/LDAall/chr6")

load('LDA_chr6.Rdata')

map <- read.table('chr6_paintingmap.csv',sep=',',header=TRUE)

#calculate LDA score
LDA_score <- LDAS_gd(n_snp=nrow(LDA),window=5)

write.table(LDA_score,'LDA_score_maxmin_chr6.csv',sep=',',row.names = FALSE)
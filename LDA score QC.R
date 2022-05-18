#LDA score plot after quality control
LDAS_plot_QC <- function(i){
  #setwd(paste0("C:/Ubuntu/LDAall/chr",i))
  LDA_score<-read.table(paste0('LDA_score_maxmin_chr',i,'.csv'),sep=',',header=TRUE)
  
  LDA_score <- cbind(data.frame(pd=LDA_score$pd,LDAS=LDA_score$LDAS,error=LDA_score[,5]-LDA_score[,6]),
                     LDA_score[,7:10])
  
  remove_index <- which(LDA_score[,3]>=0.5|
                          LDA_score[,4]<=10|
                          LDA_score[,5]<=10|
                          LDA_score[,6]<=10|
                          LDA_score[,7]<=10)
  
  if(length(remove_index)!=0){
    LDA_score <- LDA_score[-remove_index,]
  }
  
  setwd("C:/Ubuntu/LDAall")
  #png(paste0('LDAS_chr',i,'_QC.png'),width=1600,height=700)
  ggplot(LDA_score,aes(x=pd/1000000,y=LDAS))+geom_point(size=1)+
    xlab(paste('Chromosome',i,'position hg19 GRCh37 (Mb)')) + ylab ('LDA score')+
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=15))
  #dev.off()
}

LDAS_plot_QC(1)
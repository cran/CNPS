
siegel_tukey=function(x,y,adjust.median=FALSE,...)
{
  if(!is.numeric(c(x,y))) stop("The input data is not a numeric")

  data=data.frame(c(x,y),rep(c(0,1),c(length(x),length(y))))
  names(data)=c("x","y")
  data=data[order(data$x),]

  if(adjust.median==TRUE){
    data$x[data$y==0]=data$x[data$y==0]-(median(data$x[data$y==0]))
    data$x[data$y==1]=data$x[data$y==1]-(median(data$x[data$y==1]))
  }

  x<-c(x,y)
  sort.x<-sort(data$x)
  sort.id<-data$y[order(data$x)]

  data.matrix<-data.frame(sort.x,sort.id)

  base1<-c(1,4)
  iterator1<-matrix(seq(from=1,to=length(x),by=4))-1
  rank1<-apply(iterator1,1,function(x) x+base1)

  iterator2<-matrix(seq(from=2,to=length(x),by=4))
  base2<-c(0,1)
  rank2<-apply(iterator2,1,function(x) x+base2)


  if(length(rank1)==length(rank2)){
    rank<-c(rank1[1:floor(length(x)/2)],rev(rank2[1:ceiling(length(x)/2)]))
  } else{
    rank<-c(rank1[1:ceiling(length(x)/2)],rev(rank2[1:floor(length(x)/2)]))
  }


  unique.ranks<-tapply(rank,sort.x,mean)
  unique.x<-as.numeric(as.character(names(unique.ranks)))

  rank.matrix<-data.frame(unique.x,unique.ranks)

  ST.matrix<-merge(data.matrix,rank.matrix,by.x="sort.x",by.y="unique.x")



  ranks0<-ST.matrix$unique.ranks[ST.matrix$sort.id==0]
  ranks1<-ST.matrix$unique.ranks[ST.matrix$sort.id==1]

  output<-twosample_test(ranks0,ranks1,score = "wilcoxon",conf.diff = FALSE,...)
  output$method<-"Siegel-Tukey test"
  output$score<-NULL
  output$null.value[1]<-"The variance of the first"
  output$null.value[2]<-"the variance of the second"
  output
}

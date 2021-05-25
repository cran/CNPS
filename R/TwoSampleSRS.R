
FindK<-function(x,M){
  K=M
  tmp=1
  while(tmp<=x){
    tmp=tmp*(K+1)/(K-M+1)
    K=K+1
  }
  return(K-1)
}


GetBinary<-function(x,M,N){
  B<-rep(FALSE,N)
  while(x>=0&M>0){
    K<-FindK(x,M)
    x=x-choose(K,M)
    M=M-1
    B[N-K]=TRUE
  }
  return(B)
}


ReservoirSample<-function(N,samplenum){
  Reservoir<-seq(from=1,to=samplenum)
  i=samplenum+1
  while(i<=N){
    rand<-as.integer(runif(1)*i)+1
    if(rand<=samplenum){
      Reservoir[rand]=i
    }
    i=i+1
  }
  return(Reservoir)
}

TwosampleSRS<-function(data,index,samplenum,method="R"){
  if(method=="R"){
    return(replicate(samplenum,sample(data)))
  }
  
  N<-length(data)
  M<-length(index[index==1])
  
  if(samplenum>choose(N,M))
    stop("sample number is greater than the all possible combinations")
  
  if(method=="S"){
    SampleIndex<-sample(0:(choose(N,M)-1),samplenum)
  }
  else if(method=="W"){
    SampleIndex<-ReservoirSample((choose(N,M)-1),samplenum)
  }
  Index<-matrix(rep(FALSE,samplenum*N),N,samplenum)
  replace<-matrix(rep(-1,samplenum*M),N,samplenum)
  for(i in 1:samplenum){
    Index[,i]<-GetBinary(SampleIndex[i],M,N)
    replace[,i]<-c(data[Index[,i]],data[!Index[,i]])
  }
  return(replace)
}

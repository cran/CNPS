
MultiDimen_test = function(data , stat = "HT",pair=FALSE, method_p = "sampling" ,rank = FALSE, diff = FALSE , samplenum = 1000){

  ################ warning message #########################

  if(!is.matrix(data) & !is.data.frame(data)) stop("The input data must be a matrix or dataframe")
  if(!is.numeric(samplenum)) stop("samplenum input is not true")

  if(pair==TRUE){
    if(stat %in% c("wsum"))  stop("The stat input is not suitable for paired calculation")
    if( length(unique(data[,ncol(data)-1])) != 2) stop("The last two column of the data must only contain two unique numbers if pair = 'TRUE'")
    if(method_p == "asymptotic" & stat != "HT") stop("Asymptotic method can only be used when statistic is HT when pair = 'TRUE'")
    }
  if(pair==FALSE) {
    if(stat %in% c("zmax","zmaxabs")) stop("The stat input is not suitable for not paired calculation")
    if(sum(unique(data[,ncol(data)]) == c(0,1)) != 2 ) stop("The last column of the data must only contain 0 and 1 if pair = 'FALSE'")
    if(method_p == "asymptotic" & stat %in% c("tmax","tmaxabs")) stop("Asymptotic method can not be used when statistic is tmax or tmaxabs and pair = 'FALSE'")
    if(method_p == "asymptotic" & rank == TRUE & stat == "HT") stop("Asymptotic method can not be used when statistic is ranked HT and pair = 'FALSE'")
     }


  parameter_value = c(stat , method_p , rank , diff ,  pair)
  parameter_name = c("stat" , "method_p" , "rank" , "diff","pair")
  parameter_content = c("HT tmax tmaxabs wsum zmax zmaxabs" , "sampling exact asymptotic" , "TRUE FALSE 1 0" , "TRUE FALSE 1 0", "TRUE FALSE 1 0")
  parameter = data.frame(parameter_value , parameter_name , parameter_content)

  for(i in 1:nrow(parameter)){
    if(!parameter[i,1] %in% strsplit(parameter[i,3]," ")[[1]]) stop(paste(parameter[i,2] , "input is not true"))
  }

  if(stat == "wsum" & diff == TRUE & pair == FALSE) warning("The diff is not calculating when stat = 'wsum'")
  if(method_p != "sampling" & !missing(samplenum)) warning("\n samplenum is not working")

  if(stat == "wsum" & pair == FALSE) {
    rank = TRUE
    if(rank == FALSE) warning("Changing the rank to TRUE, because stat = 'wsum' is chosen")
  }

  ##############################################
  if(pair == FALSE){
    HT <- function(dat,index,rank){
      x = dat[index==0,]
      y = dat[index==1,]
      m = dim(x)[1]
      k = dim(x)[2]
      n = dim(y)[1]
      mx = apply(x,2,mean)
      vx = var(x)
      my = apply(y,2,mean)
      vy = var(y)
      v = (m-1)*vx+(n-1)*vy
      v = v/(m+n-2)

      diff = mx-my
      t2 = t(diff)%*%solve(v)%*%diff
      t2 = t2*m*n/(m+n)
      fs = (m+n-k-1)*t2/((m+n-2)*k)
      as.numeric(fs)
    }
    wilcox.sd <- function(x, index){
      rd = rank(x)
      m = sum(index==0)
      nt = length(index)
      n = nt-m
      mw = m*(nt+1)/2
      vw1 = mw*n/6

      ux = unique(x)
      freq = ux
      for(i in 1:length(ux)){
        freq[i] = sum(x==ux[i])
      }
      af = sum(freq^3-freq)*m*n
      tmp = 12*nt*(nt-1)
      af = af/tmp

      vw = vw1-af

      w = sum(rd[index==0])
      (w-mw)/sqrt(vw)
    }

    tall <- function(dat,index=NULL , rank=NULL){
      if(rank == FALSE){
        x = dat[index==0,]
        y = dat[index==1,]
        m = dim(x)[1]
        k = dim(x)[2]
        n = dim(y)[1]
        stata = rep(0,k)
        for(i in 1:k){
          stata[i] = t.test(x[,i],y[,i],var.equal=TRUE)$statistic
        }
        stata
      }
      else if(rank == TRUE){
        k = dim(dat)[2]
        stat = rep(0,k)
        for(i in 1:k){
          stat[i] = wilcox.sd(dat[,i],index)
        }
        stat
      }

    }
    tmax <- function(dat,index=NULL,rank=NULL){
      max(tall(dat,index,rank))
    }
    tmaxabs <- function(dat,index=NULL,rank=NULL){
      max(abs(tall(dat,index,rank)))
    }

    wsum <- function(rd,index , rank){
      sum(rd[index==0,])
    }



  nc = ncol(data)
  index = data[,nc]
  data = data[,-nc]

  STAT = stat
  if(stat == "HT") stat = HT
  else if(stat == "tmax") stat = tmax
  else if(stat == "tmaxabs") stat = tmaxabs
  else if(stat == "wsum") stat = wsum

  if(rank == TRUE){
    for(i in 1:ncol(data)){
      data[,i] = rank(data[,i])
    }
  }
  else if(rank == FALSE) {
  }


  x = data[index==0,]
  y = data[index==1,]
  m = dim(x)[1]
  k = dim(x)[2]
  n = dim(y)[1]
  nt=n+m

  data = as.matrix(data)
  Fobs = stat(data, index , rank)


  if( method_p == "exact"){
    all.comb = combn(nt, m)
    N = dim(all.comb)[2]
    replace = rep(1, N)
    for(i in 1:N){
      ind = all.comb[,i]
      index = rep(0,nt)
      index[ind] = 1
      replace[i] = stat(data,index,rank)
    }
    p.value <- length(replace[replace >= Fobs]) / N
  }

  else if(method_p == "sampling"){
    replace <- replicate(samplenum, stat(data, sample(index) , rank))
    p.value <- length(replace[replace >= Fobs]) / samplenum
  }

  else if( method_p == "asymptotic"){
    if(STAT == "HT" & rank == FALSE){
      p.value = 1-pf(Fobs, k, nt-k-1)
    }

    else if(STAT == "wsum"){
      S =  apply(data,1,sum)
      sigmaz = var(S)*(nt-1)/nt
      varw = m*n*sigmaz/(nt-1)
      mw = k*m*(nt+1)/2
      z = (Fobs - mw) / sqrt(varw)
      p.value = 1 - pnorm(z)
    }
  }


  if(diff == TRUE & STAT != "wsum"){
    replace <- replicate(samplenum, tmaxabs(data, sample(index) , rank))
    quan = quantile(replace , 0.95)
    large = which(abs(tall(data , index , rank)) >= quan)
    larger=NULL
    for(i in 1:length(large)){
      if(i<length(large)){larger<-c(larger,paste(large[i],", "))}
      else{larger<-c(larger,paste(large[i]))}
    }
    larger<-paste(larger,collapse = "")
  }
  }

  else {
    HT1<-function(x){
      mx = apply(x,2,mean)
      vx = var(x)
      n = dim(x)[1]
      k = dim(x)[2]
      t2 = n*t(mx)%*%solve(vx)%*%mx
      t2[1,1]
    }
    HT2<-function(x){
      mx=apply(x,2,mean)
      vx=var(x)
      n=dim(x)[1]
      k=dim(x)[2]
      t2=n*t(mx)%*%solve(vx)%*%mx
      fs = (n-k)*t2/((n-1)*k)
      fs[1,1]
    }
    statistic <- function(x, fun){
      n = dim(x)[1]
      bin = rbinom(n,1,0.5)
      bin = 2*bin-1
      y = diag(bin)%*%x
      fun(y)
    }
    binary <- function(x, k){
      tmp = NULL
      y = x
      if(x < 2^k) {
        for(i in k-1:k){
          a = floor(y/2^i)
          tmp = c(tmp, a)
          y = y-a*2^i
        }
      }
      2*(tmp-0.5)
    }

    tall <- function(dat,index=NULL , rank=NULL){
      x=dat
      mx=apply(x,2,mean)
      vx=apply(x,2,var)
      n=dim(x)[1]
      t=mx*sqrt(n)/sqrt(vx)
      t
    }
    tmax <- function(dat,index=NULL , rank=NULL){   max(tall(dat))}
    tmaxabs <- function(dat,index=NULL , rank=NULL){  max(abs(tall(dat)))}

    zall<-function(x){
      mx=apply(abs(x),2,sum)/2
      vx=apply(x^2,2,sum)/4
      sr=apply(x*(x>0),2,sum)
      z=(sr-mx)/sqrt(vx)
      z
    }
    zmax<-function(x){ max(zall(x)) }
    zmaxabs<-function(x){ max(abs(zall(x))) }

    STAT = stat
    if(stat == "HT") stat = HT1
    else if(stat == "tmax") stat = tmax
    else if(stat == "tmaxabs") stat = tmaxabs
    else if(stat == "zmax") stat = zmax
    else if(stat == "zmaxabs") stat = zmaxabs


    m = nrow(data) / 2
    n = nrow(data) / 2
    k = ncol(data) - 2

    nc = ncol(data)
    uni_index = unique(data[,nc-1])
    uni_pair = unique(data[,nc])
    x = data[data[,nc-1] == uni_index[1],]
    y = data[data[,nc-1] == uni_index[2],]

    x = x[x[,nc] == uni_pair,]
    y = y[y[,nc] == uni_pair,]

    x = x[,1:k]
    y = y[,1:k]
    D = x - y

    rank=FALSE

    if(STAT == "zmax" | STAT == "zmaxabs"){
      rank=TRUE
      sd = sign(D)
      SRD = abs(D)
      for(i in 1:k){
        SRD[,i] = rank(SRD[,i])
      }
      D = SRD*sd
    }

    Fobs = stat(D)


    if(method_p == "exact"){
      l = dim(D)[1]
      ppm = NULL
      N = 2^l-1
      for(i in 0:N) {
        condition = diag(binary(i, l)) %*% D
        ppm = c(ppm, stat(condition))
      }
      p.value = length(ppm[ppm >= Fobs]) / N

    }

    else if(method_p == "sampling"){
      results <- replicate(samplenum, statistic(D,stat))
      p.value = length(results[results >= Fobs]) / samplenum
    }

    if(method_p == "asymptotic" & STAT == "HT"){
      f = HT2(D)
      p.value = 1 - pf(f,k,n-k)
    }
  }

  ############# output ################################
  names(p.value) = method_p
  names(Fobs) = STAT
  if(method_p == "asymptotic") {
    attr(p.value , "type") = switch(STAT,"HT"="F distribution","wsum"="normal")
  }

  alternative=if(STAT=="wsum"){"Each components in sample0 is greater than sample1."}else{"The means are different."}

  output <- list(method = if(pair){"Multiple Dimensional paired test"}else("Multiple Dimensional test") , score = if(rank){"Wilcoxon"}else{"original"} , stat =Fobs ,pval = p.value ,  alternative=alternative)


  if(diff == TRUE & STAT != "wsum"){
    output<-c(output,addition=paste(larger , 'dimensions are considered having significant differences\n'))
  }

  class(output) = "nonp"

  output
  ##################################

}

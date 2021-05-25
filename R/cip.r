cip <- function(x, conf.level = 0.95, p = 0.5){

  ################### warning ##################

  if(!is.vector(x)) stop("The input must be a vector")
  if(!is.numeric(x)) stop("The input must be numerical")

  cts_value = c(conf.level , p)
  cts_name = c("conf.level" , "p")
  cts = data.frame(cts_name , cts_value)
  for(i in 1:nrow(cts)){
    if(!is.numeric(cts[i,2])) stop(paste(cts[i,1] , "input is not true"))
    if(i %in% c(2,1)){
      if(cts[i,2] < 0 | cts[i,2] > 1) stop(paste(cts[i,1] , "input is out of range"))
    }
  }

  ###########################################

  n=length(x)
  y=sort(x)
  beta = 1-(1-conf.level)/2
  d=qnorm(beta)*sqrt(n*p*(1-p))
  a=round(p*n-d)
  b=round(p*n+1+d)

  if (a>=1 & b<=n){
    cat("The ",conf.level*100,"% confidence interval for the ",(p*100),"th percentile is ","(",y[floor(a)],",",y[ceiling(b)],").\n",sep = "")
    interval <- list(lower.rank=floor(a), upper.rank=ceiling(b), lower=y[floor(a)], upper=y[upper.rank=ceiling(b)])
  }

  else if (a<1){
    a <- 1
    b <- 1+n*p+sqrt(n*p*(1-p))*qnorm(conf.level+pnorm((a-n*p)/sqrt(n*p*(1-p))))
    warning("Can't to find the symmetric shortest interval because the resulting a(the order of lower bound) is less than 1")
    cat("The ",conf.level*100,"% confidence interval for the ",(p*100),"th percentile is ","(",y[floor(a)],",",y[ceiling(b)],").\n",sep = "")
    interval <- list(lower.rank=floor(a), upper.rank=ceiling(b), lower=y[floor(a)], upper=y[ceiling(b)])
  }

  else if (b>length(x)){
    b <- length(x)
    a <- n*p+sqrt(n*p*(1-p))*qnorm(pnorm((b-1-n*p)/sqrt(n*p*(1-p)))-conf.level)
    warning("Can't to find the symmetric shortest interval because the resulting b(the order of upper bound) is greater than length of given data")
    cat("The ",conf.level*100,"% confidence interval for the ",(p*100),"th percentile is ","(",y[floor(a)],",",y[ceiling(b)],").\n",sep = "")
    interval <- list(lower.rank=floor(a), upper.rank=ceiling(b), lower=y[floor(a)], upper=y[ceiling(b)])
  }
  invisible(interval)
}

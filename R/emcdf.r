emcdf <- function(x, conf.level=0.05){

  #################### warning #######################

  if(!is.vector(x)) stop("The input must be a vector")
  if(!is.numeric(x)) stop("The input must be numerical")
  if(!is.numeric(conf.level)) stop("The conf.level must be numerical")
  if(conf.level < 0 | conf.level >1) stop("The conf.level is out of range")

  ###########################################

  u <- sort(x)
  n=length(x)
  for(i in 1:length(x)){
    u[i] = sum(x<=u[i])/n
  }

  max_x <- max(x)
  min_x <- min(x)
  rangex <- max_x - min_x
  start <- min_x - rangex/n
  final <- max_x + rangex/n

  x1=c(start, sort(x), final)
  y1=c(0, u, 1)

  deta <- (2*u*n+qnorm(1-conf.level/2)^2)^2-4*(n+qnorm(1-conf.level/2)^2)*(u^2)*n
  lower <- ((2*u*n+qnorm(1-conf.level/2)^2) - sqrt(deta))/(2*(n+qnorm(1-conf.level/2)^2))
  upper <- ((2*u*n+qnorm(1-conf.level/2)^2) + sqrt(deta))/(2*(n+qnorm(1-conf.level/2)^2))

  y <- data.frame(sample=sort(x), empirical.cdf = u, lower = lower, upper = upper)
  class(y)<-c("emcdf","data.frame")
  return(y)
}

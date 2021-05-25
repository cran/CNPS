
RMD_test = function(x , y , alternative = "greater" , mu1=median(x) , mu2=median(y), method_p="exact" , samplenum = 2000 , samplemethod = "R" , conf.level.sample = 0.95 ){

  ################## warning message #############################

    if(!is.numeric(c(x,y))) stop("The input data is not a numeric")

    parameter_value = c(alternative,method_p,samplemethod)
    parameter_name = c("alternative","method_p","samplemethod")
    parameter_content = c("greater less two.sided" , "sampling exact" , "S R W")
    parameter = data.frame(parameter_value , parameter_name , parameter_content)

    for(i in 1:nrow(parameter)){
      if(!parameter[i,1] %in% strsplit(parameter[i,3]," ")[[1]]) stop(paste(parameter[i,2] , "input is not true"))
    }

    cts_value = c(mu1 , mu2 , samplenum , conf.level.sample)
    cts_name = c("mu1" , "mu2" , "samplenum" , "conf.level.sample")
    cts = data.frame(cts_name , cts_value)
    for(i in 1:nrow(cts)){
      if(!is.numeric(cts[i,2])) stop(paste(cts[i,1] , "input is not true"))
      if(i %in% c(4)){
        if(cts[i,2] < 0 | cts[i,2] > 1) stop(paste(cts[i,1] , "input is out of range"))
      }
    }

    if(method_p!="sampling" & !missing(samplenum)) warning("\n samplenum is not working")
    if(method_p!="sampling" & !missing(samplemethod)) warning("\n samplemethod is not working")
    if(method_p!="sampling" & !missing(conf.level.sample)) warning("\n conf.level.sample is not working")

  ################################################
    mean1.diff = function(x , nx){
      mean(x[1:nx]) / mean(x[-(1:nx)])
    }

    mean2.diff = function(x , nx){
      max(mean(x[1:nx]) , mean(x[-(1:nx)])) / min(mean(x[1:nx]) , mean(x[-(1:nx)]))
    }

    nx = length(x)
    ny = length(y)
    n = nx + ny
    dev_x = abs(x - mu1)
    dev_y = abs(y - mu2)
    data  <- c(dev_x, dev_y)
    index <- c(rep(1 , nx) , rep(2 , ny))

    if(alternative != "two.sided"){
      RMD_obs =  mean(dev_x) / mean(dev_y)

      if(method_p=="exact"){
        all.comb = combn(n, nx)
        N = dim(all.comb)[2]
        replace = rep(1, N)
        for(i in 1:N){
          ind  = all.comb[,i]
          replace[i] =  mean(data[ind]) / mean(data[-ind])
        }

        if(alternative == "greater"){
          p.value = length(replace[replace >= RMD_obs]) / N
        }
        else if(alternative == "less"){
          p.value = length(replace[replace <= RMD_obs]) / N
        }
      }

      else if(method_p=="sampling"){
        replace = apply(TwosampleSRS(data,index,samplenum,method = samplemethod),2, function(x) mean1.diff(x,nx) )
        if(alternative == "greater"){
          p.value = length(replace[replace >= RMD_obs]) / samplenum
          p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
          p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        }
        else if(alternative == "less"){
          p.value = length(replace[replace <= RMD_obs]) / samplenum
          p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
          p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        }
      }
    }

    else if(alternative == "two.sided"){
      RMD_obs =  max(mean(dev_x) , mean(dev_y)) / min(mean(dev_x) , mean(dev_y))
      if(method_p=="exact"){
        all.comb = combn(n, nx)
        N = dim(all.comb)[2]
        replace = rep(1, N)
        for(i in 1:N){
          ind  = all.comb[,i]
          replace[i] =  max(mean(data[ind]) , mean(data[-ind])) / min(mean(data[ind]) , mean(data[-ind]))
        }
          p.value = length(replace[replace >= RMD_obs]) / N
      }

      else if(method_p=="sampling"){
          replace = apply(TwosampleSRS(data,index,samplenum,method = samplemethod),2, function(x) mean2.diff(x,nx) )
          p.value = length(replace[replace >= RMD_obs]) / samplenum
          p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
          p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
      }
    }
    names(p.value) = method_p
    names(RMD_obs) = "RMDobs"

    conf=NULL
    if(method_p == "sampling"){
      attr(p.value , "type") = samplemethod
      conf = c(p.down , p.up)
      attr(conf , "conf.level") = conf.level.sample
    }
    null.value = c("The variance of the first" , "the variance of the second")
    attr(null.value , "direction") = alternative
    output = list(method = "RMD test" , stat = RMD_obs , conf.int=conf ,pval = p.value ,  null.value = null.value)

    class(output) = "nonp"
    output
}

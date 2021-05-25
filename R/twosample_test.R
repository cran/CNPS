
twosample_test = function(x , y , alternative = "greater" , score = "wilcoxon" , method_p = "sampling" , samplenum = 2000 ,samplemethod="R", conf.level.sample = 0.95 , conf.diff = TRUE, conf.level.diff = 0.95){

    ################ warning #########################
    if(!is.numeric(c(x,y))) stop("The input data is not a numeric")

    parameter_value = c(alternative,score,method_p,samplemethod,conf.diff)
    parameter_name = c("alternative","score","method_p","samplemethod","conf.diff")
    parameter_content = c("greater less two.sided" , "original wilcoxon van exp" , "sampling exact asymptotic" , "S R W" , "TRUE FALSE 1 0")
    parameter = data.frame(parameter_value , parameter_name , parameter_content)

    for(i in 1:nrow(parameter)){
      if(!parameter[i,1] %in% strsplit(parameter[i,3]," ")[[1]]) stop(paste(parameter[i,2] , "input is not true"))
    }

    cts_value = c(samplenum , conf.level.sample , conf.level.diff)
    cts_name = c("samplenum" , "conf.level.sample" , "conf.level.diff")
    cts = data.frame(cts_name , cts_value)
    for(i in 1:nrow(cts)){
      if(!is.numeric(cts[i,2])) stop(paste(cts[i,1] , "input is not true"))
      if(i %in% c(2,3)){
        if(cts[i,2] < 0 | cts[i,2] > 1) stop(paste(cts[i,1] , "input is out of range"))
      }
    }


    if(method_p != "sampling" & !missing(samplenum)) warning("\n samplenum is not working")
    if(method_p != "sampling" & !missing(samplemethod)) warning("\n samplemethod is not working")
    if(method_p != "sampling" & !missing(conf.level.sample)) warning("\n conf.level.sample is not working")
    if(conf.diff==FALSE & !missing(conf.level.diff)) warning("\n conf.level.diff is not working")


    ##############################################

    expsc<-function(x)
    {
      y=rank(x)
      n=length(x)
      sc=rep(0, n)
      for(i in 1:n)
      {
        for(j in 1:y[i])
        { sc[i] =  sc[i]+1/(n+1-j) }
      }
      sc
    }

    nx = length(x)
    ny = length(y)
    n = nx + ny
    data = c(x , y)
    index = c(rep(1 , nx) , rep(2 , ny))

    if(conf.diff){
      diff = rep(1, nx*ny)
      for (i in 1: nx){
        for(j in 1: ny){
          k = (i-1)*ny+j
          diff[k] = x[i] - y[j]
        }
      }
      diff = sort(diff)
      mu = nx*ny/2
      va = mu*(n+1)/6
      ka = round(mu - qnorm(0.5+0.5*conf.level.diff)*sqrt(va))
      kb = round(mu + qnorm(0.5+0.5*conf.level.diff)*sqrt(va) + 1)
      diff.up = diff[kb]
      diff.down = diff[ka]
      Hodges = median(diff)
    }


    if(score == "original"){
    }
    else if(score == "wilcoxon"){
      data = rank(data)
    }
    else if(score == "van"){
      data = qnorm( rank(data)/(n+1))
    }
    else if(score == "exp"){
      data = expsc(data)
    }

    Dobs = sum(data[index==1])


    if(method_p == "exact"){
      all.comb = combn(n , nx)
      N = choose(n , nx)
      replace = rep(1 , N)
      for(i in 1:N){
        ind = all.comb[,i]
        replace[i] = sum(data[ind])
      }

      if(alternative == "greater"){
        p.value = length(replace[replace >= Dobs]) / N
      }
      else if(alternative == "less"){
        p.value = length(replace[replace <= Dobs]) / N
      }
      else if(alternative == "two.sided"){
        p.greater = length(replace[replace >= Dobs]) / N
        p.less = length(replace[replace <= Dobs]) / N
        p.value = min(2*p.less , 2*p.greater,1)
      }
    }

    else if(method_p == "sampling"){
      replace = apply(TwosampleSRS(data,index,samplenum,method = samplemethod),2, function(x) sum(x[1:nx]))
      if(alternative == "greater"){
        p.value = length(replace[replace >= Dobs]) / samplenum
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        }
      else if(alternative == "less"){
        p.value = length(replace[replace <= Dobs]) / samplenum
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
      }
      else if(alternative == "two.sided"){
        p.greater = length(replace[replace >= Dobs]) / samplenum
        p.less = length(replace[replace <= Dobs]) / samplenum
        p.value = min(2*p.less , 2*p.greater,1)
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
      }
    }

    else if(method_p == "asymptotic"){
      E = nx * mean(data)
      V = nx * ny * var(data) / n
      if(alternative == "greater"){
        p.value = 1 - pnorm((Dobs - E - 0.5)/sqrt(V))
      }
      else if(alternative == "less"){
        p.value = pnorm((Dobs - E + 0.5)/sqrt(V))
      }
      else if(alternative == "two.sided"){
        p.value = min(2 - 2 * pnorm((abs(Dobs - E) - 0.5)/sqrt(V)),1)
      }
    }


    ############# output ################################

    names(p.value) = method_p
    names(Dobs) = "Dobs"
    if(method_p == "asymptotic") attr(p.value , "type") = "normal"

    conf=NULL
    if(method_p == "sampling"){
      attr(p.value , "type") = samplemethod
      conf = c(p.down , p.up)
      attr(conf , "conf.level") = conf.level.sample
    }
    null.value = c("The mean of the first" , "the mean of the second")
    attr(null.value , "direction") = alternative
    output = list(method = "Two sample test" , score = score , stat = Dobs , conf.int=conf ,pval = p.value ,  null.value = null.value)


    if(conf.diff){
      output<-c(output,addition=paste("The Hodges-Lehmann statistic =" , Hodges , "\nThe" , conf.level.diff*100 , "% CI for mean difference is [" , diff.down , "," , diff.up , "]"))
    }

    class(output) = "nonp"
    output
    ##################################

}

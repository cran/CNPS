
ksample_test <- function(x , group , score = "kruskal", method_p = "sampling" , type = "normal" , samplenum = 1000 , conf.level.sample = 0.95){

    ################ warning message #########################
    if(!is.numeric(x)) stop("The input data is not a numeric")
    if(length(x) != length(group)) stop("The length of the data is not equal to the length of the group")

    parameter_value = c(score , method_p , type)
    parameter_name = c("score" , "method_p" ,"type")
    parameter_content = c("original kruskal van exp" , "sampling asymptotic exact" , "normal JT")
    parameter = data.frame(parameter_value , parameter_name , parameter_content)

    for(i in 1:nrow(parameter)){
      if(!parameter[i,1] %in% strsplit(parameter[i,3]," ")[[1]]) stop(paste(parameter[i,2] , "input is not true"))
    }

    cts_value = c(samplenum , conf.level.sample)
    cts_name = c("samplenum" , "conf.level.sample")
    cts = data.frame(cts_name , cts_value)
    for(i in 1:nrow(cts)){
      if(!is.numeric(cts[i,2])) stop(paste(cts[i,1] , "input is not true"))
      if(i %in% c(2)){
        if(cts[i,2] < 0 | cts[i,2] > 1) stop(paste(cts[i,1] , "input is out of range"))
      }
    }

    if(method_p != "sampling" & !missing(samplenum)) warning("\n samplenum is not working")
    if(method_p != "sampling" & !missing(conf.level.sample)) warning("\n conf.level.sample is not working")
    if(type == "JT" & !missing(score)) warning("\n score is not working")

    ##############################################
    expsc<-function(x){
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

    kw = function(x , group){
      group = factor(group)
      m = mean(x)
      sum(tapply(x , group , function(x) length(x) * (mean(x) - m)^2)) / var(x)
    }

    mw<-function(x1,x2){
      u = 0
      for(i in 1:length(x2)){
        u = u + sum(x1 > x2[i])
      }
      u
    }

    jt<-function(x, group){
      ug = unique(group)
      k = length(ug)
      jt = 0
      for(j in 1:(k-1)){
        for(i in (j+1):k){
          xa = x[group==ug[i]]
          xb = x[group==ug[j]]
          tmp = mw(xa,xb)
          jt = jt + tmp
        }
      }
      jt
    }

    n = length(x)
    k = length(unique(group))
    group = factor(group)

    if(type == "normal"){
      if(score == "original"){
        method_name = "F-test"
      }
      else if(score == "kruskal"){
        x = rank(x)
        method_name = "Kruskal-Wallis rank sum test"
      }
      else if(score == "van"){
        x = qnorm( rank(x)/(n+1))
        method_name = "Van der Waerden test"
      }
      else if(score == "exp"){
        x = expsc(x)
        method_name = "Savage exponential test"
      }

      Dobs = kw(x , group)


      if(method_p == "asymptotic"){
        p.value = 1 - pchisq(Dobs , k-1)
        C = pchisq(Dobs , k-1)
      }
      else if(method_p == "sampling"){
        replace = replicate(samplenum ,kw(sample(x) , group))
        p.value = length(replace[replace >= Dobs]) / samplenum
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
      }
      else if(method_p=="exact"){
        print("under preparing...")
      }
    }

    else if(type == "JT"){
      Dobs = jt(x,group)


      if(method_p == "asymptotic"){
        size = table(group)
        mu = (n^2-sum(size^2))/4
        sigma = (n^2*(2*n+3)-sum(size^2*(2*size+3)))/72
        p.value = 1 - pnorm((Dobs-mu)/sqrt(sigma))

      }
      else if(method_p == "sampling"){
        replace = replicate(samplenum, jt(sample(x), group))
        p.value = length(replace[replace >= Dobs]) / samplenum
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5

      }
      else if(method_p=="exact"){
        print("under preparing...")
      }
    }
    ############# output ################################

    method=paste(switch(type,"normal"=method_name,"JT"="Jonckheere-Terpstra test"),"with",k,"groups")

    names(p.value) = method_p
    names(Dobs) = "Dobs"
    if(method_p == "asymptotic") {
      attr(p.value , "type") = switch(type,"normal"="chisq","JT"="normal")
      if(type=="normal") attr(p.value , "df")=k-1
    }
    conf=NULL
    if(method_p == "sampling"){
      #attr(p.value , "type") = samplemethod
      conf = c(p.down , p.up)
      attr(conf , "conf.level") = conf.level.sample
    }
    alternative=switch(type,"normal"="There exists at least one pair of groups such that the means differ.","JT"="The grouped means are ordered.")
    output = list(method = method  , stat = Dobs , conf.int=conf ,pval = p.value ,  alternative=alternative)

    class(output) = "nonp"

    output
    ##################################

}

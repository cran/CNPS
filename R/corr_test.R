
corr_test = function(x , y , alternative = "greater" , measure = "pearson" , method_p = "sampling" , samplenum = 1000 , conf.level.sample = 0.95){


  ################ warning message #########################
  if(!is.numeric(c(x,y))) stop("The input data is not a numeric")
  if(length(x) != length(y)) stop("The input pair of datas have different length")

  parameter_value = c(alternative,measure,method_p)
  parameter_name = c("alternative","measure","method_p")
  parameter_content = c("greater less two.sided" , "pearson spearman kendall" , "sampling asymptotic exact")
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

  ##############################################
  kendtau<-function(x,y){
    count = 0
    n = length(x)
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        tmp = ((x[i]-x[j])*(y[i]-y[j]) >0)
        count=count+tmp
      }
    }
    tau = count / choose(n,2)
    2*tau - 1
  }

    nx = length(x)
    ny = length(y)
    n = nx
    if(measure == "pearson"){
      Dobs = cor(x , y)
    }
    else if(measure == "spearman"){
      x = rank(x) ; y = rank(y)
      Dobs = cor(x , y)
    }
    else if(measure == "kendall"){
      Dobs = kendtau(x , y)
    }

    if(method_p == "sampling"){
      if(measure != "kendall")  replace = replicate(samplenum , cor(x , sample(y)))
      else  replace = replicate(samplenum , kendtau(x , sample(y)))

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
         replace = abs(replace)
         p.value = length(replace[replace >= abs(Dobs)]) / samplenum
         p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
         p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
       }
    }
    else if(method_p == "asymptotic"){
      if(measure != "kendall")  z = Dobs * sqrt(n-1)
      else {
        V = (4*n+10)/(9*(n^2-n))
        z = Dobs / sqrt(V)
      }
       if(alternative == "greater"){
         p.value = 1 - pnorm(z)
       }
       else if(alternative == "less"){
         p.value = pnorm(z)
       }
       else if(alternative == "two.sided"){
         p.value = min(1,2 - 2 * pnorm(abs(z)))
       }
    }
    else if(method_p == "exact"){
      if(!requireNamespace("e1071")){stop("Need R-Package e1071 to finish test.")}
      if(measure != "kendall"){
        z <- e1071::permutations(ny)
        N <- nrow(z)
        replace <- apply(z,MARGIN = 1,function(z){cor(x , y[z])})
      }
      else{
        z <- e1071::permutations(ny)
        N <- nrow(z)
        replace <- apply(z,MARGIN = 1,function(z){kendtau(x , y[z])})
      }

      if(alternative == "greater"){
        p.value = length(replace[replace >= Dobs]) / length(replace)
      }
      else if(alternative == "less"){
        p.value = length(replace[replace <= Dobs]) / length(replace)
      }
      else if(alternative == "two.sided"){
        p.greater = length(replace[replace >= Dobs]) / length(replace)
        p.less = length(replace[replace <= Dobs]) / length(replace)
        p.value = min(2*p.greater , 2*p.less,1)
      }
    }
    method="Correlation test"
    names(p.value) = method_p
    names(Dobs) = "Dobs"
    if(method_p == "asymptotic") {
      attr(p.value , "type") <- "normal"
    }
    conf=NULL
    if(method_p == "sampling"){
      conf = c(p.down , p.up)
      attr(conf , "conf.level") = conf.level.sample
    }
    null.value = c(paste(measure,"correlation") , "0")
    attr(null.value , "direction") = alternative
    output = list(method = method  ,score=measure, stat = Dobs , conf.int=conf ,pval = p.value ,  null.value=null.value)

    class(output) = "nonp"
    output
}


pairwise_test = function(x , y , alternative = "greater" , score = "wilcoxon" , method_p = "asymptotic", method_asymptotic = "norm", method_wilcoxon = "type1"  , samplenum = 1000 , conf.level.sample = 0.95,samplemethod="R"){

    ################ warning message #########################
    if(!is.numeric(c(x,y))) stop("The input data is not a numeric")
    if(length(x) != length(y)) stop("The input pair of datas have different length")

    parameter_value = c(alternative,score,method_p,method_asymptotic, method_wilcoxon,samplemethod)
    parameter_name = c("alternative","score","method_p","method_asymptotic","method_wilcoxon","samplemethod")
    parameter_content = c("greater less two.sided" , "original wilcoxon sign" , "sampling exact asymptotic" , "norm binomal" , "type1 type2","R S")
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

    if(score == "sign" & method_p != "asymptotic") stop(" 'score = sign' only use asymptotic method to calculate P.value")
    if(score != "sign" & method_asymptotic == "binomal") stop("binomal asymptotic method can only be used when 'score = sign' ")


    if(score != "wilcoxon" & !missing(method_wilcoxon)) warning("\n method_wilcoxon is not working")
    if(method_p != "sampling" & !missing(samplenum)) warning("\n samplenum is not working")
    if(method_p != "sampling" & !missing(conf.level.sample)) warning("\n conf.level.sample is not working")


    ##############################################


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
    ppm <- function(x){
      k = length(x)
      ppm = NULL
      N = 2^k-1
      for(i in 0:N) {
        condition = x*binary(i, k)
        ppm = c(ppm, sum(condition[condition > 0]))
      }
      ppm
    }
    ppmrr<-function(x){
      ind = rbinom(length(x), 1, 0.5)
      ind = (ind-0.5)*2
      condition = x*ind
      sum(condition[condition > 0])
    }
    ppmr<-function(x,samplenum){
      N<-length(x)
      s_index<-sample(1:2^N,samplenum)-1
      ppmr<-NULL
      for(i in 1:samplenum){
        tmp=x*binary(s_index[i],N)
        ppmr=c(ppmr,sum(tmp[tmp>0]))
      }
      ppmr
    }
    sign_binomal_greater = function(k , n){
      total = 0
      for(i in k:n){
        total = total + 0.5^n * choose(n , i)
      }
      total
    }
    sign_binomal_less = function(k , n){
      total = 0
      for(i in 0:k){
        total = total + 0.5^n * choose(n , i)
      }
      total
    }



    data = x - y

    if(score == "original"){
    E = sum(abs(data)) / 2
    V = sum(abs(data) ^ 2) / 4
  }
    else if(score == "wilcoxon"){
    if(method_wilcoxon == "type1"){
      transit = data
      data = rank(abs(data)) * sign(data)
      data[transit == 0] = 0
    }
    else if(method_wilcoxon == "type2"){
      data = data[data != 0]
      data = rank(abs(data)) * sign(data)
    }
    n = length(data)
    E = n*(n+1)/4
    V = E*(2*n+1)/6
  }
    else if(score == "sign"){
      data[data > 0] = 1
      data[data < 0] = -1
      n = length(data)
      E = n / 2 ; V = n / 4
    }

    Dobs = sum(data[data > 0])


    if(method_p == "exact"){
      replace = ppm(data)

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

    else if(method_p == "sampling"){
      if(samplemethod=="R")
        replace = replicate(samplenum , ppmrr(data))
      else if(samplemethod=="S"){
        if(samplenum>2^length(data)){stop(paste("Samplenum of sampling without replacement should be smaller than",2^length(data)))}
        replace = ppmr(data,samplenum)
      }
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
        p.value = min(2*p.greater , 2*p.less,1)
        p.up = p.value + qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
        p.down = p.value - qnorm(0.5+0.5*conf.level.sample) * (p.value*(1 - p.value) / samplenum)^0.5
      }
      }

    else if(method_p == "asymptotic"){
      if(method_asymptotic == "norm"){
        if(alternative == "greater"){
          p.value = 1 - pnorm((Dobs - E - 0.5)/sqrt(V))
        }
        else if(alternative == "less"){
          p.value = pnorm((Dobs - E + 0.5)/sqrt(V))
        }
        else if(alternative == "two.sided"){
          p.value = min(1,2 - 2 * pnorm((abs(Dobs - E) - 0.5)/sqrt(V)))
        }
      }

      else if(method_asymptotic == "binomal"){
        if(alternative == "greater"){
          p.value = sign_binomal_greater(Dobs , n)
        }
        else if(alternative == "less"){
          p.value = sign_binomal_less(Dobs , n)
        }
        else if(alternative == "two.sided"){
          p.value = 2 * min(sign_binomal_greater(Dobs , n) , sign_binomal_less(Dobs , n))
        }
      }
    }

    ############# output ################################
    names(p.value) = method_p
    names(Dobs) = "Dobs"
    if(method_p == "asymptotic") attr(p.value , "type") = method_asymptotic

    conf=NULL
    if(method_p == "sampling"){
      attr(p.value , "type") = samplemethod
      conf = c(p.down , p.up)
      attr(conf , "conf.level") = conf.level.sample
    }
    null.value = c("The mean of the first" , "the mean of the second")
    attr(null.value , "direction") = alternative
    output = list(method = "Pairwise comparision test" , score = score , stat = Dobs , conf.int=conf ,pval = p.value ,  null.value = null.value)

    class(output) = "nonp"

    output
    ##################################

}




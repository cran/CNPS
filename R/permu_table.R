
permu_table = function(data , permu = "row" , row = NULL , col = NULL , fix = "row" , samplenum = 1000){

  ########### warning #####################
  if(!is.matrix(data) & !is.data.frame(data)) stop("The input data must be a matrix or dataframe")
  if(!is.numeric(samplenum)) stop("The samplenum must be a numeric")

  parameter_value = c(permu , fix)
  parameter_name = c("permu" , "fix")
  parameter_content = c("row col" , "row col")
  parameter = data.frame(parameter_value , parameter_name , parameter_content)

  for(i in 1:nrow(parameter)){
    if(!parameter[i,1] %in% strsplit(parameter[i,3]," ")[[1]]) stop(paste(parameter[i,2] , "input is not true"))
  }

  if(!missing(row)){
    if(!is.numeric(row)) stop("The row parameter must be numeric")
    if(!is.vector(row)) stop("The row parameter must be a vector")
    if(length(row) != nrow(data)) stop("The length of 'row' must equal to the row length of data")
  }
  if(!missing(col)){
    if(!is.numeric(col)) stop("The col parameter must be numeric")
    if(!is.vector(col)) stop("The col parameter must be a vector")
    if(length(col) != ncol(data)) stop("The length of 'col' must equal to the col length of data")
  }

  if((missing(row) | missing(col)) & !missing(fix)) warning("\n fix is not working")

  ###################################



  perm0 <- function(row,col){
    n = sum(col)
    r = length(row)
    c = length(col)
    x = NULL
    for(i in 1:c){
      x = c(x, rep(i,col[i]))
    }
    y = sample(x)

    freq = matrix(rep(0, r*c), r,c)
    count = 0
    for(j in 1:r){
      tmpy = y[(count+1):(count+row[j])]
      for(i in 1:c){
        freq[j,i] = sum(tmpy==i)
        }
      count = count+row[j]
    }
    freq
  }
  perm1 <- function(row,col){
    n = sum(col)
    r = length(row)
    c = length(col)
    x = NULL
    for(i in 1:r){
      x=c(x, rep(i,row[i]))
    }
    y = sample(x)

    freq = matrix(rep(0, r*c), r,c)
    count = 0
    for(j in 1:c){
      tmpy = y[(count+1):(count+col[j])]
      for(i in 1:r){
        freq[i,j] = sum(tmpy==i)
        }
      count = count+col[j]
    }
    freq
  }
  chi2 <- function(freq){
    row = apply(freq,1,sum)
    col = apply(freq,2,sum)
    row = matrix(rep(row , ncol(freq)) , ncol = ncol(freq))
    col = matrix(rep(col , nrow(freq)) , ncol = ncol(freq) , byrow = TRUE)
    fenmu = row * col / sum(freq)
    sum(freq^2 / fenmu) - sum(freq)
  }

  if(is.null(row) & is.null(col)){
    row = apply(data,1,sum)
    col = apply(data,2,sum)

    obs = chi2(data)
    if(permu == "row"){
      results <- replicate(samplenum,  chi2(perm0(row,col)))
    }
    else if(permu == "col"){
      results <- replicate(samplenum,  chi2(perm1(row,col)))
    }
    p.value = sum(results >= obs) / samplenum
    method="Pearson's Chi-squared test"
    names(p.value) = "simulated"
    names(obs) = "chi-squared_obs"
    alternative="There exists at least one pair of groups that are not independent."
    output = list(method = method  , stat = obs ,pval = p.value ,  alternative=alternative)
    class(output) = "nonp"
  }
  else {

    r = nrow(data)
    c = ncol(data)
    if(is.vector(row) & is.null(col)){
      scoring = NULL
      group = NULL
      for(i in 1:r){
        for(j in 1:c){
          scoring = c(scoring , rep(row[i] , data[i,j]))
          group = c(group , rep(j , data[i,j]))
        }
      }
      if(c > 2) {
        output<-ksample_test(scoring , group , score = "kruskal" , samplenum = samplenum, type = "normal" )
        output$method<-"Kruskal test on table"
        output$alternative<-"The rows are ordered."
        }
      else if(c==2) {
        output<-twosample_test(scoring[group == 1] , scoring[group == 2] ,score = "wilcoxon" , alternative = "two.sided", method_p = "sampling" , samplenum = samplenum)
        output$method<-"Wilcoxon test on table"
        output$alternative<-"The rows are ordered."
        }
    }
    else if(is.vector(col) & is.null(row)){
      scoring = NULL
      group = NULL
      for(j in 1:c){
        for(i in 1:r){
          scoring = c(scoring , rep(col[j] , data[i,j]))
          group = c(group , rep(i , data[i,j]))
        }
      }
      if(c > 2) {
        output<-ksample_test(scoring , group , score = "kruskal" , samplenum = samplenum, type = "normal" )
        output$method<-"Kruskal test on table"
        output$alternative<-"The columns are ordered."
        }
      else if(c==2) {
        output<-twosample_test(scoring[group == 1] , scoring[group == 2] ,score = "wilcoxon" , alternative = "two.sided", method_p = "sampling" , samplenum = samplenum)
        output$method<-"Wilcoxon test on table"
        output$alternative<-"The columns are ordered."
        }

    }
    else if(is.vector(col) & is.vector(row)){
      if(fix == "row"){
        scoring = NULL
        group = NULL
        for(j in 1:c){
          for(i in 1:r){
            scoring = c(scoring , rep(col[j] , data[i,j]))
            group = c(group , rep(i , data[i,j]))
          }
        }
      }
      else if(fix == "col"){
        scoring = NULL
        group = NULL
        for(i in 1:r){
          for(j in 1:c){
            scoring = c(scoring , rep(row[i] , data[i,j]))
            group = c(group , rep(j , data[i,j]))
          }
        }
      }


      output<-ksample_test(scoring , group , type = "JT",samplenum = samplenum)
      output$method<-"Jonckheere-Terpstra test on table"
      output$alternative<-"Both rows and columns are ordered."

    }
  }
  return(output)
}

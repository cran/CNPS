print.nonp<-function(x, digits = getOption("digits"), prefix = "\t\t",...){
  cat("\n")
  theme<-strwrap(x$method, prefix = prefix)
  if(!is.null(x$score)){
    if(x$score %in% c("pearson", "spearman", "kendall"))
      theme<-c(theme,paste("with",x$score,"measurement"))
    else
      theme<-c(theme,paste("using",x$score,"scoring"))
  }
  cat(theme,"\n")

  out <- character()
  if (!is.null(x$stat))
    out <- c(out, paste(names(x$stat), "=", format(x$stat, digits = max(1L, digits - 2L))))
  if (!is.null(x$pval)) {
    fp <- paste(format.pval(x$pval, digits = max(1L, digits - 3L)))
    if(!is.null(attr(x$pval,"type"))){
      if(attr(x$pval,"type")%in%c("W","R","S"))
        out <- c(out,paste(names(x$pval),switch(attr(x$pval,"type"),"W"="without replacement","S"="without replacement","R"="with replacement"),paste("method to calculate :\n\t",if(!is.null(attr(x$pval,"df"))){paste("df = ",attr(x$pval,"df"),",")}, "p-value = ",sep=""), if (startsWith(fp, "<")) fp else paste(fp)))
      else
        out <- c(out,paste(attr(x$pval,"type"),names(x$pval),paste("method to calculate :\n\t", if(!is.null(attr(x$pval,"df"))){paste("df = ",attr(x$pval,"df"),",")},"p-value = ",sep=""),if (startsWith(fp, "<")) fp else paste (fp)))
    }
    else
      out <- c(out,paste(names(x$pval),paste("method to calculate :\n\t", if(!is.null(attr(x$pval,"df"))){paste("df = ",attr(x$pval,"df"),",")}, "p-value = ",sep=""), if (startsWith(fp, "<")) fp else paste(fp)))
  }
  cat(paste(out, collapse = "\n"),sep = "\n")

  if (is.null(x$alternative)) {
    cat("alternative hypothesis: \n\t")
    alt.char <- switch(attr(x$null.value,"direction"), two.sided = "not equal to", less = "less than", greater = "greater than")
    cat(x$null.value[1], " is ", alt.char, " ", x$null.value[2], "\n", sep = "")
  }
  else if(is.null(x$null.value)) cat("alternative hypothesis: \n\t",x$alternative, "\n", sep = "")
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")),
        "% confidence interval of p-value :\n\t",
        paste("[",format(x$conf.int[1], digits = digits-3L), ",",format(x$conf.int[2], digits = digits-3L), "]",collapse = " "),
        "\n", sep = "")
  }

  if(!is.null(x$addition)) cat("\n",x$addition)

  if(!is.null(x$alternative) & !is.null(x$null.value)){warning("Alternative and null.value both exist, so output from null.value disabled.")}
  cat("\n")
  invisible(x)
}


plot.emcdf<-function(x,...){
  data<-x
  x<-data$sample
  n=length(x)
  max_x <- max(x)
  min_x <- min(x)
  rangex <- max_x - min_x
  start <- min_x - rangex/n
  final <- max_x + rangex/n

  x1=c(start, sort(x), final)
  y1=c(0, data$empirical.cdf, 1)
  lower_plot <- c(0,data$lower,1)
  upper_plot <- c(0,data$upper,1)
  plot(x1,y1, type='S', xlab='Cycles', ylab='Probability', main='Empirical Distribution',...)
  lines(x1, upper_plot, lty=2,type='S', col='red')
  lines(x1, lower_plot, lty=2,type='S', col='blue')
  legend("bottomright", lty=c(1,2,2), col=c('black', 'red', 'blue'), c('EM-CDF', 'Upper', 'Lower'), cex=0.8)
}

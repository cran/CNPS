\name{ksample_test}
\alias{ksample_test}

\title{
Multiple sample permutation test
}
\description{
Test whether there is a difference among k treatments.
}
\usage{
ksample_test(x, group, score = "kruskal", method_p = "sampling", type = "normal",
samplenum = 1000, conf.level.sample = 0.95)
}

\arguments{
  \item{x}{
numeric vector of data values
}
  \item{group}{
factor that determines the grouping of elements in x
}
  \item{score}{
a discrete value indicating the type of score. There are "original", "Wilcoxon", "van" and "exp" to be selected
}
  \item{method_p}{
a string indicating what method to use for p-value. "sampling" represents sampling; "asymptotic" represents using large sample approximations
}
  \item{type}{
"normal" refers to ordinary test, "JT" refers to ordered alternative hypothesis
}
  \item{samplenum}{
the number of SRS samples
}
  \item{conf.level.sample}{
p-value confidence level for SRS sampling
}
}
\details{
Use a permutation sample based on the F-statistic or use a large sample approximation to determine if there is a difference between the populations (treatments). If the sample is ordered, you can also use the JT test. The argument "score" allows you to choose different scoring system to do the test. If you want to use a special scoring system defined by yourself, just transform the data first and then choose score="original" in the function.
}
\value{
A list with following components
  \item{method}{
the test uesd
}
  \item{stat}{
the statistic of the data under the given scoring system
}
  \item{conf.int}{
the confidence interval for p-value(only if method_p = "sampling")
}
\item{pval}{
p-value for the test
}
\item{alternative}{
a character string describing the alternative hypothesis
}
}
\references{
Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.
}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
x1=c( 13.0, 24.1, 11.7, 16.3, 15.5, 24.5)
x2=c( 42.0, 18.0, 14.0, 36.0, 11.6, 19.0)
x3=c( 15.6, 23.8, 24.4, 24.0, 21.0, 21.1)
x4=c( 35.3, 22.5, 16.9, 25.0, 23.1, 26.0)
x <- c(x1, x2, x3, x4)
ind=c(rep(1,length(x1)), rep(2, length(x2)), rep(3, length(x3)), rep(4, length(x4)))
group=as.factor(ind)
ksample_test(x , group , type = "JT" , samplenum = 4000)
}


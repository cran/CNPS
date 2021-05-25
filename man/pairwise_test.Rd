\name{pairwise_test}
\alias{pairwise_test}

\title{
Paired Comparisons
}
\description{
Detects differences between two related samples.
}
\usage{
pairwise_test(x, y, alternative = "greater", score = "wilcoxon", method_p = "asymptotic",
method_asymptotic = "norm", method_wilcoxon = "type1", samplenum = 1000,
conf.level.sample = 0.95, samplemethod = "R")
}

\arguments{
  \item{x}{
numeric vectors of data values and should have the same length
}
  \item{y}{
numeric vectors of data values and should have the same length
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided", "greater"(default) or "less"
}
  \item{score}{
determines scoring systems and must be one of "original", "wilcoxon" or "sign"
}
  \item{method_p}{
a string indicating what method to use for p-value. "sampling" represents sampling; "asymptotic" represents using large sample approximations; "exact" represents Iterate through all combinations
}
  \item{method_asymptotic}{
determines the asymptotic distribution and should be one of "norm" or "binomial"(only for method_p="sign")
}
  \item{method_wilcoxon}{
indicates the way to compute wilcoxon ranks when the ties are 0 and could be one of "type1" or "type2"
}
  \item{samplenum}{
the number of SRS samples
}
  \item{conf.level.sample}{
p-value confidence level for SRS sampling
}
  \item{samplemethod}{
a discrete value indicating the method of sampling. "S" represents sample function sampling; "R" represents Put-back sampling
}
}
\details{
If score="sign", then method_p must be "asymptotic". Three scoring systems can use the normal approximation but only "sign" can use binomial approximation. Namely, the argument method_asymptotic can be selected as "binomial" only if method_p="sign". And method_wilcoxon indicates the method to deal with ties. "type1" means ranking with zeros and "type2" means ranking without zeros.
}
\value{
A list with following components
  \item{method}{
the test uesd
}
  \item{score}{
the score which is used
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
\item{null.value}{
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
x1=c(1530, 2130,2940,1960,2270)
x2=c(1290, 2250,2430,1900,2120)
pairwise_test(x1 , x2)
pairwise_test(x1 , x2 , method_p = "sampling" , samplenum = 4000)
pairwise_test(x1 , x2 , method_p = "asymptotic" , method_asymptotic = "norm")
}



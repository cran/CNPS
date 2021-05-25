\name{corr_test}
\alias{corr_test}

\title{
Correlation test
}
\description{
Test the correlation coefficient of the sample.
}
\usage{
corr_test(x, y, alternative = "greater", measure = "pearson",
method_p = "sampling", samplenum = 1000, conf.level.sample = 0.95)
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
  \item{measure}{
the way to measure the correlation coefficient and must be one of "pearson", "spearman" or "kendall"
}
  \item{method_p}{
a string indicating what method to use for p-value. "sampling" represents sampling; "asymptotic" represents using large sample approximations
}
  \item{samplenum}{
the number of SRS samples
}
  \item{conf.level.sample}{
p-value confidence level for SRS sampling
}
}
\details{
All procedures and methods of the correlation coefficient test based on the Spearman Correlation Coefficient are the same as for the Pearson Correlation Coefficient. But pay attention to that the correlation coefficient test based on Kendall Correlation Coefficient is a little different from the above two due to its definition.
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
x=c(68,70,71,72)
y=c(153,155,140,180)
corr_test(x , y , measure = "kendall" , method = "asymptotic")
corr_test(x , y , measure = "kendall" , method = "sampling")
}



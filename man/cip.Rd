\name{cip}
\alias{cip}

\title{
confidence interval for percentiles in the one-sample case
}
\description{
Finding confidence interval for (100p)th percentile in the one-sample case.
}
\usage{
cip(x, conf.level = 0.95, p = 0.5)
}
\arguments{
  \item{x}{
numeric vector of data values
}
  \item{conf.level}{
confidence level for the returned confidence interval
}
  \item{p}{
an arbitrary value from 0 to 1 which indicates the percentile
}
}
\details{
Usually we take the set with the shortest interval. But if the upper bound is greater than the maximum value(or the lower bound is less than the minimum value), we will choose the maximum value as the upper bound(the minimum value as the lower bound).
}
\value{
A list with following components
  \item{Lower.rank}{
rank of the lower bound of the confidence interval in the order statistic
}
  \item{Upper.rank}{
rank of the upper bound of the confidence interval in the order statistic
}
  \item{Lower}{
the lower bound of the confidence interval
}
\item{Upper}{
the upper bound of the confidence interval
}
}
\references{
Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.
}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
x <- c(72.1, 72.8, 72.9, 73.3, 76.1, 76.5, 78.8, 78.9, 79.7, 80.3, 80.5, 81.0)
cip(x)
cip(x, conf.level =0.9, p=0.7)
}



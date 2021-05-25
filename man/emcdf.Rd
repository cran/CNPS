\name{emcdf}
\alias{emcdf}

\title{
Estimating the population cdf
}
\description{
Finding confidence interval for the population cdf.
}
\usage{
emcdf(x, conf.level = 0.05)
}

\arguments{
  \item{x}{
numeric vector of data values
}
  \item{conf.level}{
confidence level for the returned confidence interval
}
}
\details{
This "emcdf" constructs the approximation interval according to the central limit theorem. And use "plot(emcdf(data))" will help us draw a plot conveniently.
}
\value{
A list with following components
  \item{sample}{
the given vector
}
  \item{empirical.cdf}{
the value of the empirical cdf
}
  \item{Lower}{
the lower bound of the confidence interval of the empirical cdf
}
\item{Upper}{
the upper bound of the confidence interval of the empirical cdf
}
}
\references{
Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.
}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
x <- c(7,11,15, 16, 20, 22, 24, 25, 29, 33, 34, 37, 41, 42, 49, 57, 66, 71, 84, 90)
em <- emcdf(x)
plot(em)
}



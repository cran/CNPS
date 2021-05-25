\name{siegel_tukey}
\alias{siegel_tukey}

\title{Siegel-Tukey Test

}
\description{Performs two-sample Siegel-Tukey test on vectors of data.

}
\usage{

siegel_tukey (x,y,adjust.median=FALSE,...)

}

\arguments{
  \item{x}{numeric vector of data values.

}
  \item{y}{numeric vector of data values.

}
  \item{adjust.median}{a logical indicating whether you want the adjusted median, which represents difference of each value with median.

}
  \item{...}{The input from twosample_test function.

}
}
\details{
The test is about the scale parameter, to test the deviance. The arguments you want to modify is almost same with \code{\link{twosample_test}}.
}
\value{

  \item{method}{ the test used.}
  \item{stat}{ the statistic of the original data.}
  \item{conf.int }{ the confidence interval for p-value(only if method_p = "sampling")}
  \item{pval}{ p-value for the test}
  \item{alternative}{ a character string describing the alternative hypothesis.

}

}
\references{Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.

}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
## A simple example
x <- c(33, 62, 84, 85, 88, 93, 97)
y<-c(4, 16, 48, 51, 66, 98)
siegel_tukey(x,y,adjust.median=FALSE)

}


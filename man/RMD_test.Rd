\name{RMD_test}
\alias{RMD_test}

\title{RMD Test

}
\description{Perform two-sample RMD test on vectors of data.

}
\usage{

RMD_test(x , y , alternative = "greater" , mu1=median(x) , mu2=median(y),
method_p="exact" , samplenum = 2000 , samplemethod = "R" , conf.level.sample = 0.95 )

}

\arguments{
  \item{x}{numeric vector of data values.

}
  \item{y}{numeric vector of data values.

}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".

}
  \item{mu1}{mean value of x, or median if not given

}
  \item{mu2}{mean value of y, or median if not given

}
  \item{method_p}{a character string specifying the method of calculating p-value, must be one of "exact" (default), "sampling".

}
  \item{samplenum}{The number of samples

}
  \item{samplemethod}{a discrete value indicating the method of sampling. "S" represents sample function sampling; "W" represents Cistern Sampling; "R" represents Put-back sampling.

}
  \item{conf.level.sample}{p-value confidence level for SRS sampling

}
}
\details{
The test is to test the difference of deviance of two samples, which has different sample means.

The arguments \code{samplenum} and \code{samplemethod} only work when \code{method_p="sampling"}.
}
\value{

  \item{method}{the test uesd}
  \item{stat }{RMD of the original data.}
  \item{pval}{p-value for the test.}
  \item{conf.int}{the confidence interval for p-value(only if method_p = "sampling")}
  \item{alternative}{a character string describing the alternative hypothesis}

}
\references{Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.

}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
## A simple example
x=c(16.55, 15.36, 15.94, 16.43, 16.01)
y=c(16.05, 15.98, 16.10, 15.88, 15.91)
RMD_test(x , y , alternative = "greater" )
}

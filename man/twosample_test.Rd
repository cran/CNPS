\name{twosample_test}
\alias{twosample_test}

\title{Comprehensive two-sample permutation tests

}
\description{Perform two-sample permutation test on vectors of data.

}
\usage{

twosample_test (x , y , alternative = "greater" , score = "wilcoxon" ,
method_p = "sampling" , samplenum = 2000 ,samplemethod="R",
conf.level.sample = 0.95 , conf.diff = TRUE, conf.level.diff = 0.95)
}

\arguments{
  \item{x}{numeric vector of data values.

}
  \item{y}{numeric vector of data values.

}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided", "greater"(default) or "less".

}
  \item{score}{a discrete value indicating the type of score. There are "original", "Wilcoxon", "van" and "exp" to be selected.

}
  \item{method_p}{a string indicating what method to use for p-value. "sampling" represents sampling; "asymptotic" represents using large sample approximations; "exact" represents Iterate through all combinations.

}
  \item{samplenum}{The number of samples

}
  \item{samplemethod}{a discrete value indicating the method of sampling. "S" represents sample function sampling; "W" represents Cistern Sampling; "R" represents Put-back sampling.

}
  \item{conf.level.sample}{p-value confidence level for SRS sampling

}
  \item{conf.diff}{a logical indicating whether to calculate the confidence interval of drift parameters.

}
  \item{conf.level.diff}{the level of confidence of drift parameters.

}
}
\details{
  \code{score} has 4 options: "original", "Wilcoxon", "van" and "exp". When choosing "original", the test is based on the original data; if \code{score = "Wilcoxon"}, the test is baesd on rank-sum; if \code{score = "van"}, the test is based on Van der Waerden score; if \code{score = "exp"}, the test is based on exponential score.

\code{samplenum} and \code{samplemethod} only work when \code{method_p="sampling"}. Similarly, \code{conf.level.diff} only works when \code{conf.diff =TRUE}.
}
\value{

  \item{method }{ the test used.}
  \item{score}{ the score which is used.}
  \item{stat }{   the statistic of the original data.}
  \item{conf.int }{the confidence interval for p-value(only if method_p = "sampling")}
  \item{pval }{p-value for the test}
  \item{alternative }{ a character string describing the alternative hypothesis.}
  \item{addition    }{a character string describing the Hodges-Lehmann estimate and the confidence interval of the drift parameter.}

}
\references{Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.

}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
## A simple example
x = c(1,2,3,4,5)
y = c(2,3,4,5,6)
twosample_test(x,y,samplemethod = "R" )

}

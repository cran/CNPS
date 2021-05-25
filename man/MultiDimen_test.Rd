\name{MultiDimen_test}
\alias{MultiDimen_test}

\title{Multivariate Permutation Test and Paired Comparisons

}
\description{Performs multivariate permutation tests, including paired tests.

}
\usage{

MultiDimen_test (data , stat = "HT",pair=FALSE, method_p = "sampling",rank = FALSE,
diff = FALSE , samplenum = 1000)

}

\arguments{
  \item{data}{
  a matrix or data frame of data values.

}
  \item{stat}{
  a character string specifying the statistic, must be one of "HT" (default), "tmax", "tmaxabs", "wsum", "zmax", "zmaxabs".

}
  \item{pair}{
   a logical indicating whether you want a paired test.

}
  \item{method_p}{
  a character string specifying the method of calculating p-value, must be one of "sampling" (default), " asymptotic", "exact".

}
  \item{rank}{
   a logical indicating whether you want Wilcoxon test.

}
  \item{diff}{
  a logical indicating whether you want to present which variables are different.

}
  \item{samplenum}{
      a number specifying the number of sampling.

}
}
\details{
The test can be used for multivariate permutation test and multivariate paired comparisons.

When doing multivariate paired comparisons, that is \code{pair = TRUE}, the statistic \code{wsum} is not suitable. Meanwhile, asymptotic method can only be used when statistic is \code{HT}. Besides, the second last column of the data must only contain two unique numbers to represent the two samples; the last column represents different pairs.

When doing multivariate permutation test, that is \code{pair = FALSE}, the statistic \code{zmax} and \code{zmaxabs} are not suitable. Meanwhile, the last column of the data must only contain 0 and 1 to represent the two samples. Besides, asymptotic method can not be used when statistic is \code{tmax} or \code{tmaxabs}.
}
\value{

  \item{method}{the test which is used.}
  \item{score }{ a character string describing the score used for test.}
  \item{stat}{ the test statistic.}
  \item{pval}{ p-value for the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{addition}{ a character string describing which variable is different in two samples.(presents only if pair = FALSE)}

}
\references{
Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.

}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}





\examples{
## Multivariate permutation test
data = matrix(c(6.81, 6.16, 5.92, 5.86, 5.80, 5.39,
              6.68, 6.30, 6.12, 5.71, 6.09, 5.28,
              6.34, 6.22, 5.90, 5.38, 5.20, 5.46,
              6.68, 5.24, 5.83, 5.49, 5.37, 5.43,
              6.79, 6.28, 6.23, 5.85, 5.56, 5.38,
              6.85, 6.51, 5.95, 6.06, 6.31, 5.39,
              6.64, 5.91, 5.59, 5.41, 5.24, 5.23,
              6.57, 5.89, 5.32, 5.41, 5.32, 5.30,
              6.84, 6.01, 5.34, 5.31, 5.38, 5.45,
              6.71, 5.60, 5.29, 5.37, 5.26, 5.41,
              6.58, 5.63, 5.38, 5.44, 5.17, 6.62,
              6.68, 6.04, 5.62, 5.31, 5.41, 5.44),
              nrow = 12,ncol = 6,byrow = TRUE
)
data=as.matrix(data)
index=c(rep(0,6),rep(1,6))
data = cbind(data,index)
x = MultiDimen_test(data ,  rank = FALSE ,  method_p = "sampling", samplenum = 100
, stat = "HT",diff = TRUE )
y = MultiDimen_test(data ,  rank = FALSE ,  method_p = "sampling", samplenum = 100
, stat = "tmax",diff = TRUE)
z = MultiDimen_test(data ,  rank = TRUE , method_p = "sampling"  , stat = "HT"
, samplenum = 100,diff = TRUE)

## Multivaraite paired comparisons
data = matrix(c(82, 60,  72, 62,
                75, 71,  70, 68,
                85, 59,  87, 64,
                90, 77,  87, 78),
              nrow = 4,ncol = 4,byrow = TRUE
)
x = data[,c(1,2)]
y = data[,c(3,4)]
data = cbind(rbind(x,y) , c(0,0,1,1) , c(1,2,1,2))
MultiDimen_test(data , method_p = "exact" , pair = TRUE)

}


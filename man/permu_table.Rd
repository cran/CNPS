\name{permu_table}
\alias{permu_table}

\title{Permutation Tests for Contingency Tables

}
\description{Performs permutation tests on contingency tables, including tables with ordered or disordered categories.

}
\usage{

permu_table(data , permu = "row" , row = NULL , col = NULL , fix = "row" ,
samplenum = 1000)
}

\arguments{
  \item{data}{ a matrix or data frame of data values.

}
  \item{permu}{a character string specifying the method of generating permutation samples, must be one of "row" (default), "col".

}
  \item{row}{ a numeric vector of the order of row categories. "row =NULL" indicates the categories are disordered.


}
  \item{col}{a numeric vector of the order of colnum categories."col =NULL" indicates the categories are disordered.


}
  \item{fix}{a character string specifying the group characteristic, must be one of "row" (default), "col". This argument is used for JT test when both "row" and "col" arguments are not "NULL".

}
  \item{samplenum}{ a number specifying the number of sampling.

}
}
\details{
The test can deal with the contingency tables with or without ordered categories. if both \code{row} and \code{col} are \code{null}, the data is viewed as normal contingency table with disordered categories. The test will use chi-square statistic.

If one of \code{row} and \code{col} is not \code{null}, the data is viewed as contingency table with one ordered characteristic. The vector you put in represents the order. For another characteristic, the Wilcoxon test can be performed if it has two classes; if the characteristic has several different values, the Kruskal-Wallis test can be used.

if both \code{row} and \code{col} are not \code{null}, that is both characteristics are ordered, the different values of one of the characteristic can be treated as observations and the JT test used for the other characteristic.
}

\value{

  \item{method}{a character string describing the type of test.}
  \item{stat}{ the value of the test statistic with a name describing it.}
  \item{pval}{the p-value for the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{conf.int}{ 95\% confidence interval of p-value(presents if either "row" and "col" arguments are not "NULL".)}

}
\references{Higgins, J. J. (2004). An introduction to modern nonparametric statistics. Pacific Grove, CA: Brooks/Cole.

}
\author{
Jiasheng Zhang, Feng Yu, Yangyang Zhang, Siwei Deng. Tutored by YuKun Liu and Dongdong Xiang.
}

\examples{
## generate a contingency table
x1=c(10,12, 17, 30)
x2=c( 9, 9, 11, 35)
x3=c( 7, 8, 12, 43)
data = rbind(x1,x2,x3)

## without ordered categories
permu_table(data)

## with ordered column categories
permu_table(data , col = c(1,2,3,4) )

## with ordered row categories
permu_table(data , row = c(1,2,3))

## with ordered row and column categories
permu_table(data , col = c(1,2,3,4),row = c(1,2,3),fix = "row")

}


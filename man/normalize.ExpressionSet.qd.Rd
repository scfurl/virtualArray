\name{normalize.ExpressionSet.qd}
\alias{normalize.ExpressionSet.qd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
quantile discretization normalization
}
\description{
Performs a discretization of continous expression levels to discrete steps.
}
\usage{
normalize.ExpressionSet.qd(ExpressionSet = NULL, nbin = 8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
An ExpressionSet object.
}
  \item{nbin}{
%%     ~~Describe \code{nbin} here~~
Numeric vector of length 1; indicates the number of binnings for the discretization.
}
  \item{\dots}{
Can be used to pass on parameters to underlying functions.
}
}
\details{
Please see the references for details.
}
\value{
An ExpressionSet object with adjusted expression matrix.
}
\references{
%% ~put references to the literature/web site here ~
Warnat P, Eils R, Brors B: Cross-platform analysis of cancer microarray data improves gene expression based classification of phenotypes. BMC bioinformatics 2005, 6:265.
}
\author{
Andreas Heider
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArrayExpressionSets
}
\examples{
# library(affydata)
# data(Dilution)
# Dilution <- rma(Dilution)
# normalize.ExpressionSet.qd(Dilution)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtualArray }
\keyword{ batch effect removal }
\keyword{ normalization }
\keyword{ cross-platform }

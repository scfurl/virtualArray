\name{normalize.ExpressionSet.mrs}
\alias{normalize.ExpressionSet.mrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
median rank scores normalization
}
\description{
The median rank scores algorithm is a modification of the quantile normalization for multi-platform microarray data.
}
\usage{
normalize.ExpressionSet.mrs(ExpressionSet = NULL, Batch = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionSet}{
An ExpressionSet object; if it contains a "Batch" column in its pData slot, this information can be used instead of the "Batch" parameter.
}
  \item{Batch}{
A character or numeric vector containing batch contribution of the samples. 
The order of the items should be the same as the samples (column) in the supplied ExpressionSet.
Defining this parameter overrides the information of the pData slot.
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
# normalize.ExpressionSet.mrs(Dilution,Batch=c("A","B","A","B"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtualArray }
\keyword{ batch effect removal }
\keyword{ normalization }
\keyword{ cross-platform }

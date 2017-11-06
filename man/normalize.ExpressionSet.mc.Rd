\name{normalize.ExpressionSet.mc}
\alias{normalize.ExpressionSet.mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mean centering normalization
}
\description{
Performs mean centering on the supplied ExpressionSet object.
}
\usage{
normalize.ExpressionSet.mc(ExpressionSet = NULL, Batch = NULL, ...)
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
Tibshirani R, Hastie T, Narasimhan B, Chu G: Diagnosis of multiple cancer types by shrunken centroids of gene expression. Proceedings of the National Academy of Sciences of the United States of America 2002, 99:6567-6572.
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
# normalize.ExpressionSet.mc(Dilution,Batch=c("A","B","A","B"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtualArray }
\keyword{ batch effect removal }
\keyword{ normalization }
\keyword{ cross-platform }

\name{normalize.ExpressionSet.nordi}
\alias{normalize.ExpressionSet.nordi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normal discretization normalization
}
\description{
Performs a discretization of continous expression levels to discrete steps.
}
\usage{
normalize.ExpressionSet.nordi(ExpressionSet = NULL, pvalue = 0.01, alpha = 0.05, ...)
}
\arguments{
  \item{ExpressionSet}{
An ExpressionSet object
}
  \item{pvalue}{
The pvalue determines the sensitivity for detecting outliers in each column of the gene expression matrix.
}
  \item{alpha}{
The alpha value determines the size of the tails of the normal distributions considered to be over or under expressed.
}
  \item{\dots}{
Can be used to pass on parameters to underlying functions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Please see the references for details.
}
\value{
An ExpressionSet object with adjusted expression matrix.
}
\references{
"GenMiner: Mining informative association rules from genomic data." 
Ricardo Martinez, Claude Pasquier and Nicolas Pasquier, Proceedings of the IEEE BIBM international conference on Bioinformatics and Biomedecine, pages 15-22, IEEE Computer Society, 2007.

Martinez R, Pasquier N, Pasquier C: GenMiner: mining non-redundant association rules from integrated gene expression data and annotations. Bioinformatics 2008, 24:2643-2644.}
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
# normalize.ExpressionSet.nordi(Dilution)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtualArray }
\keyword{ batch effect removal }
\keyword{ normalization }
\keyword{ cross-platform }

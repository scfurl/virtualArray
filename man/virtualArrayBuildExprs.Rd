\name{virtualArrayBuildExprs}
\alias{virtualArrayBuildExprs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Format ExpressionSets as data.frames with identifiers in first column
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function takes ExpressionSets as input and formats the expression matrix into a data.frame whose identifiers are put into the first column
}
\usage{
virtualArrayBuildExprs(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
An ExpressionSet
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is normally only called by "virtualArray.ExpressionSet".
It can be used, however, to extract the "exprs" slot of an ExpressionSet as a data.frame.
At the same time the identifiers (rownames) are saved in column "1" of the data.frame.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The value returned is a data.frame which is based on the expression matrix of the input ExpressionSet.
}
% \references{%% ~put references to the literature/web site here ~}
\author{
%%  ~~who you are~~
Andreas Heider
}
% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
## first we need to load dummy data
library(affydata)
data(Dilution)
## we apply RMA to get an ExpressionSet
Dilution <- rma(Dilution,normalize=FALSE)
## now we derive a data.frame from the expression values
new_exprs_data.frame <- virtualArrayBuildExprs(Dilution)
head(new_exprs_data.frame,n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ virtual }
\keyword{ virtualArray }% __ONLY ONE__ keyword per line

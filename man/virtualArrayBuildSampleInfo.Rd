\name{virtualArrayBuildSampleInfo}
\alias{virtualArrayBuildSampleInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Build data.frame suitable as "pData" for several ExpressionSets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A data.frame is created that holds the sample names and names of the ExpressionSets of the supplied list of ExpressionSets.
}
\usage{
virtualArrayBuildSampleInfo(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A list of ExpressionSets. Each entry must be represented by a character vector.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is normally only called by "virtualArray.ExpressionSet".
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.frame is returned. It consists of the columns "Array.name", "Sample.name", "Batch" and "Covariate.1".
The "Batch" column represents the name of the originating ExpressionSet. 
The column "Covariate.1" is only a dummy holding running numbers, one for each sample.
It has to be edited and filled with meaningful information when "virtualArray" is run in supervised mode.
Finally the "sample_info" data.frame can be used to fill the "pData" slot for a merged ExpressionSet derived from several ones.
}

% \references{%% ~put references to the literature/web site here ~}

\author{
%%  ~~who you are~~
Andreas Heider (2011)
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
data(sample.ExpressionSet)
## we apply RMA to get an ExpressionSet
Dilution <- rma(Dilution,normalize=FALSE)
## we store the names of the ExpressionSets in a list
all_ExpressionSets <- list("Dilution","sample.ExpressionSet")
## now we generate a new data.frame suitable for "pData"
virtualArrayBuildSampleInfo(all_ExpressionSets)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ virtualArray }
\keyword{ virtual }% __ONLY ONE__ keyword per line

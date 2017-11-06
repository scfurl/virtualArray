\name{GPLs}
\alias{GPLs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Provides precompiled mappings of NCBI GEO platform IDs (GPLs) to Bioconductor annotation packages.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This data object can be used to translate GPL accession IDs like "GPL570" to a Bioconductor annotation package like "hgu133plus2".
}
\usage{
GPLs
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
The "virtualArray" package provides tools to merge raw data from multiple microarray platforms.
The source of this data might be an entry from the NCBI GEO web service using e.g. "getGEO" from the "GEOquery" package.
The annotation slot of the retreived ExpressionSets is filled with a GPL entry like "GPL570" which corresponds to the Bioconductor annotation package "hgu133plus2".
The "GPLs" list can be used to map between these IDs and replace the annotation slot. A reverse map can be achieved using "revmap" from the annotate package.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of Bioconductor annotation packages (redundant) named by their corresponding NCBI GEO GPL.
}

% \references{%% ~put references to the literature/web site here ~}

\author{
%%  ~~who you are~~
Andreas Heider (2013)
}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
GPLs["GPL570"]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ GPL }
\keyword{ GEO }
\keyword{ mapping }
\keyword{ virtual }
\keyword{ virtualArray }% __ONLY ONE__ keyword per line

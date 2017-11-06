\name{virtualArrayLoadRaw}
\alias{virtualArrayLoadRaw}
\title{
Load user provided raw data files in the virtualArray directory structure
}
\description{
This function loads all raw data files provided by the user in the virtualArray directory structure, 
that was created by virtualArrayDirs in the first place.
}
\usage{
virtualArrayLoadRaw(root_dir = getwd(), affy_sum = "rma")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root_dir}{
The top level directory holding the rawdata and the "sample_info.txt" file.
}
  \item{affy_sum}{
Select the summerization method to use for Affymetrix chips; defaults to "RMA".
}
}
\details{
The function first reads in the virtualArray directory including raw data files. 
Based on this a "sample_info" data.frame is created holding information about chip type and the corresponding file name. 
Then "sample_info" is checked for the presence of certain microarray platforms. 
If present, the function loads raw data of all platforms. The "sample_info" data.frame is exported to "sample_info.txt". 
It is up to the user to change the values of the "Covariate 1" column to something meaningful, thus providing the additional information for the batch effect removal in the final step of the process.
}
\value{
This function does not return a value, but creates objects in the current environment holding raw data of every provided chip type. The objects are named "rawdata_*". Where "*" is the chip type.
}

% \references{%% ~put references to the literature/web site here ~}

\author{
Andreas heider (2011)
}

%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
## EMPTY
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtual }
\keyword{ combine }% __ONLY ONE__ keyword per line

\name{virtualArrayDirs}
\alias{virtualArrayDirs}
\title{
%%  I don'T know.
Build virtualArray directory structure to hold raw data of multiple microarray chips
}
\description{
This function simply creates a directory tree in the current working directory or a user specified one. 
This directory tree helps organizing raw data files and at the same time gives R/Bioconductor the clue who to handle the provided data.
}
\usage{
virtualArrayDirs(root_dir = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root_dir}{
If declared, this variable sets the top level directory where the whole tree will be placed.
}
}
\details{
The directory tree consists of a top level folder called "rawdata". 
Inside this folder there are directories for "Affymetrix", "Agilent" and "Illumina" respectively. 
Each folder holds sub-folders for the specific chips of the respective manufacturer. 
The user has to supply copies of the raw data into the appropriate directories for the package to work as expected.
}
\value{
This function returns no value. It just creates a directory tree.
}
% \references{%% }

\author{
Andreas Heider (2011)
}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
## we use the function in the current working directory
## please note that this will write to your hard drive
virtualArrayDirs()
## now we can see what happened
dir(recursive=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtual }
\keyword{ combine }% __ONLY ONE__ keyword per line
\keyword{ merge }
\keyword{ array }

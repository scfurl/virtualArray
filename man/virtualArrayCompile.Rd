\name{virtualArrayCompile}
\alias{virtualArrayCompile}
\title{
%%  I don'T know.
Compiles a single virtual array holding user provided microarray data from multiple platforms.
}
\description{
The function takes ExpressionSets of rawdata provided by the "virtualArrayLoadRaw" function together with annotation data.frames provided by the function "virtualArray.load.annot" and matches the rawdata on the basis of the given identifier. Currently only "Gene.Symbol" is supported by "virtualArrayCompile" but this will be extendend in future versions. Meanwhile you can use "virtualArray.ExpressionSet" to use other identifiers.
}
\usage{
virtualArrayCompile(root_dir = getwd(), identifier = "Gene.Symbol", covars= 3, virtualArray = virtualArray)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{root_dir}{
The top level directory holding the rawdata and the "sample_info.txt" file.
}
  \item{identifier}{
The annotation identifier based on which the comparison between the different platforms is made.
}
  \item{covars}{
Numer of columns to read from "sample_info" data.frame or ".txt". 
Defaults to 3, which denotes only to use the sample/array naming and the "Batch" column.
A 4 expects the "Covariate.1" column to be filled with meaningful information regarding grouping of samples.
When the value is the character vector "all" even more columns can be used.
}
  \item{virtualArray}{
A list which holds the intermediate output of the "virtArray" package.
It contains both raw data as well as annotation data.frames.
It defaults to "virtualArray" but can be any other name.
}
}

% \details{}

\value{
An ExpressionSet is created, holding data of all supplied microarray platforms.
}
% \references{%% ~put references to the literature/web site here ~}

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
# This is example is commented out, because it would slow down building the package!
# virtualArrayDirs()
# getGEOSuppFiles(baseDir="rawdata/Affymetrix/U133Plus2/",makeDirectory=F,GEO="GSM589512")
# gunzip(filename="rawdata/Affymetrix/U133Plus2/GSM589512.CEL.gz",remove=T)
# getGEOSuppFiles(baseDir="rawdata/Affymetrix/U133A/",makeDirectory=F,GEO="GSM589506")
# getGEOSuppFiles(baseDir="rawdata/Affymetrix/U133A/",makeDirectory=F,GEO="GSM589509")
# gunzip(filename="rawdata/Affymetrix/U133A/GSM589506.CEL.gz",remove=T)
# virtualArray <- virtualArrayLoadRaw()
# my_ExpressionSet <- virtualArrayCompile(virtualArray=virtualArray)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtual }
\keyword{ compare }% __ONLY ONE__ keyword per line

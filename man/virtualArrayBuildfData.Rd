\name{virtualArrayBuildfData}
\alias{virtualArrayBuildfData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapses expression values according to given identifiers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Expression values are collapsed according to given identifiers, e.g. gene symbols.
In the same run, a data.frame with 2 columns to fill the "fData" slot of the ExpressionSet is build.
The data is pulled from the Bioconductor annotation package defined in the "annotation" slot of the ExpressionSet.
}
\usage{
virtualArrayBuildfData(x, identifier = "SYMBOL", collapse_fun = median)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Name of an ExpressionSet as a character vector
}
  \item{identifier}{
%%     ~~Describe \code{identifier} here~~
Wich identifier to pull from the annotation (default="SYMBOL").
Theoretically it can be anything that is supported by the annotation package.
Practically it is limited to identifiers giving a 1 to 1 mapping, such as ENTREZID, SYMBOL, GENENAME, UNIPROT, UNIGENE, etc.
Identifiers giving 1 to 1+ mappings will be supported in future versions of the package.
}
  \item{collapse_fun}{
%%     ~~Describe \code{collapse_fun} here~~
Which function to use to treat multiple lines targetting the same gene (default=median)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is normally only called from "virtualArray.ExpressionSet".
You can use it, however, to collapse expression values referring to the same identifier together with their corresponding 
annotation in the "exprs" and "fData" slots by means of a user specified function (default is "median").
Note, that it is critical to define the correct Bioconductor annotation package in the "annotation" slot of the ExpressionSet
and that the name of the ExpressionSet must be supplied as a character vector, due to the implementation in "virtualArray.ExpressionSet".
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The value returned is an ExpressionSet whose "fData" slot has been filled with the selected identifiers.
The rows in the expression matrix ("exprs" slot) have been collapsed to the selected identifiers.
Therefore the size of the expression matrix has decreased during the process.
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
## now we collapse the expression values as stated in "Details"
Dilution_genesymbols <- virtualArrayBuildfData(x="Dilution")
## as you can see, we now have only one row per gene symbol
Dilution_genesymbols
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ ExpressionSet }
\keyword{ virtualArray }
\keyword{ virtual }% __ONLY ONE__ keyword per line

\name{virtualArrayComBat}
\alias{virtualArrayComBat}
\alias{virtualArrayComBat,ExpressionSet-method}
\alias{virtualArrayComBat,data.frame-method}
\alias{virtualArrayComBat,character-method}
\title{
Removes batch effects from microarray derived expression matrices. Modified version.
}
\description{
This is a modified version of the R script "ComBat.R" (see references). It is used to adjust for batch effects in microarray data. The modification is restricted to make the script accept expression matrices and data.frames instead of plain text files.
}
\usage{
virtualArrayComBat(expression_xls, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_xls}{
%%     The expression matrix to adjust.
The expression matrix to adjust.
}
  \item{sample_info_file}{
%%     The sample information data.frame regarding batch contribution and possibly covariates.
The sample information data.frame regarding batch contribution and possibly covariates.
}
  \item{type}{
The type of input; Defaults to "txt".
}
  \item{write}{
Write output to external file or provide new expression matrix.
}
  \item{covariates}{
Describe which Covariates to use in the process and which to dismiss. The default is to use only "Batch".
}
  \item{par.prior}{
Logical; set prior parameters or not; Use prespecified values for the variables ("TRUE") or start a priori ("FALSE").
}
  \item{filter}{
Filter for genes not present in a given percentage of the samples. Requires present/absent calls in the data. Can be either "FALSE" or a numeric between "0" and "1". Recommended is "0.8" or "FALSE".
}
  \item{skip}{
Columns to skip in the input "expression_xls" matrix.
}
  \item{prior.plots}{
Create quantile-quantile and kernel density plots including prior estimates to assess the quality of the estimation.
}
}

% \details{%  ~~ If necessary, more details than the description above ~~}

\value{
%%  Returns a matrix holding adjusted expression values.
Returns a matrix holding adjusted expression values.
}
\references{
%% 
Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.
}
\author{
%%  
Original author: Johnson, WE, Rabinovic, A, and Li, C (2007)

Modified by: Andreas Heider (2011)
}
\note{
%%  ~~further notes~~
Original code by Johnson, WE, Rabinovic, A, and Li, C, made available in this package by Andreas Heider
}

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
\keyword{ batch effects }
\keyword{ batch }% __ONLY ONE__ keyword per line


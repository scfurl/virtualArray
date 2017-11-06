\name{virtualArrayExpressionSets}
\alias{virtualArrayExpressionSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Combine different ExpressionSets into one
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function selects all ExpressionSets in the current environment and builds a new single ExpressionSet of the raw data included in the input.
This is done by annotating the expression values with the selected identifiers, that are pulled from Bioconductor annotation packages.
Then lines targetting the same gene are collapsed by the specified function. In the next step compatible rows of the expression matrices are merged. 
As a final step batch effects resulting from different platforms or labs can be removed in a supervised or non-supervised mode.
}
\usage{
virtualArrayExpressionSets(all_expression_sets=FALSE, identifier = "SYMBOL", covars = "Batch", collapse_fun = median, removeBatcheffect = "EB", sampleinfo = FALSE, parallel="BiocParallel", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_expression_sets}{
Logical or a character vector. If "FALSE", "virtualArray" tries to catch all ExpressionSets in the current environment.
If set to a character vector holding names of ExpressionSets, these are used instead of all available ones.
}
  \item{identifier}{
%%     ~~Describe \code{identifier} here~~
annotation identifier by which the expression values are combined
}
  \item{covars}{
%%     ~~Describe \code{covars} here~~
Character vector of length "1" or longer. Used to define non-supervised (length = "1") or supervised (length > "1") mode.
THe default is "Batch". See details for more info.
}
  \item{collapse_fun}{
%%     ~~Describe \code{collapse_fun} here~~
The function to be used to collapse expression values targetting the same gene/identifier.
Defaults to "median".
}
  \item{removeBatcheffect}{
%%     ~~Describe \code{removeBatcheffect} here~~
Logical or character vector. "FALSE" will lead to just a combined ExpressionSet, you will then have to use other functions to remove the batch effects.
You can set it to "EB", "GQ", "MRS", QD" ,"NORDI" or "MC" to use empirical Bayes methods, gene quantiles, median rank scores, quantile discretization ,normal discretization or mean centering to remove batch effects, respectively.
}
  \item{sampleinfo}{
%%     ~~Describe \code{sampleinfo} here~~
This parameter selects in which way the information of the relationships between batches and samples/datasets will be supplied.
The default (FALSE) uses a sample_info data.frame that is generated on the fly from the pData slots of the supplied ExpressionSets with an additional "Batch" column.
If you run in non-interactive mode, you can specify a data.frame to be used as the input "sample_info".
Another option is to hand over a file name so the preconfigured text file can be fed into the procedure.
Note that a text file "sample_info.txt" can be created on the fly, so you can set it up manually in this case. You can select this option using "create".
}
  \item{parallel}{
  A character string or a logical to select which package to use for parallel processing.
  Defaults to "BiocParallel", but can be "multicore", "none" or FALSE as well.
}
  \item{\dots}{
Can be used to pass on parameters to underlying functions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The "covars" argument determines the mode of batch removal. 
It refers to the columns in the sample_info data.frame which contains information about all ExpressionSets, their samples and relations thereof.
The default value "Batch" will use only the different ExpressionSets for batch effect removal, this is referred to as the non-supervised mode.
The supervised mode is to be accessed by using a character vector with a length > 1 e.g. c("Batch","celltype"). 
In this case a column "celltype" must be common to all pData slots of all datasets prior to invoking the package.
The default name of the batch column ("Batch") can even be replaced to match another column.
The sample_info data.frame is generated on the fly from the pData slots of the supplied ExpressionSets, during this procedure a "Batch" column is generated.
All columns are preserved and common ones joined. These can be used as additional covariates during batch effect removal.
In this case the sample_info data.frame has to be modified manually to contain more information on the batches in additional columns.
Please note, that during computation you will be notified that "sample_info.txt" has been written to your current working directory for you to modify and save it.
If you do so, please select "y" to use the additional columns. Also note that you can not provide a covariate that is distributed only in one batch, this way the procedure will fail.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A new ExpressionSet is returned that combines all ExpressionSets from the current environment.
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
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile, normalize.ExpressionSet.nordi, normalize.ExpressionSet.mrs, normalize.ExpressionSet.qd, normalize.ExpressionSet.gq
}
\examples{
# Due to the flexibility of this function and the time 
# it takes to get meaningful results, please see the 
# vignette for a comprehensive example, governing 
# several modes of usage. Thanks.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ virtual }
\keyword{ combine }
\keyword{ virtualArray }% __ONLY ONE__ keyword per line

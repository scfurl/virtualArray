\name{virtualArrayMergeRecurse}
\alias{virtualArrayMergeRecurse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merge a list of data.frames by common rows
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A list of data.frames is merged together by rows common to all data.frames.
The column by which the merging is performed can be selected with "by".
}
\usage{
virtualArrayMergeRecurse(dfs, by, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfs}{
%%     ~~Describe \code{dfs} here~~
A list of data.frames
}
  \item{by}{
%%     ~~Describe \code{dfs} here~~
Select the columns to be used for the comparison
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Parameters passed on to the "merge" function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function was adopted from the "reshape" package and modified to be suitable to handle problems arising when combining data from different microarray platforms in the "virtualArray" package.
This function is normally only called by "virtualArray.ExpressionSet".
It can be used, however, to merge several data.frames by common rows, whereas rows that do not match between the data.frames are discarded.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.frame consisting of the matching parts of the input data.frames.
}

% \references{%% ~put references to the literature/web site here ~}

\author{
%%  ~~who you are~~
Andreas heider (2011)
}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
## we set up 2 example data.frames
authors <- data.frame(
name = I(c("Tukey", "Venables", "Tierney", "Ripley", "McNeil")),
nationality = c("US", "Australia", "US", "UK", "Australia"),
deceased = c("yes", rep("no", 4)))
books <- data.frame(
name = I(c("Tukey", "Venables", "Tierney",
"Ripley", "Ripley", "McNeil", "R Core")),
title = c("Exploratory Data Analysis",
"Modern Applied Statistics ...",
"LISP-STAT",
"Spatial Statistics", "Stochastic Simulation",
"Interactive Data Analysis",
"An Introduction to R"),
other.author = c(NA, "Ripley", NA, NA, NA, NA,
	  "Venables & Smith"))
## let's have a look at the data
authors
books
## we store the data.frames in a list
my_dfs <- list(authors,books)
## now we can combine the 2 data.frames
virtualArrayMergeRecurse(my_dfs,by="name")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ virtual }
\keyword{ virtualArray }% __ONLY ONE__ keyword per line

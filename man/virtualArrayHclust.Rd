\name{virtualArrayHclust}
\alias{virtualArrayHclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot a hclust object with colored labels
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function takes a hclust object and a vector of colors as its input and plots the hclust object with colored labels.
}
\usage{
virtualArrayHclust(hclust, lab = hclust$labels, lab.col = rep(1, length(hclust$labels)), hang = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hclust}{
%%     ~~Describe \code{hclust} here~~
a hclust object
}
  \item{lab}{
%%     ~~Describe \code{lab} here~~
labels to use
}
  \item{lab.col}{
%%     ~~Describe \code{lab.col} here~~
colors to use for each label
}
  \item{hang}{
%%     ~~Describe \code{hang} here~~
see "hclust"
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
see "hclust"
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
A modified hclust object is returned.
}

% \references{%% ~put references to the literature/web site here ~}

\author{
%%  ~~who you are~~
Original author: Eva Chan

Modified by: Andreas Heider (2011)
}
\note{
%%  ~~further notes~~
The function was written by Eva Chan and incorporated into this package by Andreas Heider
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
virtualArray-package, virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
# we use the USArrests dataset 
new.arrests <- USArrests
# now we need to specify a function to create colors
color_me <- function(x){
					if(x >= 100){"red"}
					else{"blue"}
					}
# lets build the colors in relation to the "Assaults" statistic
new.arrests[,5] <- sapply(new.arrests[,2],color_me)
# we need a distence matrix to build our tree
dist1 <- dist(new.arrests[,1:4],method="euclidian")
# and we need a conventional hclust object
hc_nocolor <- hclust(dist1,method="average")
# now we can hook up to that object and change the label colors
# note that this call will automatically plot the tree
virtualArrayHclust(hc_nocolor,lab=rownames(new.arrests),lab.col=new.arrests[,5])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ combine }
\keyword{ virtual }
\keyword{ virtualArray }% __ONLY ONE__ keyword per line

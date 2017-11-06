\name{virtualArray-package}
\alias{virtualArray-package}
\alias{virtualArray}
\docType{package}
\title{
Combine raw data of several ExpressionSets into a "virtual array".
}
\description{
Using this package one can create a new ExpressionSet object, that holds both 
phenoData and expression values of different microarray platforms.
}
\details{

\tabular{ll}{
Package: \tab virtualArray\cr
Type: \tab Package\cr
Version: \tab 1.3.1\cr
Date: \tab 2012-12-19\cr
License: \tab GPL V3\cr
LazyLoad: \tab yes\cr
}
Current tools for the analysis of microarray data only allow the comparison of datasets generated on the same platform and chip generation, 
and restrict meta-analysis studies to the evaluation of study results from different groups, 
rather than the direct comparison of the available raw datasets.
With the help of the virtualArray package it becomes easy to combine raw expression data of different microarray platforms into one "virtual array".
Thus, the user may compare his own data to other datasets including public sources, 
regardless of the platform and chip generation used, or perform meta-analysis directly based on available raw datasets.
The package generates a combined virtual array as a "ExpressionSet" object from different datasets by matching raw data entries based on probe,
transcript, gene or protein identifiers. 
Redundancies, gaps, and batch effects are removed before proceeding with data analysis.

virtualArray consists of several subsequent functions, requiring minimal user input. 
Briefly, (1) raw data are loaded into R, (2) probe sets for each platform are annotated, (3) genes, 
proteins or transcripts common to all platforms are matched, 
including checks for redundancy or missing values, data are (4) compiled into a new "virtual array", 
(5) normalized and (6) subjected to batch effect removal using empirical Bayes methods [1], or other batch effect removal methods. Six of which have been implemented into this package. 
The generated "virtual array" can than be directly analyzed in R/Bioconducor or exported for use in other suitable software (e.g. MeV).

There are essentially four modes of operation: 

Firstly, the "virtualArrayCompile" function can integrate the major (but not all) human microarray platforms in a default mode requiring 
only minimal user input.

The second way is built around the "virtualArray.ExpressionSet" function. 
This approach allows to integrate any kind of raw expression data that can be loaded into an ExpressionSet object in R/BioC. 
The downside of this mode is that the user will have to deal with details such as log2-transformations, 
16 bit - 20 bit transformations, assignment of correct annotations, etc.

Finally, each of these two approaches can be used in a supervised or non-supervised mode. 
The non-supervised mode uses empirical Bayes methods (implemented through "ComBat.R", [1]) to adjust for batch effects between the single datasets. 
In the supervised mode the user assigns additional covariates next to the batch assignment. 
This means e.g. "treated" and "untreated", or "non-differentiated" and "differentiated". 
Note that this information has to be valid, as it impacts the results you will get.

Last but not least it is possible to use the package to integrate data without batch effect removal.
This way you can use other methods of batch effect removal than the empirical Bayes methods approach.

The combined data is presented as a regular Bioconductor "ExpressionSet" object, which permits using all of R/Bioconductor's 
power on the dataset as a whole.

To load the package type:

> library(virtualArray)
}
\author{
Andreas Heider

Maintainer: Andreas Heider \email{aheider@trm.uni-leipzig.de}
}
% \references{}
\keyword{ package }
\seealso{
virtualArray.ExpressionSet, virtualArrayCompile
}
\examples{
## Please see the vignette for a comprehensive example
}

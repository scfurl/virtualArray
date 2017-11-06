virtualArrayBuildExprs <-
function(x) {
							#dtfrm <- data.frame(exprs(eval(as.symbol(x))))
							dtfrm <- data.frame(exprs(x),stringsAsFactors=FALSE,check.names=FALSE)
							dtfrm <- cbind(identifier=as.character(rownames(dtfrm)),dtfrm)
							dtfrm[,1] <- as.character(dtfrm[,1])
							#dtfrm[,1] <- as.factor(dtfrm[,1])
							#print(head(dtfrm))
							#print(dim(dtfrm))
							#x <- data.frame(exprs(eval(as.symbol(x))))
							#cat("\nFirst 5 rows of dtfrm\n")
							#print(head(dtfrm))
							#cat("Size of dtfrm")
							#print(dim(dtfrm))
							return(dtfrm)
							}


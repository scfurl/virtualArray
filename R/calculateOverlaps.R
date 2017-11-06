calculateOverlaps <- function(listOfIdentifiers){
mtrx1 <- mapply(
	Y=names(listOfIdentifiers),
	FUN=function(Y){
		mapply(	x=names(listOfIdentifiers),
				FUN=function(x){
					length(
						intersect(listOfIdentifiers[[x]],listOfIdentifiers[[Y]])
						)
						}
				)
					}
		)
colnames(mtrx1) <- apply(X=as.array(colnames(mtrx1)),FUN=function(X){paste("vs.",X,sep="")},1)
mtrx2 <- cbind(full=as.numeric(lapply(listOfIdentifiers,length)),mtrx1)
return(mtrx2)
	}

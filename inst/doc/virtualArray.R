### R code from vignette source 'virtualArray.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: virtualArray.Rnw:68-69
###################################################
library(virtualArray)


###################################################
### code chunk number 3: virtualArray.Rnw:133-136
###################################################
library("GEOquery")
GSE23402 <- getGEO("GSE23402",GSEMatrix=T,AnnotGPL=FALSE)
GSE26428 <- getGEO("GSE26428",GSEMatrix=T,AnnotGPL=FALSE)


###################################################
### code chunk number 4: virtualArray.Rnw:141-143
###################################################
GSE23402 <- GSE23402[[1]][,1:24]
GSE26428 <- GSE26428[[1]]


###################################################
### code chunk number 5: virtualArray.Rnw:146-150 (eval = FALSE)
###################################################
## library(virtualArray)
## library(affy)
## GSE23402 <- GSE23402
## GSE26428 <- GSE26428


###################################################
### code chunk number 6: virtualArray.Rnw:156-158
###################################################
summary(exprs(GSE23402)[,1:3])
summary(exprs(GSE26428))


###################################################
### code chunk number 7: virtualArray.Rnw:168-170
###################################################
exprs(GSE23402) <- log2(exprs(GSE23402))
exprs(GSE26428) <- (exprs(GSE26428)/20*16)


###################################################
### code chunk number 8: virtualArray.Rnw:175-177
###################################################
summary(exprs(GSE23402)[,1:4])
summary(exprs(GSE26428))


###################################################
### code chunk number 9: virtualArray.Rnw:183-185
###################################################
annotation(GSE23402) <- "hgu133plus2"
annotation(GSE26428) <- "hgug4112a"


###################################################
### code chunk number 10: virtualArray.Rnw:190-191
###################################################
my_virtualArrays <- NULL


###################################################
### code chunk number 11: options
###################################################
options(width=60)
if(require(BiocParallel))
    register(MulticoreParam(verbose=TRUE))


###################################################
### code chunk number 12: virtualArray.Rnw:209-210
###################################################
my_virtualArrays$iPSC_hESC_noBatchEffect <- virtualArrayExpressionSets()


###################################################
### code chunk number 13: virtualArray.Rnw:215-216
###################################################
my_virtualArrays$iPSC_hESC_withBatchEffect <- virtualArrayExpressionSets(removeBatcheffect=FALSE)


###################################################
### code chunk number 14: virtualArray.Rnw:221-225
###################################################
pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[5] <- 
	c(as.character(pData(GSE23402)[,8]),as.character(pData(GSE26428)[,1]))
pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[6] <- 
	c(rep("red",24),rep("blue1",3))


###################################################
### code chunk number 15: virtualArray.Rnw:243-246
###################################################
dist_iPSC_hESC_noBatchEffect <- 
	dist(t(exprs(my_virtualArrays$iPSC_hESC_noBatchEffect)), 
	method="euclidian")


###################################################
### code chunk number 16: virtualArray.Rnw:258-261
###################################################
hc_iPSC_hESC_noBatchEffect <- 
	hclust(dist_iPSC_hESC_noBatchEffect, method="average")
hc_iPSC_hESC_noBatchEffect$call <- NULL


###################################################
### code chunk number 17: virtualArray.Rnw:284-289
###################################################
virtualArrayHclust(hc_iPSC_hESC_noBatchEffect,
	lab.col=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,6],
	lab=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,5],
	main="batch effect removed",cex=0.7,
	xlab="sample names")


###################################################
### code chunk number 18: virtualArray.Rnw:346-347 (eval = FALSE)
###################################################
## my_virtualArrays$iPSC_hESC_supervised <- virtualArrayExpressionSets(sampleinfo="create")


###################################################
### code chunk number 19: virtualArray.Rnw:350-351
###################################################
my_virtualArrays$iPSC_hESC_supervised <- virtualArrayExpressionSets(sampleinfo=sample_info_imported)


###################################################
### code chunk number 20: virtualArray.Rnw:354-357
###################################################
dist_iPSC_hESC_supervised <- 
	dist(t(exprs(my_virtualArrays$iPSC_hESC_supervised)), 
	method="euclidian")


###################################################
### code chunk number 21: virtualArray.Rnw:359-362
###################################################
hc_iPSC_hESC_supervised <<- 
	hclust(dist_iPSC_hESC_supervised, method="average")
hc_iPSC_hESC_supervised$call <- NULL


###################################################
### code chunk number 22: virtualArray.Rnw:365-370
###################################################
virtualArrayHclust(hc_iPSC_hESC_supervised,
	lab.col=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,6],
	lab=pData(my_virtualArrays$iPSC_hESC_noBatchEffect)[,5],
	main="batch effect removed - supervised mode",cex=0.7,
	xlab="sample names")


###################################################
### code chunk number 23: virtualArray.Rnw:375-378
###################################################
pca_supervised <- prcomp(t(exprs(my_virtualArrays$iPSC_hESC_supervised)))
plot(pca_supervised$x, pch=19, cex=2, col=c(rep("red",24),rep("blue",3),pch=17))
legend("topleft",c("GSE23402","GSE26428"),col=c("red","blue"),pch=19,cex=1)



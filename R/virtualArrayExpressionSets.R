virtualArrayExpressionSets <- 
function(all_expression_sets=FALSE,identifier="SYMBOL",covars="Batch",collapse_fun=median,removeBatcheffect="EB",sampleinfo=FALSE, parallel="BiocParallel", ...){
#require(affy)
if(class(try(utils::win.version(),silent=T))=="try-error" && parallel=="multicore")
	{require(multicore)
	lapply=mclapply}
loadedPackages <- loadedNamespaces()
names(loadedPackages) <- loadedNamespaces()
if(parallel=="BiocParallel"){
if(try(installed.packages()["BiocParallel",1],silent=T)=="BiocParallel"){require(BiocParallel)}
	else
		{
		source("http://www.bioconductor.org/biocLite.R")
		biocContribUrl <- sapply(biocinstallRepos(), contrib.url)
		try(install.packages(as.character("BiocParallel"),method="internal",dependencies=c("Depends","Imports"),repos= biocinstallRepos(),contriburl=biocContribUrl),silent=T)
		require(BiocParallel)
		}}
if(is.na(loadedPackages["BiocParallel"])!=TRUE){register(MulticoreParam(verbose=T));lapply <- bplapply}
# create list with loaded ExpressionSets
 if(class(all_expression_sets) != "character")
	{
	all_expression_sets <- sapply(X=ls(envir=.GlobalEnv),FUN=function(x){class(eval(as.symbol(x)))})
	all_expression_sets <- grep(all_expression_sets,pattern="ExpressionSet",value=TRUE)
	all_expression_sets <- names(grep(all_expression_sets,pattern="ExpressionSet",value=TRUE))
	}
all_expression_sets <- sapply(X=all_expression_sets,FUN=function(x){as.symbol(x)},USE.NAMES=TRUE)
# build new ExpressionSets with valid fData and remap features to selected identiifier in both exprs and fData, return list of ExpressionSets
expsts <- lapply(all_expression_sets,virtualArrayBuildfData,collapse_fun=collapse_fun,identifier=identifier)
sample_info <- virtualArrayBuildSampleInfo(all_expression_sets)
rm(all_expression_sets)
# build data.frames of the exprs and add features as column 1
dtfrms <- lapply(expsts,virtualArrayBuildExprs)
rm(expsts)
names_dtfrms <- names(dtfrms)
rownames_dtfrms <- lapply(dtfrms,FUN=function(X){rownames(X)})
names(rownames_dtfrms) <- names_dtfrms
rm(names_dtfrms)
message("Matching and merging all data sets. This could take some time...",appendLF=TRUE)
#merged <- na.omit(merge_recurse(dtfrms,by="identifier"))
merged <- virtualArrayMergeRecurse(dtfrms,by="identifier",incomparables=NA)
rm(dtfrms)
rownames(merged) <- merged[,1]
merged <- merged[2:dim(merged)[2]]
message("Size of expression matrix of whole dataset: ",dim(merged)[1]," rows and ",dim(merged)[2]," columns.",appendLF=TRUE)
rownames_dtfrms[["result"]] <- rownames(merged)
merged <- new("ExpressionSet",exprs=as.matrix(merged))
#require(affyPLM)
merged_norm <- normalize.ExpressionSet.quantiles(merged)
sample_info <- as.data.frame(sample_info)
#print(sample_info[,1])
#print(sampleNames(merged_norm))
sampleNames(merged_norm) <- sample_info[,1]
#sample_info[1] <- sampleNames(merged_norm)
#sample_info[2] <- sampleNames(merged_norm)

if(sampleinfo=="create" && length(covars)==1 ){
	write.table(sample_info,file="sample_info.txt",sep="\t",row.names=FALSE)
	message("The file 'sample_info.txt' has been written to your current working directory. Please modify it appropriately!",appendLF=TRUE)
	answer <- readline(prompt="Did you modify sample_info.txt? [y] or [n] ")
	sample_info_old <- sample_info
	sample_info <- read.table(file="sample_info.txt",sep="\t",header=TRUE)
	if(identical(sample_info_old,sample_info)==TRUE){warning("WARNING: You did not modify the created sample_info.txt file!\n",call.=F)}
	rm(sample_info_old)}

if(sampleinfo!=FALSE && class(sampleinfo)=="character"){
	sample_info <- read.table(file=sampleinfo,sep="\t",header=TRUE)}

if(sampleinfo!=FALSE && class(sampleinfo)=="data.frame"){
	sample_info <- sampleinfo}

message("Using ",cat(colnames(sample_info[,covars]),sep=" and "),"columns as information for batch effect removal.",appendLF=TRUE)

merged_combat <- merged_norm
#print(sample_info)
if (removeBatcheffect == "EB") {exprs(merged_combat) <- virtualArrayComBat(	expression_xls=exprs(merged_combat),
																			sample_info_file=sample_info,
																			covariates=covars,
																			write=FALSE,
																			prior.plots=FALSE,
																			...)}
rownames(sample_info) <- sample_info[,1]
pData(merged_combat) <- as.data.frame(sample_info)
annotation(merged_combat) <- "org.Hs.eg"
if (removeBatcheffect == "QD") {merged_combat <- normalize.ExpressionSet.qd(merged_combat,...)}
if (removeBatcheffect == "GQ") {merged_combat <- normalize.ExpressionSet.gq(merged_combat,Batch=pData(merged_combat)[,covars[1]],...)}
if (removeBatcheffect == "NORDI") {merged_combat <- normalize.ExpressionSet.nordi(merged_combat,...)}
if (removeBatcheffect == "MRS") {normalize.ExpressionSet.mrs(merged_combat,Batch=pData(merged_combat)[,covars[1]],...)}
if (removeBatcheffect == "MC") {normalize.ExpressionSet.mc(merged_combat,Batch=pData(merged_combat)[,covars[1]],...)}
merged_combat_norm <- normalize.ExpressionSet.quantiles(merged_combat)
#require(org.Hs.eg.db)
#fData(virtualArray$combo_combat_norm.exst) <- virtualArray$annot_longest[match(x=rownames(fData(virtualArray$combo_combat_norm.exst)),table=virtualArray$annot_longest$Gene.Symbol),]
#varMetadata(merged_combat_norm)[,1] <- c("name of the array in the data set; in most cases this corresponds to the file name","further description of the array, e.g the biological source","this designates the array platform or chip type","additional information to group samples together so that biological information isn't lost during batch removal")
message("The numbers of overlapping genes/identifiers between each pair of datasets and the final result were as follows:",appendLF=TRUE)
overlaps <- as.data.frame(calculateOverlaps(rownames_dtfrms),stringsAsFactors=F)
print(overlaps)
if(summary(overlaps[["vs.result"]]<0.5)[3]!=0)
	{
	badDataset <- rownames( overlaps[overlaps[["vs.result"]]<0.5,])
	warning("Datasets",paste(badDataset,sep=" "),"show only little overlap with other datasets and/or result!",call.=F)
	}
return(merged_combat_norm)
}

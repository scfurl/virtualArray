virtualArrayBuildfData <-
function(x,identifier="SYMBOL",collapse_fun=median) {
message("Now preprocessing raw data of ",as.character(x),": Loading annotations...")
identifier <- identifier
z <- NULL
z <- eval(as.symbol(x))
annot <- annotation(z)
if(is.na(GPLs[as.character(annot)])==TRUE)
	{
		message("Using",annot,"as Bioconductor annotation package for dataset",as.character(x),".",appendLF=TRUE)
	}
	else{
		annot2<-GPLs[as.character(annot)]
		message("Detected",annot2,"as Bioconductor annotation package for GEO",annot,"and dataset",as.character(x),".",appendLF=TRUE)
		annot <- annot2
		rm(annot2)
		}
annot_identifier <- paste(annot,identifier,sep="")
annot_package <- sub(" ",".",paste(annot,"db"))
if(is.na((grep(pattern=as.character(annot_package),installed.packages())[1])))
	{
		message("Package",annot_package,"not available. Atempting to install it.",appendLF=TRUE)
		#try(rm(biocLite))
		source("http://www.bioconductor.org/biocLite.R")
		biocContribUrl <- sapply(biocinstallRepos(), contrib.url)
		try(install.packages(as.character(annot_package),method="internal",dependencies=c("Depends","Imports"),repos= biocinstallRepos(),contriburl=biocContribUrl))
		if(try(installed.packages()[as.character(annot_package),1],silent=T)==as.character(annot_package)){message("Package ",as.character(annot_package)," has been successfully installed!",appendLF=TRUE)}else{warning("Package ",as.character(annot_package)," could not be installed! Please install it manually and check for errors!",.call=F)}
	}
require(package=annot_package,character.only=TRUE)
if(class(annot_identifier)[1]=="ProbeAnnDbBimap"){
	annot_identifier <- toggleProbes(value="all",x=eval(as.symbol(annot_identifier)))}
else{annot_identifier <- eval(as.symbol(annot_identifier))}
fData(z) <- toTable(annot_identifier)
#fData(z) <- toTable(eval(as.symbol(annot_identifier)))
#rownames(fData(z)) <- fData(z)[,1]
match_fData_y <- match(fData(z)[,1],rownames(exprs(z)))
# ATTENTION: uncomment lines 21,22 and 34 for some MOUSE ILLUMINA chips
#match_fData_y <- na.omit(match_fData_y)
#fData(z) <- fData(z)[match_fData_y,]
exprs(z) <- as.matrix(exprs(z)[match_fData_y,])
exprs_z <- as.data.frame(exprs(z))
#fData_z <- fData(z)
#rownames(fData_z) <- fData_z[,2]
#rownames(exprs_z) <- fData_z[,2]
#rownames(fData_z) <- rownames(exprs_z)
# ATTENTION: following lines 26,27 are not needed anymore but are left in just in case
#match_fData_y3 <- match(rownames(exprs_z), fData(z)[,2])
#fData(z) <- fData(z)[match_fData_y3,]
exprs_z <- cbind(identifier=fData(z)[,2],exprs_z,row.names=NULL)
#require(reshape)
message("Now preprocessing raw data of ",as.character(x),": Collapsing expression values to their ",as.character(as.list(collapse_fun)[[3]])[2],"...")
# following line because of illumina: some rows with same gene conatin NAs some not, but all are thrown to NAs instead!
exprs_z <- na.omit(exprs_z)
exprs_z <- recast(data=exprs_z,formula=identifier~variable,fun.aggregate=collapse_fun)
message("Now preprocessing raw data of ",as.character(x),": Annotating expression values with ",as.character(identifier),"...")
# ATTENTION: uncomment lines 21,22 and 34 for some MOUSE ILLUMINA chips
#exps_z <- na.omit(exprs_z)
rownames(exprs_z$data) <- exprs_z$labels[[1]][,1]
colnames(exprs_z$data) <- exprs_z$labels[[2]][,1]
exprs_z <- exprs_z$data
exprs(z) <- data.matrix(exprs_z,rownames.force=TRUE)
match_fData_y2 <- match(rownames(exprs(z)),fData(z)[,2])
fData(z) <- fData(z)[match_fData_y2,]
rownames(fData(z)) <- fData(z)[,2]
rownames(exprs(z)) <- fData(z)[,2]
message("Size of expression matrix of ",as.character(x),": ",dim(exprs(z))[1]," rows and ",dim(exprs(z))[2]," columns.")
x <- z
#x <- as.data.frame(exprs(x))
return(x)
}

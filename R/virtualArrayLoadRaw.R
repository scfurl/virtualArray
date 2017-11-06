virtualArrayLoadRaw <-
function(root_dir=getwd(),affy_sum="rma") {
	# PREREQUISITES
	# one folder for each different platform/chip/manufacturer holding .CEL files or other RAW data files		<= CHECK
	# one annotation file for each of these with at least 1 common identifier									<= CHECK
	# their should at least 1 sample be present for all chips/platforms for efficient batch effect removal		<= sample_info.txt automatically created
	virtualArray <- new.env(parent=.GlobalEnv)
	setwd(root_dir)
	cwd <- getwd()
	# generate "sample_info.txt" for use with "virtualArrayComBat" at finish
	cat('\nNow starting to compile "sample_info.txt" for later use with the "virtualArrayComBat" function to remove batch effects.\n\n')
	virtualArray$sample_info <- dir(recursive=TRUE,path="rawdata")
	virtualArray$sample_info <- strsplit(virtualArray$sample_info,"/")
	virtualArray$sample_info <- unlist(virtualArray$sample_info)
	virtualArray$sample_info <- matrix(virtualArray$sample_info,,ncol=3,byrow=TRUE)
	virtualArray$sample_info <- cbind(virtualArray$sample_info[,3],virtualArray$sample_info[,3],virtualArray$sample_info[,2],array(1:(dim(virtualArray$sample_info)[1])))
	colnames(virtualArray$sample_info) <- c("Array name","Sample name","Batch","Covariate 1")
	write.table(virtualArray$sample_info,file="sample_info.txt",sep="\t",row.names=FALSE)
	cat('\nThe file "sample_info.txt" has been written to',root_dir,'.\nPlease edit the "Covariate 1" column to your needs and save before proceeding.\n')
	# load rawdata of underlying directory structure
	cat("\nStarting to load raw data files stored in",root_dir,".\nThis could take some time, please be patient.\n")
	# save current directory for later use
	if (is.loaded("mc_fork", PACKAGE="multicore")){lapply = mclapply}
	#require(GEOquery,quietly=TRUE)
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133Plus2")) > 0) {
		#require(affy,quietly=TRUE)
		setwd("rawdata/Affymetrix/U133Plus2")
		cat("\nLoading Affymetrix HG-U133Plus2 data with justPlier ...\n")
		if (affy_sum == "plier") {
			#require(plier,quietly=TRUE)
			virtualArray$rawdata_U133Plus2 <- ReadAffy()
			virtualArray$rawdata_U133Plus2 <- justPlier(normalize=FALSE,usemm=FALSE,eset=rawdata_U133Plus2)
			}
		if (affy_sum == "rma") {
			virtualArray$rawdata_U133Plus2 <- justRMA(normalize=FALSE)
			}
		if (affy_sum == "gcrma") {
			#require(gcrma,quietly=TRUE)
			virtualArray$rawdata_U133Plus2 <- justGCRMA(normalize=FALSE,type="affinities")
			}
		cat("Done.\n\n")
		setwd(cwd)
		cat('\nNow loading annotation for Affymetrix HG-U133Plus2.\n\n')
		virtualArray$annot_U133Plus2 <- (Table(object=dataTable(getGEO(GEO="GPL570"))))
		virtualArray$annot_U133Plus2$Gene.Symbol <- as.character(virtualArray$annot_U133Plus2$Gene.Symbol)
		virtualArray$annot_U133Plus2$Gene.Symbol <-  lapply(virtualArray$annot_U133Plus2$Gene.Symbol,function(x) {strsplit(x," /// ")})
		virtualArray$annot_U133Plus2$Gene.Symbol <- sapply(virtualArray$annot_U133Plus2$Gene.Symbol,function(x) unlist(x)[1])
	}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U219")) > 0) {
		#require(affy,quietly=TRUE)
		setwd("rawdata/Affymetrix/U219")
		cat("\nLoading Affymetrix HG-U219 data with justPlier ...\n")
		if (affy_sum == "plier") {
			#require(plier,quietly=TRUE)
			virtualArray$rawdata_U219 <- ReadAffy()
			virtualArray$rawdata_U219 <- justPlier(normalize=FALSE,usemm=FALSE,eset=rawdata_U219)
			}
		if (affy_sum == "rma") {
			virtualArray$rawdata_U219 <- justRMA(normalize=FALSE)
			}
		if (affy_sum == "gcrma") {
			#require(gcrma,quietly=TRUE)
			virtualArray$rawdata_U219 <- justGCRMA(normalize=FALSE,type="affinities")
			}
		cat("Done.\n\n")
		setwd(cwd)
		cat('\nNow loading annotation for Affymetrix HG-U219.\n\n')
		virtualArray$annot_U219 <- (Table(object=dataTable(getGEO(GEO="GPL13667"))))
		virtualArray$annot_U219$Gene.Symbol <- as.character(virtualArray$annot_U219$Gene.Symbol)
		virtualArray$annot_U219$Gene.Symbol <-  lapply(virtualArray$annot_U219$Gene.Symbol,function(x) {strsplit(x," /// ")})
		virtualArray$annot_U219$Gene.Symbol <- sapply(virtualArray$annot_U219$Gene.Symbol,function(x) unlist(x)[1])
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133A")) > 0) {
		#require(affy,quietly=TRUE)
		setwd("rawdata/Affymetrix/U133A")
		cat("\nLoading Affymetrix HG-U133A data with justPlier ...\n")
		if (affy_sum == "plier") {
			#require(plier,quietly=TRUE)
			virtualArray$rawdata_U133A <- ReadAffy()
			virtualArray$rawdata_U133A <- justPlier(normalize=FALSE,usemm=FALSE,eset=rawdata_U133A)
			}
		if (affy_sum == "rma") {
			virtualArray$rawdata_U133A <- justRMA(normalize=FALSE)
			}
		if (affy_sum == "gcrma") {
			#require(gcrma,quietly=TRUE)
			virtualArray$rawdata_U133A <- justGCRMA(normalize=FALSE,type="affinities")
			}
		cat("Done.\n\n")
		setwd(cwd)
		cat('\nNow loading annotation for Affymetrix HG-U133A.\n\n')
		virtualArray$annot_U133A <- (Table(object=dataTable(getGEO(GEO="GPL96"))))
		virtualArray$annot_U133A$Gene.Symbol <- as.character(virtualArray$annot_U133A$Gene.Symbol)
		virtualArray$annot_U133A$Gene.Symbol <-  lapply(virtualArray$annot_U133A$Gene.Symbol,function(x) {strsplit(x," /// ")})
		virtualArray$annot_U133A$Gene.Symbol <- sapply(virtualArray$annot_U133A$Gene.Symbol,function(x) unlist(x)[1])
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112F")) > 0) {
		#require(limma,quietly=TRUE)
		setwd("rawdata/Agilent/G4112F")
		virtualArray$rawdata_G4112F_targets <- data.frame(dir())
		colnames(virtualArray$rawdata_G4112F_targets) <- "FileName"
		cat("\nLoading Agilent 4x44k G4112F data ...\n")
		virtualArray$rawdata_G4112F.Elst <- read.maimages(files=virtualArray$rawdata_G4112F_targets$FileName,names=virtualArray$rawdata_G4112F_targets$FileName,source="agilent.median",green.only=TRUE)
		virtualArray$rawdata_G4112F <- new("ExpressionSet",exprs=virtualArray$rawdata_G4112F.Elst$E)
		cat("Done.\n\n")
		setwd(cwd)
		cat('\nNow loading annotation for Agilent 4x44k human G4112F.\n\n')
		virtualArray$annot_G4112F <- (Table(object=dataTable(getGEO(GEO="GPL6480"))))
		virtualArray$annot_G4112F$GENE_SYMBOL <- as.character(virtualArray$annot_G4112F$GENE_SYMBOL)
		virtualArray$annot_G4112F$GENE_SYMBOL <-  lapply(virtualArray$annot_G4112F$GENE_SYMBOL,function(x) {strsplit(x," /// ")})
		virtualArray$annot_G4112F$GENE_SYMBOL <- sapply(virtualArray$annot_G4112F$GENE_SYMBOL,function(x) unlist(x)[1])
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112A")) > 0) {
		#require(limma,quietly=TRUE)
		setwd("rawdata/Agilent/G4112A")
		virtualArray$rawdata_G4112A_targets <- data.frame(dir())
		colnames(virtualArray$rawdata_G4112A_targets) <- "FileName"
		cat("\nLoading Agilent 4x44k G4112A data ...\n")
		virtualArray$rawdata_G4112A.Elst <- read.maimages(files=virtualArray$rawdata_G4112A_targets$FileName,names=virtualArray$rawdata_G4112A_targets$FileName,source="agilent.median",green.only=TRUE)
		virtualArray$rawdata_G4112A <- new("ExpressionSet",exprs=virtualArray$rawdata_G4112A.Elst$E)
		cat("Done.\n\n")
		setwd(cwd)
		cat('\nNow loading annotation for Agilent 4x44k human G4112A.\n\n')
		virtualArray$annot_G4112A <- (Table(object=dataTable(getGEO(GEO="GPL6480"))))
		virtualArray$annot_G4112A$GENE_SYMBOL <- as.character(virtualArray$annot_G4112A$GENE_SYMBOL)
		virtualArray$annot_G4112A$GENE_SYMBOL <-  lapply(virtualArray$annot_G4112A$GENE_SYMBOL,function(x) {strsplit(x," /// ")})
		virtualArray$annot_G4112A$GENE_SYMBOL <- sapply(virtualArray$annot_G4112A$GENE_SYMBOL,function(x) unlist(x)[1])
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanHT12orWG6v3")) > 0) {
		#require(lumi,quietly=TRUE)
		setwd("rawdata/Illumina/HumanHT12orWG6v3")
		cat("\nLoading Illumina Human HT-12 or WG-6 v3.0 data ...\n")
		virtualArray$rawdata_HT12WG6_V3.0_48803_targets <- data.frame(dir())
		colnames(virtualArray$rawdata_HT12WG6_V3.0_48803_targets) <- "FileName"
		virtualArray$rawdata_HT12WG6_V3.0_48803 <- lumiR.batch(as.character(virtualArray$rawdata_HT12WG6_V3.0_48803_targets$FileName),convertNuID=FALSE)
		setwd(cwd)
		virtualArray$sample_info_temp <- read.table(file="sample_info.txt",sep="\t",header=TRUE,colClasses="character")
		virtualArray$sample_info_temp <- virtualArray$sample_info_temp[!virtualArray$sample_info_temp$Batch=="HumanHT12orWG6v3",]
		rownames(virtualArray$sample_info_temp) <- array(1:(dim(virtualArray$sample_info_temp)[1]))
		virtualArray$sample_info_temp[4] <- array(1:(dim(virtualArray$sample_info_temp)[1]))
		virtualArray$sample_info_temp <- rbind(virtualArray$sample_info_temp,data.frame(cbind(Array.name=sampleNames(virtualArray$rawdata_HT12WG6_V3.0_48803),Sample.name=sampleNames(virtualArray$rawdata_HT12WG6_V3.0_48803),Batch="HumanHT12orWG6v3",Covariate.1=c((length(virtualArray$sample_info_temp[,1])+1):(length(virtualArray$sample_info_temp[,1])+length(sampleNames(virtualArray$rawdata_HT12WG6_V3.0_48803)))))))
		write.table(virtualArray$sample_info_temp,file="sample_info.txt",sep="\t",row.names=FALSE)
		#rm(sample_info_tmp)
		cat('\nNow loading annotation for Illumina HumanHT-12_V3.0.\n\n')
		virtualArray$annot_HT12WG6_V3.0_48803 <- (Table(object=dataTable(getGEO(GEO="GPL6947"))))
		virtualArray$annot_HT12WG6_V3.0_48803$Symbol <- as.character(virtualArray$annot_HT12WG6_V3.0_48803$Symbol)
		virtualArray$annot_HT12WG6_V3.0_48803$Symbol <-  lapply(virtualArray$annot_HT12WG6_V3.0_48803$Symbol,function(x) {strsplit(x," /// ")})
		virtualArray$annot_HT12WG6_V3.0_48803$Symbol <- sapply(virtualArray$annot_HT12WG6_V3.0_48803$Symbol,function(x) unlist(x)[1])
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanWG6v2")) > 0) {
		#require(lumi,quietly=TRUE)
		setwd("rawdata/Illumina/HumanWG6v2")
		cat("\nLoading Illumina Human-WG-6 v2.0 ...\n")
		virtualArray$rawdata_WG6_v2.0_48701_targets <- data.frame(dir())
		colnames(virtualArray$rawdata_WG6_v2.0_48701_targets) <- "FileName"
		virtualArray$rawdata_WG6_v2.0_48701 <- lumiR.batch(as.character(virtualArray$rawdata_WG6_v2.0_48701_targets$FileName),convertNuID=FALSE)
		setwd(cwd)
		virtualArray$sample_info_temp <- read.table(file="sample_info.txt",sep="\t",header=TRUE,colClasses="character")
		virtualArray$sample_info_temp <- virtualArray$sample_info_temp[!virtualArray$sample_info_temp$Batch=="HumanWG6v2",]
		rownames(virtualArray$sample_info_temp) <- array(1:(dim(virtualArray$sample_info_temp)[1]))
		virtualArray$sample_info_temp[4] <- array(1:(dim(virtualArray$sample_info_temp)[1]))
		virtualArray$sample_info_temp <- rbind(virtualArray$sample_info_temp,data.frame(cbind(Array.name=sampleNames(virtualArray$rawdata_WG6_v2.0_48701),Sample.name=sampleNames(virtualArray$rawdata_WG6_v2.0_48701),Batch="HumanWG6v2",Covariate.1=c((length(virtualArray$sample_info_temp[,1])+1):(length(virtualArray$sample_info_temp[,1])+length(sampleNames(virtualArray$rawdata_WG6_v2.0_48701)))))))
		write.table(virtualArray$sample_info_temp,file="sample_info.txt",sep="\t",row.names=FALSE)
		#rm(sample_info_tmp)
		cat('\nNow loading annotation for Illumina HumanWG-6_v2.0.\n\n')
		virtualArray$annot_WG6_v2.0_48701 <- (Table(object=dataTable(getGEO(GEO="GPL13376"))))
		virtualArray$annot_WG6_v2.0_48701$Symbol <- as.character(virtualArray$annot_WG6_v2.0_48701$Symbol)
		virtualArray$annot_WG6_v2.0_48701$Symbol <-  lapply(virtualArray$annot_WG6_v2.0_48701$Symbol,function(x) {strsplit(x," /// ")})
		virtualArray$annot_WG6_v2.0_48701$Symbol <- sapply(virtualArray$annot_WG6_v2.0_48701$Symbol,function(x) unlist(x)[1])
		}
	# now all data read and represented as single ExpressionSets and of course NON-NORMALIZED
	cat("\nAll raw data has been loaded.\n\n")
	return(virtualArray)
	}


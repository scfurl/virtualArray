virtualArrayCompile <-
function(root_dir=getwd(),identifier="Gene.Symbol",covars=3,virtualArray=virtualArray) {
	#print(ls(envir=.GlobalEnv,pattern="annot_"))
	virtualArray <- virtualArray
	virtualArray$annot_longest <- eval(as.name(names(sort(sapply(X=ls(envir=virtualArray,pattern="annot_"),FUN=function(x) {dim(eval(as.symbol(x),envir=virtualArray))[1]}),decreasing=TRUE)[1])),envir=virtualArray)
	virtualArray$annot_longest_genesymbols <- virtualArray$annot_longest$Gene.Symbol
	virtualArray$sample_info <- read.table(file="sample_info.txt",sep="\t",header=TRUE)
	cat('\nReading in user supplied "sample_info.txt".\n\n')
	#sample_info <- read.table(file="sample_info.txt",sep="\t",header=TRUE,colClasses="character")
	# transform expression matrices into data.frames
	# collapse rows which target the same gene to their median
	# match common Gene.Symbols of all chips
	#require(reshape,quietly=TRUE)
	if (is.loaded("mc_fork", PACKAGE="multicore")){lapply = mclapply}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133Plus2")) > 0) {
		cat('\nPreparing Affymetrix HG-U133Plus2 data.\n\n')
		virtualArray$rawdata_U133Plus2.dtfrm <- as.data.frame(exprs(virtualArray$rawdata_U133Plus2))
		virtualArray$rawdata_U133Plus2.dtfrm <- cbind(rownames(virtualArray$rawdata_U133Plus2.dtfrm),virtualArray$annot_U133Plus2[11],virtualArray$rawdata_U133Plus2.dtfrm)
		virtualArray$rawdata_U133Plus2_recast.dtfrm <<- recast(data=virtualArray$rawdata_U133Plus2.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_U133Plus2_genesymbols <- virtualArray$rawdata_U133Plus2_recast.dtfrm$Gene.Symbol
		virtualArray$rawdata_U133Plus2_dim.dtfrm <<- dim(virtualArray$rawdata_U133Plus2_recast.dtfrm)[2]
		virtualArray$annot_U133Plus2_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_U133Plus2_recast.dtfrm[, 1])
		virtualArray$rawdata_U133Plus2_subset.dtfrm <<- virtualArray$rawdata_U133Plus2_recast.dtfrm[virtualArray$annot_U133Plus2_genesymbols, 1:virtualArray$rawdata_U133Plus2_dim.dtfrm]
		virtualArray$annot_U133Plus2_genesymbols <- virtualArray$rawdata_U133Plus2_recast.dtfrm$Gene.Symbol
		#colnames(rawdata_U133Plus2.dtfrm)[1] <<- "ID"
	}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133A")) > 0) {
		cat('\nPreparing Affymetrix HG-U133A data.\n\n')
		virtualArray$rawdata_U133A.dtfrm <- as.data.frame(exprs(virtualArray$rawdata_U133A))
		virtualArray$rawdata_U133A.dtfrm <- cbind(rownames(virtualArray$rawdata_U133A.dtfrm),virtualArray$annot_U133A[11],virtualArray$rawdata_U133A.dtfrm)
		virtualArray$rawdata_U133A_recast.dtfrm <<- recast(data=virtualArray$rawdata_U133A.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_U133A_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_U133A_recast.dtfrm[, 1])
		virtualArray$rawdata_U133A_dim.dtfrm <<- dim(virtualArray$rawdata_U133A_recast.dtfrm)[2]
		virtualArray$rawdata_U133A_subset.dtfrm <<- virtualArray$rawdata_U133A_recast.dtfrm[virtualArray$annot_U133A_genesymbols, 1:virtualArray$rawdata_U133A_dim.dtfrm]
		#colnames(rawdata_U133A.dtfrm)[1] <<- "ID"
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U219")) > 0) {
		cat('\nPreparing Affymetrix HG-U219 data.\n\n')
		virtualArray$rawdata_U219.dtfrm <- as.data.frame(exprs(virtualArray$rawdata_U219))
		virtualArray$rawdata_U219.dtfrm <- cbind(rownames(virtualArray$rawdata_U219.dtfrm),virtualArray$annot_U219[15],virtualArray$rawdata_U219.dtfrm)
		virtualArray$rawdata_U219_recast.dtfrm <<- recast(data=virtualArray$rawdata_U219.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_U219_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_U219_recast.dtfrm[, 1])
		virtualArray$rawdata_U219_dim.dtfrm <<- dim(virtualArray$rawdata_U219_recast.dtfrm)[2]
		virtualArray$rawdata_U219_subset.dtfrm <<- virtualArray$rawdata_U219_recast.dtfrm[virtualArray$annot_U219_genesymbols, 1:virtualArray$rawdata_U219_dim.dtfrm]
		#colnames(virtualArray$rawdata_U219.dtfrm)[1] <<- "ID"
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112F")) > 0) {
		cat('\nPreparing Agilent G4112F data.\n\n')
		#virtualArray$rownames(exprs(rawdata_G4112F.exst)) <- virtualArray$rawdata_G4112F$genes$ProbeName
		virtualArray$rawdata_G4112F.dtfrm <- as.data.frame(log2(exprs(virtualArray$rawdata_G4112F))/20*16)
		virtualArray$rawdata_G4112F.dtfrm <- cbind(rownames(virtualArray$rawdata_G4112F.dtfrm),rownames(virtualArray$rawdata_G4112F.dtfrm),virtualArray$rawdata_G4112F.dtfrm)
		virtualArray$rawdata_G4112F.dtfrm[1] <- virtualArray$rawdata_G4112F.Elst$genes$ProbeName
		virtualArray$temp_ids <- match(unlist(virtualArray$annot_G4112F[,1]),virtualArray$rawdata_G4112F.dtfrm[,1])
		virtualArray$temp_dim <- dim(virtualArray$rawdata_G4112F.dtfrm)[2]
		virtualArray$rawdata_G4112F.dtfrm <- virtualArray$rawdata_G4112F.dtfrm[virtualArray$temp_ids, 1:virtualArray$temp_dim]
		virtualArray$rawdata_G4112F.dtfrm[2] <- virtualArray$annot_G4112F[6]
		colnames(virtualArray$rawdata_G4112F.dtfrm)[1:2] <- c("ID","Gene.Symbol")
		#rawdata_G4112F.dtfrm <- cbind(rownames(rawdata_G4112F.dtfrm),annot_G4112F[6],rawdata_G4112F.dtfrm)
		virtualArray$rawdata_G4112F_recast.dtfrm <<- recast(data=virtualArray$rawdata_G4112F.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_G4112F_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_G4112F_recast.dtfrm[, 1])
		virtualArray$rawdata_G4112F_dim.dtfrm <<- dim(virtualArray$rawdata_G4112F_recast.dtfrm)[2]
		virtualArray$rawdata_G4112F_subset.dtfrm <<- virtualArray$rawdata_G4112F_recast.dtfrm[virtualArray$annot_G4112F_genesymbols, 1:virtualArray$rawdata_G4112F_dim.dtfrm]
		#colnames(rawdata_G4112F.dtfrm)[1] <<- "ID"
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112A")) > 0) {
		cat('\nPreparing Agilent G4112A data.\n\n')
		#rownames(exprs(rawdata_G4112A.exst)) <- rawdata_G4112A$genes$ProbeName
		virtualArray$rawdata_G4112A.dtfrm <- as.data.frame(log2(exprs(virtualArray$rawdata_G4112A))/20*16)
		virtualArray$rawdata_G4112A.dtfrm <- cbind(rownames(virtualArray$rawdata_G4112A.dtfrm),rownames(virtualArray$rawdata_G4112A.dtfrm),virtualArray$rawdata_G4112A.dtfrm)
		virtualArray$rawdata_G4112A.dtfrm[1] <- virtualArray$rawdata_G4112A.Elst$genes$ProbeName
		virtualArray$temp_ids <- match(unlist(virtualArray$annot_G4112A[,1]),virtualArray$rawdata_G4112A.dtfrm[,1])
		virtualArray$temp_dim <- dim(virtualArray$rawdata_G4112A.dtfrm)[2]
		virtualArray$rawdata_G4112A.dtfrm <- virtualArray$rawdata_G4112A.dtfrm[virtualArray$temp_ids, 1:virtualArray$temp_dim]
		virtualArray$rawdata_G4112A.dtfrm[2] <- virtualArray$annot_G4112A[6]
		colnames(virtualArray$rawdata_G4112A.dtfrm)[1:2] <- c("ID","Gene.Symbol")
		#rawdata_G4112A.dtfrm <- cbind(rownames(rawdata_G4112A.dtfrm),annot_G4112A[6],rawdata_G4112A.dtfrm)
		virtualArray$rawdata_G4112A_recast.dtfrm <<- recast(data=virtualArray$rawdata_G4112A.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_G4112A_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_G4112A_recast.dtfrm[, 1])
		virtualArray$rawdata_G4112A_dim.dtfrm <<- dim(virtualArray$rawdata_G4112A_recast.dtfrm)[2]
		virtualArray$rawdata_G4112A_subset.dtfrm <<- virtualArray$rawdata_G4112A_recast.dtfrm[virtualArray$annot_G4112A_genesymbols, 1:virtualArray$rawdata_G4112A_dim.dtfrm]
		#colnames(rawdata_G4112A.dtfrm)[1] <<- "ID"
		}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanHT12orWG6v3")) > 0) {
		cat('\nPreparing Illumina Human HT-12 or WG-6 v3.0 data.\n\n')
		#rawdata_HT12WG6_V3.0_48803.dtfrm <- as.data.frame(exprs(rawdata_HT12WG6_V3.0_48803))
		virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm <- as.data.frame(log2(exprs(virtualArray$rawdata_HT12WG6_V3.0_48803)))
		virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm <- cbind(rownames(virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm),virtualArray$annot_HT12WG6_V3.0_48803[15],virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm)
		virtualArray$temp_ids <- match(unlist(virtualArray$annot_HT12WG6_V3.0_48803[,1]),virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm[,1])
		virtualArray$temp_dim <- dim(virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm)[2]
		virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm <- virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm[virtualArray$temp_ids, 1:virtualArray$temp_dim]
		virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm[2] <- virtualArray$annot_HT12WG6_V3.0_48803[13]
		colnames(virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm)[1:2] <- c("ID","Gene.Symbol")
		virtualArray$rawdata_HT12WG6_V3.0_48803_recast.dtfrm <<- recast(data=virtualArray$rawdata_HT12WG6_V3.0_48803.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_HT12WG6_V3.0_48803_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_HT12WG6_V3.0_48803_recast.dtfrm[, 1])
		virtualArray$rawdata_HT12WG6_V3.0_48803_dim.dtfrm <<- dim(virtualArray$rawdata_HT12WG6_V3.0_48803_recast.dtfrm)[2]
		virtualArray$rawdata_HT12WG6_V3.0_48803_subset.dtfrm <<- virtualArray$rawdata_HT12WG6_V3.0_48803_recast.dtfrm[virtualArray$annot_HT12WG6_V3.0_48803_genesymbols, 1:virtualArray$rawdata_HT12WG6_V3.0_48803_dim.dtfrm]
	}
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanWG6v2")) > 0) {
		cat('\nPreparing Illumina Human WG-6 v2.0 data.\n\n')
		#rawdata_WG6_v2.0_48701.dtfrm <- as.data.frame(exprs(rawdata_WG6_v2.0_48701))
		virtualArray$rawdata_WG6_v2.0_48701.dtfrm <- as.data.frame(log2(x=exprs(virtualArray$rawdata_WG6_v2.0_48701)-min(exprs(virtualArray$rawdata_WG6_v2.0_48701))))
		virtualArray$rawdata_WG6_v2.0_48701.dtfrm <- cbind(rownames(virtualArray$rawdata_WG6_v2.0_48701.dtfrm),virtualArray$annot_WG6_v2.0_48701[15],virtualArray$rawdata_WG6_v2.0_48701.dtfrm)
		virtualArray$temp_ids <- match(unlist(virtualArray$annot_WG6_v2.0_48701[,1]),virtualArray$rawdata_WG6_v2.0_48701.dtfrm[,1])
		virtualArray$temp_dim <- dim(virtualArray$rawdata_WG6_v2.0_48701.dtfrm)[2]
		virtualArray$rawdata_WG6_v2.0_48701.dtfrm <- virtualArray$rawdata_WG6_v2.0_48701.dtfrm[virtualArray$temp_ids, 1:virtualArray$temp_dim]
		virtualArray$rawdata_WG6_v2.0_48701.dtfrm[2] <- virtualArray$annot_WG6_v2.0_48701[13]
		colnames(virtualArray$rawdata_WG6_v2.0_48701.dtfrm)[1:2] <- c("ID","Gene.Symbol")
		virtualArray$rawdata_WG6_v2.0_48701_recast.dtfrm <<- recast(data=virtualArray$rawdata_WG6_v2.0_48701.dtfrm,formula=Gene.Symbol~variable,median)
		virtualArray$annot_WG6_v2.0_48701_genesymbols <- match(virtualArray$annot_longest_genesymbols,virtualArray$rawdata_WG6_v2.0_48701_recast.dtfrm[, 1])
		virtualArray$rawdata_WG6_v2.0_48701_dim.dtfrm <<- dim(virtualArray$rawdata_WG6_v2.0_48701_recast.dtfrm)[2]
		virtualArray$rawdata_WG6_v2.0_48701_subset.dtfrm <<- virtualArray$rawdata_WG6_v2.0_48701_recast.dtfrm[virtualArray$annot_WG6_v2.0_48701_genesymbols, 1:virtualArray$rawdata_WG6_v2.0_48701_dim.dtfrm]
	}
	cat('\nPutting single chip data together.\n\n')
	virtualArray$combo.dtfrm_nrow <- dim(eval(as.name(ls(pattern="subset",name=virtualArray)[1]),envir=virtualArray))[1]
	virtualArray$combo.dtfrm <- data.frame(virtualArray$annot_longest$Gene.Symbol)
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133Plus2")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_U133Plus2_subset.dtfrm[2:virtualArray$rawdata_U133Plus2_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U133A")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_U133A_subset.dtfrm[2:virtualArray$rawdata_U133A_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="U219")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_U219_subset.dtfrm[2:virtualArray$rawdata_U219_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112F")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_G4112F_subset.dtfrm[2:virtualArray$rawdata_G4112F_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="G4112A")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_G4112A_subset.dtfrm[2:virtualArray$rawdata_G4112A_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanHT12orWG6v3")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_HT12WG6_V3.0_48803_subset.dtfrm[2:virtualArray$rawdata_HT12WG6_V3.0_48803_dim.dtfrm]) }
	if (length(grep(x=as.character(virtualArray$sample_info[,3]),pattern="HumanWG6v2")) > 0) {
			virtualArray$combo.dtfrm <- cbind(virtualArray$combo.dtfrm,virtualArray$rawdata_WG6_v2.0_48701_subset.dtfrm[2:virtualArray$rawdata_WG6_v2.0_48701_dim.dtfrm]) }
	cat('\nDone putting single chip data together.\n\n')
	virtualArray$combo_noNA.dtfrm <- na.omit(virtualArray$combo.dtfrm)
	virtualArray$combo_noNA.dtfrm <- unique(virtualArray$combo_noNA.dtfrm)
	rownames(virtualArray$combo_noNA.dtfrm) <- virtualArray$combo_noNA.dtfrm$annot_longest.Gene.Symbol
	virtualArray$combo_noNA.dtfrm <- virtualArray$combo_noNA.dtfrm[2:dim(virtualArray$combo_noNA.dtfrm)[2]]
	virtualArray$combo.matrix <- data.matrix(virtualArray$combo_noNA.dtfrm,rownames.force=TRUE)
	virtualArray$combo.exst <- new("ExpressionSet",exprs=virtualArray$combo.matrix)
	cat('\nApplying quantile normalization to whole dataset.\n\n')
	#require(affyPLM,quietly=TRUE)
	virtualArray$combo_norm.exst <- normalize.ExpressionSet.quantiles(virtualArray$combo.exst)
	#boxplot(exprs(combo_norm.exst))
	#dev.off()
	virtualArray$combo_combat.exst <- virtualArray$combo_norm.exst
	cat('\nRemoving batch effects with "virtualArrayComBat".\n\n')
	exprs(virtualArray$combo_combat.exst) <- virtualArrayComBat(covariates=covars,expression_xls=exprs(virtualArray$combo_norm.exst),par.prior=TRUE,prior.plots=FALSE,sample_info_file=virtualArray$sample_info,skip=0,write=FALSE)
	virtualArray$combo_combat_norm.exst <- normalize.ExpressionSet.quantiles(virtualArray$combo_combat.exst)
	rownames(virtualArray$sample_info) <- virtualArray$sample_info[,1]
	pData(virtualArray$combo_combat_norm.exst) <- as.data.frame(virtualArray$sample_info)
	fData(virtualArray$combo_combat_norm.exst) <- virtualArray$annot_longest[match(x=rownames(fData(virtualArray$combo_combat_norm.exst)),table=virtualArray$annot_longest$Gene.Symbol),]
	varMetadata(virtualArray$combo_combat_norm.exst)[,1] <- c("name of the array in the data set; in most cases this corresponds to the file name","further description of the array, e.g the biological source","this designates the array platform or chip type","additional information to group samples together so that biological information isn't lost during batch removal")
	cat('\nAll done!\n\n')
	return(virtualArray$combo_combat_norm.exst)
	}


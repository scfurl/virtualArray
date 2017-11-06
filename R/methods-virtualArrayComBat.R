setMethod("virtualArrayComBat", "matrix",
{
virtualArrayComBat <- function(expression_xls, sample_info_file, type='txt', write=FALSE, covariates="Batch", par.prior=TRUE, filter=FALSE, skip=0, prior.plots=FALSE){
	# apply ComBat.R for batch effect removal and write result to file
	#write.table(exprs(combo_norm.exst),file="110720_test_AFX+AG1.txt",sep="\t")
	#source(file="http://jlab.byu.edu/ComBat/Download_files/ComBat.R")
	#ComBat(covariates="all",expression_xls="110720_test_AFX+AG1.txt",par.prior=TRUE,prior.plots=TRUE,sample_info_file="sample_info.txt",skip=1,write=FALSE)
	# 'expression_xls' is the expression index file (e.g. outputted by dChip); 'sample_info_file' is a tab-delimited text file containing the colums: Array  name, sample name, Batch, and any other covariates to be included in the modeling; 'type' currently supports two data file types 'txt' for a tab-delimited text file and 'csv' for an Excel .csv file (sometimes R handles the .csv file better, so use this if you have problems with a .txt file!); 'write' if 'T' ComBat writes adjusted data to a file, and if 'F' and ComBat outputs the adjusted data matrix if 'F' (so assign it to an object! i.e. NewData <- ComBat('my expression.xls','Sample info file.txt', write=FALSE)); 'covariates=all' will use all of the columns in your sample info file in the modeling (except array/sample name), if you only want use a some of the columns in your sample info file, specify these columns here as a vector (you must include the Batch column in this list); 'par.prior' if 'T' uses the parametric adjustments, if 'F' uses the nonparametric adjustments--if you are unsure what to use, try the parametric adjustments (they run faster) and check the plots to see if these priors are reasonable; 'filter=value' filters the genes with absent calls in > 1-value of the samples. The defaut here (as well as in dchip) is .8. Filter if you can as the EB adjustments work better after filtering. Filter must be numeric if your expression index file contains presence/absence calls (but you can set it >1 if you don't want to filter any genes) and must be 'F' if your data doesn't have presence/absence calls; 'skip' is the number of columns that contain probe names and gene information, so 'skip=5' implies the first expression values are in column 6; 'prior.plots' if true will give prior plots with black as a kernal estimate of the empirical batch effect density and red as the parametric estimate. 
	#debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=FALSE; covariates='all'; par.prior=TRUE; filter=FALSE; skip=0; prior.plots=TRUE
	cat('Reading Sample Information File\n')
	#saminfo <- read.table(sample_info_file, header=TRUE, sep='\t',comment.char='')
	saminfo <- sample_info_file
	if(sum(colnames(saminfo)=="Batch")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
		
	cat('Reading Expression Data File\n')
	dat <- expression_xls
#	if(type=='csv'){
#		dat <- read.csv(expression_xls,header=TRUE,as.is=TRUE)
#                #print(dat[1:2,])
#	#	dat <- dat[,trim.dat(dat)]  
#                #print(colnames(dat))
#		colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=TRUE)[1:ncol(dat)]
#                #print(colnames(dat))
#		}
#         else{
#		dat <- read.table(expression_xls,header=TRUE,comment.char='',fill=TRUE,sep='\t', as.is=TRUE)
#		dat <- dat[,trim.dat(dat)]
#		colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=TRUE)[1:ncol(dat)]
#		}
	if (skip>0){
              geneinfo <- as.matrix(dat[,1:skip])
              dat <- dat[,-c(1:skip)]}
        else{geneinfo=NULL}
        #print(geneinfo[1:4])
        #print(dat[1:2,])
	
	if(filter){
		ngenes <- nrow(dat)
		col <- ncol(dat)/2
		present <- apply(dat, 1, filter.absent, filter)
		dat <- dat[present, -(2*(1:col))]
		if (skip>0){geneinfo <- geneinfo[present,]}
		cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
		}

	if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
	tmp <- match(colnames(dat),saminfo[,1])
	if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
	#tmp1 <- match(saminfo[,1],colnames(dat))
	#saminfo <- saminfo[tmp1[!is.na(tmp1)],]
	saminfo <- saminfo[tmp,]  ## Bug fixed 01/04/2011		

	#if(any(covariates != 'all')){  }
	saminfo <- saminfo[,c("Array.name","Sample.name",covariates)]
	design <- design.mat(saminfo)	


	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(dat))
	if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
        #print(dat[1:2,])
	##Standardize Data across genes
	cat('Standardizing Data across genes\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=TRUE)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=TRUE))
		}

	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}


	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
	if(write){
		output_file <- paste('Adjusted',expression_xls,'.xls',sep='_')
                 #print(geneinfo[1:2])
                 #print(bayesdata[1:2,1:4])
		 #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
		#suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE))
                outdata <- cbind(ProbeID=geneinfo, bayesdata); write.table(outdata, file=output_file, sep="\t")
		cat("Adjusted data saved in file:",output_file,"\n")
		}else{return(cbind(geneinfo,bayesdata))}
	}
})

setMethod("virtualArrayComBat", "character",
function(expression_xls, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
{
	dat <- read.csv(expression_xls,header=TRUE,as.is=TRUE)
	callGeneric(expression_xls = dat, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
})

setMethod("virtualArrayComBat", "ExpressionSet",
function(expression_xls, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
{
	dat <- exprs(expression_xls)
	callGeneric(expression_xls = dat, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
})

setMethod("virtualArrayComBat", "data.frame",
function(expression_xls, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
{
	dat <- as.data.matrix(expression_xls)
	callGeneric(expression_xls = dat, sample_info_file, type = "txt", write = FALSE, covariates = "Batch", par.prior = TRUE, filter = FALSE, skip = 0, prior.plots = FALSE)
})

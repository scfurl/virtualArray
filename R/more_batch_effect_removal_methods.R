# more batch effect removal methods from WebArrayDB and CONOR R package

normalize.ExpressionSet.mc <- function(ExpressionSet=NULL, Batch=NULL, ...) {
	# ExpressionSet: an ExpressionSet object, needs a pData slot including "Batch" column!
	Batch_internal <- NULL
	Batch <- Batch
	M <- exprs(ExpressionSet)
	if (length(Batch)>=0) 
		{
			Batch_internal <- Batch
		}
	if (length(Batch)==0) 
		{
		if (length(pData(ExpressionSet)$Batch)>=0) 
		{
			Batch_internal <- pData(ExpressionSet)$Batch
		}
		else
		{
			cat("No batches found! Quitting.\n")
			return(NULL)
		}
		}
  mdmtrx <- model.matrix( ~ factor(Batch_internal) - 1)
  M <- M - mc.misreg.simple(mdmtrx, M)
  exprs(ExpressionSet) <- M
  return(ExpressionSet)
}


mc.misreg.simple <- function(mdmtrx, M) {
# M is a microarray data matrix
# mdmtrx is the model or response matrix
  naM <- is.na(M)
  nsamples <- (!naM)%*%mdmtrx
  M[naM] <- 0
  Msum <- M%*%mdmtrx
  Mbar <- Msum/nsamples
  Mbar %*% t(mdmtrx)
}


# QD: quantile discretization normalization; 
normalize.ExpressionSet.qd <- function(ExpressionSet=NULL, nbin=8, ...) { # quantile discretization algorithm
		#lapply(seq(ncol(M)), function(i) {ra <- rank(M[,i]); rlt<-as.integer(ra/(length(ra)/nbin)); M[,i] <- rlt - median(rlt)} )
		#lapply(seq(ncol(M)), function(i) {ra <- rank(M[,i]); rlt<-ceiling(ra/(length(ra)/nbin)); M[,i] <- rlt - median(rlt)} )
		# ExpressionSet: an ExpressionSet object
		M <- exprs(ExpressionSet)
		ref <- c(-1*((nbin/2)-1):0, 0:((nbin/2)-1))
		ref <- ref[ceiling(seq(nrow(M))/(nrow(M)/nbin))]
		lapply(seq(ncol(M)), function(i) M[,i] <<- ref[rank(M[,i])])
		exprs(ExpressionSet) <- M
		return(ExpressionSet)
	}

# NorDi: NORmal DIscretization; http://keia.i3s.unice.fr/?Logiciels:Nordi
nordi1p_getoutliers<-function(D,pvalue){
		Out=TRUE;
		Noa=TRUE;
		while (Out && Noa){
			J<-jarque.bera.test(D);
			J0<-J$statistic;
			G<-grubbs.test(D,type=10);
			DBack<-D;
			if (G$p.value<pvalue) {
				D<-rm.outlier(D, fill = FALSE, median = FALSE)
			}else{
				Out<-FALSE
			};
		
			J1<-jarque.bera.test(D);
			J1<-J1$statistic;
			if (J0-J1>=0) {
				Noa<-TRUE
			}else{
				D<-DBack; 
				Noa<-FALSE
			};
		}; 
		return(D)
	};
normalize.ExpressionSet.nordi = function(ExpressionSet=NULL, pvalue = .01, alpha = .05, ...){
		#The pvalue determines the sensitivity for detecting outliers 
		#in each column of the gene expression matrix.
		#The alpha value determines the size of the tails of the 
		#normal distributions considered to be over or under expressed.
		# ExpressionSet: an ExpressionSet object, needs a pData slot including "Batch" column!
		M <- exprs(ExpressionSet)
		C<-split(M,col(M));


		Ot<- c(0*1:length(C));
		Ut<- c(0*1:length(C));
		for (i in 1:length(C)){
			C[[i]]<-na.omit(C[[i]]);
			C[[i]]<-nordi1p_getoutliers(C[[i]],pvalue);
			m<-mean(C[[i]]);
			stdev<-sqrt(var(C[[i]]));
			zcutoff = -1*qnorm(alpha/2)
			Ot[i]<-m + zcutoff*stdev;
			Ut[i]<-m - zcutoff*stdev;
		}
		Mfin<-M*0;    
		for(i in 1:length(colnames(M))){
			for (j in 1:length(rownames(M))){
				ifelse(M[j,i]<Ut[i],Mfin[j,i]<-(-1),ifelse(M[j,i]>Ot[i],Mfin[j,i]<-(1),Mfin[j,i]<-0));
			}
		}
		exprs(ExpressionSet) <- Mfin
		return(ExpressionSet)
	}

# MRS: median rank scores; 
normalize.ExpressionSet.mrs <- function(ExpressionSet=NULL, Batch=NULL, ...) { # median rank  scores algorithm  - a modification of quantile for multi-platform data
	# ExpressionSet: an ExpressionSet object, needs a pData slot including "Batch" column!
    # Batch: a vector (character or numeric) containing batch contribution of individual samples when not included in ExpressionSets pData!
    # Batch=c( rep("affy",12),rep("lumi",6) )
	Batch_internal <- NULL
	Batch <- Batch
	M <- exprs(ExpressionSet)
	if (length(Batch)>=0) 
		{
			Batch_internal <- Batch
		}
	if (length(Batch)==0) 
		{
		if (length(pData(ExpressionSet)$Batch)>=0) 
		{
			Batch_internal <- pData(ExpressionSet)$Batch
		}
		else
		{
			cat("No batches found! Quitting.\n")
			return(NULL)
		}
		}
	idx <- tapply(seq(Batch_internal), Batch_internal, function(x) x)
	cat("Following batches were found:\n")
	print(idx)
	if (length(Batch_internal)<=0) return(M)
	imax <- which.max(sapply(idx, length)) # for reference
	ref_med_srt <- sort(apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE) ))
	idx[imax] <- NULL
	lapply(unlist(idx), function(i) M[,i] <<- ref_med_srt[rank(M[,i])])
    exprs(ExpressionSet) <- M
	return(ExpressionSet)
	}

# GQ: "Gene Quantile" cross platform normalization algorithm; 
normalize.ExpressionSet.gq <- function(ExpressionSet=NULL, Batch=NULL, ...) { 
	#This function was provided by Xiao-Qin Xia, one of the authors of 
	#webarraydb.
    # modified MRS
    # ExpressionSet: an ExpressionSet object, needs a pData slot including "Batch" column!
    # Batch: a vector (character or numeric) containing batch contribution of individual samples when not included in ExpressionSets pData!
    # Batch=c( rep("affy",12),rep("lumi",6) )
	Batch_internal <- NULL
	Batch <- Batch
	M <- exprs(ExpressionSet)
	if (length(Batch)>=0) 
		{
			Batch_internal <- Batch
		}
	if (length(Batch)==0) 
		{
		if (length(pData(ExpressionSet)$Batch)>=0) 
		{
			Batch_internal <- pData(ExpressionSet)$Batch
		}
		else
		{
			cat("No batches found! Quitting.\n")
			return(NULL)
		}
		}
	idx <- split(seq(Batch_internal), Batch_internal)
	cat("Following batches were found:\n")
	print(idx)
    imax <- which.max(sapply(idx, length)) # for reference
    ref_med <- apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE))
    ref_med_srt <- sort(ref_med)
    idx[imax] <- NULL
    lapply(idx, function(i) {
         MTMP <- sapply(i, function(x) ref_med_srt[rank(M[,x])]); 
         M[,i] <<- MTMP - apply(MTMP, 1, median) + ref_med 
         } )
    exprs(ExpressionSet) <- M
	return(ExpressionSet)
}

# QN: quantile normalization; 
#normalize.ExpressionSet.quantiles(ExpressionSet)

# XPN: cross platform normalization; 
##xpn = function(platform1.data,platform2.data,K=10,L=4,p1.names=0,p2.names=0,gene.cluster="kmeans",assay.cluster="kmeans",corr="pearson", iterations=30, skip.match=FALSE){
#	#If K or L is not a single value, it is taken as a list of possible values.
#	
#	#Match names
#	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
#	x <- input$x
#	y <- input$y
#	input <- NA
#	
#	#Remove the medians
#	x_demed = x - rowMedians(as.matrix(x))
#	y_demed = y - rowMedians(as.matrix(y))
#	
#	#Get the dimensions
#	nx=dim(x)[2]
#	ny=dim(y)[2]
#	mx=dim(x)[1]
#	my=dim(y)[1]
#	
#
#	#Create the combined dataframe for clustering purposes
#	combined <- cbind(x_demed,y_demed)
#	
#	#Detect K and L if necessary
#	K = detect_K(combined,K,corr)
#	L = detect_L(x_demed,y_demed,L,corr)
#
#	
#	xout = 0
#	yout = 0
#	
#
#	for(iter in 1:iterations){
#	cat("XPN iteration",iter,"out of",iterations,"\n")
#
#	#Do the assay clustering
#	assayclusters = xpnassaycluster(x_demed, y_demed, L, assay.cluster, corr)
#
#	#Do the gene clustering
#	geneclusters = xpngenecluster(combined, K, gene.cluster, corr)
#
#	#Estimate the XPN parameters
#	xparams = xpnparamsp(x,geneclusters,assayclusters[1:nx])
#	yparams = xpnparamsp(y,geneclusters,assayclusters[(nx+1):(nx+ny)])
#	
#	#Calcuate the weighted averages
#	nor = 1/(nx + ny)
#	Aav = try((1/(xparams$nblock+yparams$nblock))*(xparams$nblock*xparams$A + yparams$nblock*yparams$A))
#	bav = nor*(nx*xparams$b + ny*yparams$b)
#	cav = nor*(nx*xparams$c + ny*yparams$c)
#	sigmaav = sqrt(nor*(nx*xparams$s2 + ny*yparams$s2))
#	sigmax = sqrt(xparams$s2)
#	sigmay = sqrt(yparams$s2)
#
#	
#	#Calculate the expanded A
#	expAx = xpnclusterexpand(xparams$A,geneclusters,assayclusters[1:nx])
#	expAy = xpnclusterexpand(yparams$A,geneclusters,assayclusters[(nx+1):(nx+ny)])
#
#	#Calculate the residuals
#	epsilonx = as.matrix(x) - (as.vector(xparams$b) * expAx + kronecker(ones(1,nx),xparams$c))
#	epsilony = as.matrix(y) - (as.vector(yparams$b) * expAy + kronecker(ones(1,ny),yparams$c))
#	
#	#Calculate the expanded average A
#	expAavx = xpnclusterexpand(Aav,geneclusters,assayclusters[1:nx])
#	expAavy = xpnclusterexpand(Aav,geneclusters,assayclusters[(nx+1):(nx+ny)])
#
#	#Calculate the output values 
#	xout = xout + (1/iterations)*((as.vector(bav) * expAavx) + kronecker(ones(1,nx),cav) + as.vector(sigmaav/sigmax) * epsilonx)
#	yout = yout + (1/iterations)*((as.vector(bav) * expAavy) + kronecker(ones(1,ny),cav) + as.vector(sigmaav/sigmay) * epsilony)
#
#	}#end of the enclosing for loop
#	
#	#Put the rownames back in and convert to data frames
#	xout = as.data.frame(xout,row.names=rownames(x))
#	yout = as.data.frame(yout,row.names=rownames(y))
#	
#	#All done!
#	return(list(x=xout,y=yout))
#}
#xpnparamsp = function(x, geneclusters, assayclusters, method){
#	#x should be a dataframe or hashdataframe	
#
#	#Get K and L
#	K=max(geneclusters)
#	L=max(assayclusters)
#	
#	#Get number of assays and genes
#	numassays=length(x)
#	numgenes=length(rownames(x))
#	
#	#Number of assays in each cluster
#	nj = matrix(table(as.factor(assayclusters)),1,L)
#	
#	#Set up the output variables
#	A = matrix(nrow=K,ncol=L)
#	nblock = matrix(nrow=K,ncol=L)
#	b = matrix(nrow=numgenes, ncol=1)
#	c = matrix(nrow=numgenes, ncol=1)
#	sigma = matrix(nrow=numgenes, ncol=1)
#	mu = matrix(nrow=numgenes,ncol=L)
#	
#
#	#For each gene cluster, estimate the XPN parameters
#	for (a in 1:K){
#
#		#Get the logical indexing vector for gene cluster a 
#		geneinds = geneclusters==a
#		numgenesa = sum(as.numeric(geneinds))
#		
#		#Get the data for this gene cluster
#		xa = as.matrix(x[geneinds,])
#		
#		#Get the number of genes in this cluster
#		Ga = dim(xa)[1]
#		S = dim(xa)[2]
#		
#		#For each assay cluster and each gene, get the (relatively) unconstrained
#		#MLE for mu
#		mua = matrix(nrow=Ga,ncol=L)
#		for (bb in 1:L){
#			assayinds = assayclusters==bb
#			xab = matrix(xa[,assayinds],numgenesa,sum(as.numeric(assayinds)))
#			mua[,bb]=t(t(rowMeans(xab)))
#		}
#		
#		
#		#Calculate sigma2
#		expmua = xpnclusterexpand(mua, colclusters = assayclusters)
#		sigma2 = xa - expmua
#		sigma2 = sigma2*sigma2
#		sigma2 = xpnclustercollapse(sigma2,colclusters = assayclusters)
#		
#		solna = XPN_MLE_C(mua,sigma2,nj)
#	
#		ba = solna$b
#		Aa = solna$A
#		ca = solna$c
#		s2a = solna$s2
#
#		#Write the values for this gene group into the larger arrays
#		sigma[geneinds,] = s2a
#		mu[geneinds,] = mua
#		A[a,] = Aa
#		nblock[a,] = nj
#		b[geneinds,] = ba
#		c[geneinds,] = ca
#	}
#
#	return(list(A=A,nblock=nblock,b=b,c=c,s2=sigma,mu=mu))
#	
#}
#xpnclustercollapse = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2])){
#	#This undoes the work of xpnclusterexpand
#	
#	l = dim(aMatrix)[1]
#	w = dim(aMatrix)[2]
#	
#	#collapse the rows
#	rowidx = split(1:l,rowclusters)
#	lc = length(rowidx)
#	collapse1 = matrix(NA,lc,w)
#	for ( i in 1:lc){
#		collapse1[i,] =  colMeans(matrix(aMatrix[rowidx[[i]],],length(rowidx[[i]]),w))
#	}
#	
#	#collapse the columns
#	colidx = split(1:w,colclusters)
#	wc = length(colidx)
#	collapse2 = matrix(NA,lc,wc)
#	for ( i in 1:wc){
#		collapse2[,i] =  rowMeans(matrix(collapse1[,colidx[[i]]],lc,length(colidx[[i]])))
#	}
#	
#	return(collapse2)
#	
#}
#xpnclusterexpand = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2]), noise = 0){
#	#This is basically a fancy repmat.  rowclusters and colclusters 
#	#as produced by pam.  aMatrix must have the same number of rows
#	#and columns as there are row and column clusters.
#	
#	l = dim(aMatrix)[1]
#	w = dim(aMatrix)[2]
#	
#	#Get the dimensions of the output
#	m = length(rowclusters)
#	n = length(colclusters)
#	
#	#Get the number of row and column clusters
#	nrowclust = max(rowclusters)
#	ncolclust = max(colclusters)
#	
#	#initialize the vertically expanded matrix
#	vert = zeros(m,dim(aMatrix)[2])
#	
#	#Do the vertical expansion
#	for (i in 1:m){
#		vert[i,]=aMatrix[rowclusters[i],] + (noise*randn(1,w))
#	}
#
#	#initialize the fully expanded matrix
#	output = matrix(NA,m,n)
#	
#	#Do the horizontal expansion
#	for (j in 1:n){
#		output[,j] = vert[,colclusters[j]] + (noise*randn(m,1))
#	}
#	
#	#That is it!
#	return(output)
#}
#xpnassaycluster = function(x,y,L,cluster="kmeans",corr="pearson"){
#
#	if(is.numeric(cluster) | is.factor(cluster)){
#		return(as.numeric(cluster))
#	}
#
#	#The numbers of assays
#	nx = dim(x)[2]
#	ny = dim(y)[2]
#	
#	#The number of genes
#	m = dim(x)[1]
#	
#	#Compute the two correlation matrices if necessary
#	if (cluster == 'pam' | (length(L)>1)){
#		xdiss = 1 - cor(as.matrix(x), method = corr)
#		ydiss = 1 - cor(as.matrix(y), method = corr)
#	}
#	
#	#Determine L if a range has been given
#	if(length(L)>1){
#		Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
#		Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
#		L = max(Lx,Ly)
#		cat(as.character(L), "assay clusters detected...\n")
#	}
#	
#	#Do the clustering
#	if (cluster == "classic"){
#		
#		done = FALSE
#		while (!done){
#			xyclust <- ckmeans(t(cbind(x,y)),centers=L, iter.max=20, distance=corr)
#			#print(xyclust$cluster)
#			xclust <- xyclust$cluster[1:nx]
#			yclust <- xyclust$cluster[(nx+1):(nx+ny)]
#			#print(xclust)
##			print(yclust)
##			print(nlevels(as.factor(xclust)))
##			print(nlevels(as.factor(yclust)))
#			done <- nlevels(as.factor(xclust)) == L & nlevels(as.factor(yclust)) == L
#			if(!done){
#				cat("Not all platforms contained all assay clusters.  Trying again (only \"classic\" mode has this problem)...\n")
#			}
#		}
#		return(xyclust$cluster)
#	}
#	else if (cluster == "pam"){
#		xclust = pam(xdiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
#		yclust = pam(ydiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
#	}else if (cluster == "kmeans"){
#
#		#This loop ensures there are no empty clusters.
#		done = FALSE
#		while (!done){
#			xclust = ckmeans(t(x), centers=L, iter.max = 20, distance = corr)
#			yclust = ckmeans(t(y), centers=L, iter.max = 20, distance = corr)
#
#			done = !xclust$empty & !yclust$empty
#		}
#		xclust = xclust$cluster
#		yclust = yclust$cluster	
#
#	} else if (cluster == "flexclust"){
#
#		distCor1 = function (x, centers) 
#		{
#    		z <- matrix(0, nrow(x), ncol = nrow(centers))
#    		for (k in 1:nrow(centers)) {
#        		z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
#    		}
#    		z
#		}
#		#This loop ensures there are no empty clusters.
#		done = FALSE
#		while (!done){
#		xclust = kcca(t(x), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)
#		yclust = kcca(t(y), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)
#			
#			done = (dim(xclust@clusinfo)[1] == L) &(dim(yclust@clusinfo)[1] == L) 
#		}
#		xclust = xclust@cluster
#		yclust = yclust@cluster	
#	}else if(cluster == "random"){
#			xclust = sample(c(1:L,sample(1:L,nx-L,replace = TRUE)),nx,replace=FALSE)
#			yclust = sample(c(1:L,sample(1:L,ny-L,replace = TRUE)),ny,replace=FALSE)
#	}
#
#	#Compute cluster averages
#	xave = matrix(NA,m,L)
#	yave = matrix(NA,m,L)
#	for (i in 1:L){
#		xinds = xclust==i
#		yinds = yclust==i
#		xave[,i] = rowMeans(as.matrix(x[,xinds],nrow=m,ncol=sum(as.numeric(xinds))))
#		yave[,i] = rowMeans(as.matrix(y[,yinds],nrow=m,ncol=sum(as.numeric(yinds))))
#	}
#
#	#Compute the cluster correlation matrix
#	clustercor = cor(xave,yave, method = corr)
#	
#	#Map clusters
#	xtoymap = matrix(NA,L,1)
#	for (i in 1:L){
#		highest = which.max(clustercor)
#		xind = highest%%L
#		
#		yind = highest%/%L + 1
#		
#		if (xind==0){
#			xind = L
#			yind = yind-1
#		}
#		xtoymap[xind]=yind
#		clustercor[xind,] = -2
#		clustercor[,yind] = -2
#	}
#
#	#Change the x clusters using the map
#	newxclust = xclust
#	for (i in 1:L){
#		newxclust[xclust==i] = xtoymap[i]
#	}
#	xclust = newxclust
#
#	#Return the clustering in a combined vector
#	return(c(xclust,yclust))
#}
#xpngenecluster = function(x,K,cluster="pam",corr="pearson"){
#	
#	#Dimensions
#	m = dim(x)[1]
#	n = dim(x)[2]
#	
#	#Compute the dissimilarity matrix
#	if (cluster=="pam"|(length(K)>1)){
#		xdiss = 1 - cor(t(x),method=corr)
#	}
#	
#	#Determine K if a range was given
#	if(length(K)>1 ){
#		pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
#		K = pamklist$nc
#		genepamobj = pamklist$pamobject
#	}else{
#		genepamobj = NULL	
#	}
#	
#	
#	#Do the gene clustering
#	if(cluster=="pam"){
#		if (!is.null(genepamobj)){
#			geneclusters = genepamobj$clustering
#		}else{
#			geneclusters = pam(xdiss, k=K, diss=TRUE, keep.diss=FALSE, keep.data=FALSE, cluster.only=TRUE)
#		}
#	}else if(cluster=="kmeans"){
#		#This loop ensures there are no empty clusters.
#		done = FALSE
#		while (!done){
#			geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
#			done = !geneclusters$empty			
#		}
#		geneclusters = geneclusters$cluster
#	}else if(cluster=="flexclust"){
#		
#		
#		distCor1 = function (x, centers) #change
#		{
#    		z <- matrix(0, nrow(x), ncol = nrow(centers))
#    		for (k in 1:nrow(centers)) {
#        		z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
#    		}
#    		z
#		}
#		#This loop ensures there are no empty clusters.
#		done = FALSE
#		while (!done){
#			geneclusters = kcca(x, k=K, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE) #change
#			#geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
#			done = (dim(geneclusters@clusinfo)[1]==K)#change		
#		}
#		geneclusters = geneclusters@cluster
#	}else if (cluster == "random"){
#		geneclusters = sample(c(1:K,sample(1:K,m-K,replace = TRUE)),m,replace=FALSE)
#	}else{
#		cat("Unknown gene clustering method:",cluster,"\n")
#	}
#	
#	return(geneclusters)
#	
#}
#XPN_MLE_C = function(xbar,sigma2,nj){
#	
#	
#	n = sum(nj)
#	I = dim(xbar)[1]
#	J = dim(xbar)[2]
#	
#	A = zeros(1,J)
#	b = ones(I,1)
#	c = zeros(I,1)
#	s2 = ones(I,1)
#	
#	cout = .C("XPN_MLE_C_C",as.double(xbar),A=as.double(A),b=as.double(b),c=as.double(c),s2=as.double(s2),as.double(sigma2),as.integer(I),as.integer(J),as.integer(n),as.integer(nj))
#	
#	A = matrix(cout$A,1,J)
#	b = matrix(cout$b,I,1)
#	c = matrix(cout$c,I,1)
#	s2 = matrix(cout$s2,I,1)
#	
#	return(list(A=A,b=b,c=c,s2=s2))
#}
#XPN_MLE = function(xbar,sigma2,nj){
#	#This function is never used.  It serves to illustrate the method of maximum likelihood estimation used by XPN.  It has been replaced by XPN_MLE_C, which calls a more efficient C function.  The inputs and outputs of these functions are identical.
#	
#	n = sum(nj)
#	I = dim(xbar)[1]
#	J = dim(xbar)[2]
#	
#	A = zeros(1,J)
#	b = ones(I,1)
#	c = zeros(I,1)
#	s2 = ones(I,1)
#	
#	old = c(A, b, c, s2)*0
#	iter = 0
#	while(sum((c(A, b, c, s2)-old)^2)>(1e-16)*max(old)){
#		
#		print(iter<-iter+1)
#		print(sum((c(A, b, c, s2)-old)^2))
#		old = c(A, b, c, s2)
#		
#		c = matrix(rowSums((xbar-repmat(b,1,J)*repmat(A,I,1))*repmat(nj,I,1))/n,I,1)
#		
#		if(sum(b)<0){
#			b = -1*b;
#		}
#		
#		A = matrix(colSums(repmat(b,1,J)*(xbar-repmat(c,1,J))/repmat(s2,1,J)/sum(b^2/s2)),1,J)
#		
#		A = A - mean(A)
#		A = A * sqrt(J/sum(A^2))
#		
#		b = matrix(rowSums(repmat(A,I,1)*(xbar-repmat(c,1,J))*repmat(nj,I,1))/sum(A^2 * nj),I,1)
#		
#		
#		s2 = matrix(rowSums(((xbar-repmat(c,1,J)-repmat(A,I,1)*repmat(b,1,J))^2 + sigma2)*repmat(nj,I,1))/n,I,1)
#		
#		
#		
#		s2[s2==0] = 2.2251e-308
#		
#		
#		
#	}
#	
#	return(list(A=A,b=b,c=c,s2=s2))
#}
#detect_L = function(x,y,L,corr="pearson"){
#	
#	#Determine L if a range has been given
#	if(length(L)>1 & length(L) != (dim(x)[2]+dim(y)[2])){
#		xdiss = 1 - cor(as.matrix(x), method = corr)
#		ydiss = 1 - cor(as.matrix(y), method = corr)
#		Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
#		Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
#		L = max(Lx,Ly)
#		cat(as.character(L), "assay clusters detected...\n")
#	}
#	return(L)
#}
#detect_K = function(x,K,corr="pearson"){
#
#	#Determine K if a range was given
#	if(length(K)>1 & length(K) != dim(x)[1]){
#		xdiss = 1 - cor(t(x),method=corr)
#		pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
#		K = pamklist$nc
#		genepamobj = pamklist$pamobject
#	}else{
#		genepamobj = NULL	
#	}
#	return(K)
#}
#
## EB: empirical Bayes methods; 
## see ComBat => already implemented!
#
## DWD: Distance Weighted Discrimination;
#sepelimdwd = function(Xp,Xn,penalty, useSparse=TRUE){
#	#Xp and Xn are matrices.  penalty is a scalar.  This is an adaptation of the Matlab function of the same name, written by J. S. Marron, available at https://genome.unc.edu/pubsup/dwd/.
#	
#	flag <- 0
#	#Dimensions of the data
#	dp <- dim(Xp)[1]
#	np <- dim(Xp)[2]
#	dn <- dim(Xn)[1]
#	nn <- dim(Xn)[2]
#	if (dn != dp) {stop('The dimensions are incomapatible.')}
#	d <- dp
#
#	#Dimension reduction in HDLSS setting.
#	XpnY <- as.matrix(cbind(Xp,-1*Xn))
#	XpnY11 <- XpnY[1,1]
#	n <- np + nn
#	if(d>n){
#		qrfact = qr(XpnY)
#		Q = qr.Q(qrfact)
#		R = qr.R(qrfact)
#		RpnY = R
#		dnew = n
#	}else{
#		RpnY = XpnY
#		dnew = d
#	}
#
#	y = ones(np + nn,1)
#	y[(np+1):(np+nn),1] = -1
#	ym = y[2:n,1]
#	
#	# nv is the number of variables (eliminating beta)
#	# nc is the number of constraints
#	nv = 1 + dnew + 4*n
#	nc = 2*n
#	#Set up the block structure, constraint matrix, rhs, and cost vector
#	
#	blk = list()
#	blk$type = character()
#	blk$size = list()
#	blk$type[1] = 's'
#	blk$size[[1]] = cbind(dnew+1,3*ones(1,n))
#	blk$type[2] = 'l'
#	blk$size[[2]] = n
#	
#	
#	Avec = list()
#	A = zeros(nc,nv-n)
#	col1 = RpnY[,1]
#	A[1:(n-1),2:(dnew+1)] = t(RpnY[,2:n] - col1%*%t(ym))
#	A[1:(n-1),seq(dnew+5,dnew+1+3*n,3)]	= -1*speye(n-1)
#	A[1:(n-1),seq(dnew+6,dnew+2+3*n,3)]	= speye(n-1)
#	A[1:(n-1),dnew+2] = ym
#	A[1:(n-1),dnew+3] = -1*ym
#	A[n,1] = 1
#	A[(n+1):(n+n),seq(dnew+4,dnew+3+3*n,3)] = speye(n)
#	
#
#	Avec[[1]] = t(A)
#	Avec[[2]] = (rbind(cbind(-1*ym,speye(n-1)),zeros(1+n,n)))
#	b = rbind(zeros(n-1,1),ones(1+n,1))	
#
#	
#	C = list()
#	c = zeros(nv-n,1)
#	c[seq(dnew+2,dnew+1+3*n,3),1] = ones(n,1)
#	c[seq(dnew+3,dnew+2+3*n,3),1] = ones(n,1)
#	
#	
#	C[[1]] = c
#	C[[2]] = penalty*ones(n,1)
#	
#	
#	#Solve the SOCP problem
#
###################CLSOCP############################
#	CL_K <- unlist(blk$size)
#	CL_qlen <-sum(blk$size[[1]])
#	CL_llen <-sum(blk$size[[2]])
#	CL_type <- c(rep('q',length(blk$size[[1]])),rep('l',length(blk$size[[2]])))
#	CL_A <- cbind(t(Avec[[1]]),Avec[[2]])
#	CL_c <- rbind(C[[1]],C[[2]])
#	
#
#	soln <- socp(CL_A,b,CL_c,CL_K,CL_type,gamma_fac=.3,sigma0 = .1,use_sparse=useSparse)
#	
#	
#	X1 <- soln$x[1:CL_qlen]
#	X2 <- soln$x[(CL_qlen+1):(CL_qlen+CL_llen)]
#	lambda <- soln$y
######################################################
#
#
#	# Compute the normal vector w and constant term beta.
#
#	barw = X1[2:(dnew+1)]
#	if (d>n){
#		w = Q %*% barw
#	}else{
#		w = barw	
#	}
#	beta = X1[dnew + 2] - X1[dnew + 3] - X2[1] - t(col1)%*%barw
#	normw = norm(w)
#	if (normw < 1 - 1e-3){
#		print(normw)
#	}
#	normwm1 = 0
#	if (normw > 1 - 1e-3){
#		w = w/normw
#		normwm1 = norm(w) - 1
#		beta = beta/normw
#	}
#
#	
#	# Compute the minimum of the supposedly positive 
#	# and the maximum of the supposedly negative residuals.
#	# Refine the primal solution and print its objective value.
#	
#	residp = t(Xp) %*% w + beta[1] #optimization
#	residn = t(Xn) %*% w + beta[1] #optimization
#	minresidp = min(residp)
#	maxresidn = max(residn)
#	res = t(XpnY) %*% w + beta[1] * y
#	rsc = 1/sqrt(penalty)
#	xi = -1* res + rsc[1]
#	xi[xi<0] <- 0
#	totalviolation = sum(xi)
#	minresidpmod = min(residp + xi[1:np])
#	maxresidnmod = max(residn - xi[(np+1):n])
#	minxi = min(xi)
#	maxxi = max(xi)
#	resn = res + xi
#	rresn = 1 / resn
#	primalobj = penalty * sum(xi) + sum(rresn)
##	print(primalobj)
#	
#
#	#Compute the dual solution alp and print its objective value.
#	alp = zeros(n,1)
#	lambda1 = lambda[1:(n-1)]
#	alp[1] = -1*t(ym)%*%lambda1
#	alp[2:n] = lambda1
#	alp = alp * (as.numeric(alp>0))
#	sump = sum(alp[1:np])
#	sumn = sum(alp[(np+1):n])
#	sum2 = (sump + sumn)/2	
#	alp[1:np] = (sum2/sump)*alp[1:np]
#	alp[(np+1):n] = (sum2/sumn)*alp[(np+1):n]
#	maxalp = max(alp)
#	if (maxalp > penalty | maxxi > 1e-3){
#		alp = (penalty[1]/maxalp)*alp
#	}
#	minalp = min(alp)
#	p = RpnY%*%alp
#	eta = -1*norm(p)
#	gamma = 2*sqrt(alp)
#	dualobj = eta + sum(gamma)
##	print(dualobj)
#	
#
#	
#	#dualgap is a measure of the accuracy of the solution
#	dualgap = primalobj - dualobj
##	print(dualgap)
#	
#	if (dualgap > 1e-4){
#		flag = -1
#	}
#	
#	
#	return(list(w=w,beta=beta,residp=residp,residn=residn,alp=alp,totalviolation=totalviolation,dualgap=dualgap,flag=flag))
#	
#}
#DWD1SM = function(trainp,trainn,threshfact = 100, useSparse = TRUE){
#
#	np = dim(trainp)[2]
#	nn = dim(trainn)[2]
##	vpwdist2 = numeric(np*nn)
#	vpwdist2x <- rdist(t(trainp),t(trainn))
##	for (ip in 1:np){
##		vpwdist2[((ip-1)*nn+1):(ip*nn)] <- colSums((trainp[,ip] - trainn)^2) #optimization
##	}
##	medianpwdist2 = median(vpwdist2)
#	medianpwdist2 = median(vpwdist2x)^2
#
#	penalty = threshfact / medianpwdist2
#	sepelimout = sepelimdwd(trainp,trainn,penalty,useSparse=useSparse)
#	w = sepelimout$w
#	flag = sepelimout$flag
#	if (flag == -1){	
#		cat("Inaccurate solution!\n")
#	}
#	if (flag == -2){
#		cat("Infeasible or unbounded optimization problem!\n")
#	}
#	dirvec = w/norm(w)
#	return(dirvec)
#}
#dwd = function(platform1.data, platform2.data, platform1.train = NULL, platform2.train=NULL, p1.names=0, p2.names=0,p1.train.names=0, p2.train.names=0, skip.match=FALSE, use.sparse = TRUE) {
#
#	#Match names
#	if(is.null(platform1.train) & is.null(platform2.train)){
#		input = processplatforms(list(x=platform1.data,y=platform2.data), namesvec = c(p1.names, p2.names), skip.match=skip.match)
#
#		
#		
#		#Get normal vector for DWD adjustment
#		dirvec = DWD1SM(input[[1]],input[[2]],useSparse=use.sparse)
#		#Project the data
#		vprojp = t(input[[1]]) %*% dirvec
#		vprojn = t(input[[2]]) %*% dirvec
#		meanprojp = mean(vprojp)
#		meanprojn = mean(vprojn)
#		output = list()
#		p1.adjust <- -1 * meanprojp * dirvec 
#		p2.adjust <- -1 * meanprojn * dirvec 
#		names(p1.adjust) <- rownames(input[[1]])
#		names(p2.adjust) <- rownames(input[[2]])
#		output$x = input[[1]] + p1.adjust
#		output$y = input[[2]] + p2.adjust
#		output$p1.adjust <- p1.adjust
#		output$p2.adjust <- p2.adjust
#		
#	}else{
#		input = processplatforms(list(x=platform1.data,y=platform2.data,platform1.train,platform2.train), namesvec = c(p1.names, p2.names, p1.train.names, p2.train.names), skip.match=skip.match)
#	
#	
#		#Get normal vector for DWD adjustment
#		dirvec = DWD1SM(input[[3]],input[[4]],useSparse=use.sparse)
#		
#		#Project the data
#		vprojp = t(input[[3]]) %*% dirvec
#		vprojn = t(input[[4]]) %*% dirvec
#
#		meanprojp = mean(vprojp)
#		meanprojn = mean(vprojn)
#	
#		output = list()
#		p1.adjust <- -1 * meanprojp * dirvec 
#		p2.adjust <- -1 * meanprojn * dirvec 
#		names(p1.adjust) <- rownames(input[[1]])
#		names(p2.adjust) <- rownames(input[[2]])
#		output$x = input[[1]] + p1.adjust
#		output$y = input[[2]] + p2.adjust
#		output$p1.adjust <- p1.adjust
#		output$p2.adjust <- p2.adjust
#	}
#		
#	return(output)
#	
#}
#norm = function(aMatrix){
##Returns the largest singular value of aMatrix
#	o = svd(aMatrix,nu=0,nv=0)
#	return(o$d[1])
#	
#}
#
## distran: Distribution Transformation; 
#distran = function(platform1.data, platform2.data, p1.assayclasses=NULL, L=4, cluster="pam",corr="pearson", p1.names=0, p2.names=0, skip.match=FALSE){
#	#platform1.data is used as a reference set, which is why p1.sampleinds is present.  
#	#p1.assayclass should be a vector with the same number of elements as there are
#	#assays in platform1.data.  If not given, the same clustering procedure used by xpn
#	#will be used to automatically compute assay clusters, giving an estimate of which 
#	#assays come from the same sample or sample class.  Additional parameters apply to 
#	#the xpn sample clustering procedure.
#	
#	#Match names
#	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
#	
#	genenames <- rownames(input$x)
#	colnamesx <- colnames(input$x)
#	colnamesy <- colnames(input$y)
#
#	#dimensions
#	m = dim(input$x)[1]
#	nx = dim(input$x)[2]
#	ny = dim(input$y)[2]
#	
#	#Cluster if necessary
#	if (is.null(p1.assayclasses)){
#		xpnclusters = xpnassaycluster(input$x,input$y,L=L,cluster=cluster,corr=corr)
#		p1.assayclasses = xpnclusters[1:nx]
#	}
#
#	
#	#Construct the reference data
#	classes = unique(p1.assayclasses)
#	nclasses = length(classes)
#	refdat = numeric(m)
#	for (class in classes){
#		idx = p1.assayclasses==class
#		refdat = refdat + (1/nclasses)*rowMeans(as.matrix(input$x[,idx],m,length(idx)))
#	}
#	refdat = sort(refdat)
#	
#	#Get the ranks of the data
#	input$x = colRanks(input$x)
#	input$y = colRanks(input$y)
#	
#	#Do the substitution
#	input$x = as.matrix(input$x)
#	input$y = as.matrix(input$y)
#	for (i in 1:(m*ny)){
#		input$y[i] = refdat[input$y[i]]
#	}
#	for (i in 1:(m*nx)){
#		input$x[i] = refdat[input$x[i]]
#	}
#	
#	colnames(input$x) <- colnamesx
#	colnames(input$y) <- colnamesy
#	rownames(input$x) <- genenames
#	rownames(input$y) <- genenames
#
#	
#	#All done!
#	return(list(x=data.frame(input$x),y=data.frame(input$y)))
#	
#}
#colRanks = function(aDataFrame){return(apply(aDataFrame,2,rank))}
#
virtualArrayBuildSampleInfo <-
function(x){
								message("Now combining pData of all data sets...")
								sample_info <- NULL
								single_sample_info <- sapply(x,virtualArrayBuildSingleSampleInfo,simplify=FALSE)
								#print(single_sample_info)
								#for (i in 1:length(single_sample_info)) {
								#	sample_info <- sample_info
								#	sample_info<-rbind(sample_info,single_sample_info[[i]])}
								p1 = single_sample_info[[1]]
								for(i in 2:length(single_sample_info))
								   {
									p2 = single_sample_info[[i]]
									cp = intersect(colnames(p1),colnames(p2))
									tp = unique(union(colnames(p1),colnames(p2)))

									sp1 = setdiff(colnames(p1),cp)
									sp2 = setdiff(colnames(p2),cp)

									pheno = matrix(NA,ncol=length(tp),nrow=nrow(p1)+nrow(p2))
									rownames(pheno) = c(rownames(p1),rownames(p2))
									colnames(pheno) = tp;
									
									if(length(cp)!=0){
										pheno[1:nrow(p1),cp] = as.matrix(p1[,cp])
										pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)),cp] = as.matrix(p2[,cp])}
									
									if(length(sp1)!=0){
										pheno[1:nrow(p1),sp1] = as.matrix(p1[,sp1])}
									if(length(sp2)!=0){
										pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)), sp2] = as.matrix(p2[,sp2])}
									
									pData = as.data.frame(pheno)
									p1 <- pData
								   }
								sample_info <- p1
								#sample_info <- cbind(sample_info,Covariate.1=c(1:nrow(sample_info)))
								return(sample_info)
								}


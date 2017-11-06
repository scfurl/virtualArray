virtualArrayBuildSingleSampleInfo <-
function(x){
								message("Now setting up pData of ",as.character(x),"...")
								temp <- eval(as.symbol(x))
								sample_names <-sampleNames(temp)
								chip_names <- as.character(x)
								single_sample_info <- cbind(Array.name=sample_names,
															Sample.name=sample_names,
															Batch=chip_names,
															pData(temp))
								#print(single_sample_info)
								sample_info<-NULL
								return(single_sample_info)
								}


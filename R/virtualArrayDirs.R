virtualArrayDirs <-
function(root_dir=getwd()) {
	# create the directory structure suitable for and expected by the array comparison scripts
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Affymetrix/U219")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Affymetrix/U133A")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Affymetrix/U133Plus2")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Agilent/G4112F")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Agilent/G4112A")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Illumina/HumanHT-12_V3.0")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Illumina/HumanRef-8_v2.0")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Illumina/HumanWG-6_v2.0")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Illumina/HumanWG-6_v3.0")
	dir.create(recursive=TRUE,showWarnings=FALSE,path="rawdata/Illumina/HumanRef-8_v3.0")
	cat("\nCreated directory tree in",root_dir,". \nPlease copy your raw data files into the appropriate directories.\n\n")
	}


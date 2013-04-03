setwd("/Users/seanmaguire/Desktop/spring2013_Bioinformatics/bisulfate/f12/")
library(methylKit)
#returns list of .sam files

################################some functions####################################
#finds files that have a certain pattern in the file name for example ('.sam')
fileList.fnc<-function(pattern){
	fileList<-vector()
	for(files in list.files()){
		if(grepl(pattern,files))(fileList<-c(fileList,files))
	}
	return(as.list(fileList))
}
#here we need to sort the files using unix. THIS WILL NOT WORK ON WINDOWS. only needs to run once. pretty slow.
sortSAM<-function(fileList){
for(file in fileList){
	newFile<-paste("sorted.",file,sep="")
	unix<-paste("grep -v '^[[:space:]]*@'",file,"| sort -k3,3 -k4,4n  >",newFile,sep=" ")
	system(unix)
}
}
sampleNames<-function(SortedFileList){
	names<-vector()
	for(file in fileList){
		name<-sub("sorted.","",file)
		name<-sub("bismark_pe.sam","",name)
		names<-c(names,name)
	}
return(as.list(names))
}
#######################################################################################
fileList<-fileList.fnc(".sam")
sortSAM(fileList)

newFileList<-fileList.fnc("sorted.")
names<-sampleNames(newFileList)
read.bismark(newFileList,sample.id=names,assembly='hg18',treatment=as.list(c(rep("test",length(newFileList)))))


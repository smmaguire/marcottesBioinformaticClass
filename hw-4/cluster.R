setwd("/Users/seanmaguire/Desktop/spring2013_Bioinformatics/hw-4")
list.files()
data=read.table("protien.csv",sep='\t',header=TRUE)
names(data)[1]<-'pro_name'
# Ward Hierarchical Clustering
d <- dist(data[,2:ncol(data)], method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
#plot(fit,labels=data[,1]) # display dendogram
dendro<-as.dendrogram(fit)
plot(cut(dendro, h=3)$upper,)
clust<-cutree(fit,h=3)
data$pro_name[clust==7]

wss <- (nrow(data[,2:ncol(data)])-1)*sum(apply(data[,2:ncol(data)],2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data[,2:ncol(data)], 
  	 centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares")

fit2<-kmeans(data[,2:ncol(data)],7,iter.max=1000,nstart=1000)
data <- data.frame(data, fit2$cluster)
library(cluster)
clusplot(data[,2:21], fit2$cluster, color=TRUE, shade=TRUE, 
         labels=0, lines=0)

###phylo
phylo<-read.table("yeast_phyloprofiles.txt",sep=" ",header=TRUE)
d2 <- dist(phylo, method = "euclidean") # distance matrix
fit3 <- hclust(d2, method="ward")
dendro3<-as.dendrogram(fit3)
plot(cut(dendro3, h=100)$upper,)
clust2<-cutree(fit3,h=100)
row.names(phylo)[clust2==21]

###micro
expVec<-read.delim("yeast_microarraydata.txt",header=FALSE)
expVec<-expVec[,1:301]
#rows 1040,1281,3225,3795,5092 have ONLY NAs. Have to ommitted
badRows<-c(1040,1281,3225,3795,5092)
goodRows<-c(1:1039,1041:1280,1282:3224,3226:3794,3796:5091,5093:length(expVec[,1]))
expVec<-expVec[goodRows,]
library(cluster)
d3 <- clara(expVec[,2:301],k=10, metric = "euclidean") # distance matrix
plot(d3,which.plots=1)
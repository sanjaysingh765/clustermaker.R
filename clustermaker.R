
library(RColorBrewer)
library(gplots)   
library(ggplot2)
library(reshape2)



#upload raw data
rawdata<- read.csv("expression",header=TRUE,row.names=1)
head(rawdata)
dim(rawdata)


# Remove all gene which has 0 value in all sample
all <- apply(rawdata, 1, function(x) all(x==0) )
newdata <- rawdata[!all,]



# remove uninformative genes keep only genes that are expressed in at least 10 RPKM in 1 tissue
filter <- newdata[rowSums(newdata > 1) >= 1,]
head (filter)
dim(filter)

# Generates row and column dendrograms.
hr <- hclust(as.dist(1-cor(t(filter), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(filter, method="pearson")), method="complete") 
   

# Cuts the tree and creates color vector for clusters.
mycl <- cutree(hr, h=max(hr$height)/1.5); 
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); 
nofclust.height <-  length(unique(as.vector(mycl)));
selcol <- colorRampPalette(brewer.pal(12,"Set3"))
selcol2 <- colorRampPalette(brewer.pal(9,"Set1"))
clustcol.height = selcol2(nofclust.height);
mycolhc <- clustcol.height[mycl]

#create a heatmap with all possible cluster   
jpeg("all_cluster.jpeg", units="in", family="Times New Roman",  width=5, height=5, res=300, pointsize = 9) #pointsize is font size| increase image size to see the key
heatmap.2(as.matrix(filter), 
Rowv=as.dendrogram(hr), 
Colv=as.dendrogram(hc), 
#Colv=NULL,
col=greenred(75), 
scale="row", 
density.info="none", 
trace="none", 
#RowSideColors=mycolhc,
RowSideColors=mycolhc,
srtCol = 45,
labRow = FALSE,
margins=c(6,8),
cexRow = 0.75, 
cexCol=0.75, 
lhei = c(0.21,1),
lwid=c(2,4.5),
key.title=NA)
dev.off()


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
x11(height=6, width=2); 
names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order]))
dev.copy2pdf(file="key.pdf")

# Select sub-cluster number clid=c(1,2) and generate corresponding dendrogram.
clid <- c(1); #change cluster number
ysub <- filter[names(mycl[mycl%in%clid]),]; 
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete")
 

#create a heatmap for  selected cluster
jpeg("cluster1.jpeg", units="in", family="Times New Roman",  width=5, height=4, res=300, pointsize = 9) 
# Create heatmap for chosen sub-cluster.
heatmap.2(
as.matrix(ysub), 
Rowv=as.dendrogram(hrsub), 
#Colv=as.dendrogram(hc), 
col=greenred(75), 
#Colv=FALSE,
#Rowv=FALSE,
#dendrogram="none", 
scale="row", 
density.info="none", 
trace="none", 
RowSideColors=mycolhc[mycl%in%clid],
srtCol = 45,
#labRow = T,
margins=c(8,12),
cexRow = 0.75, 
cexCol=0.75, 
lhei = c(0.31,1),
lwid=c(2,4.5),
key.title=NA)

dev.off()



#make a boxplot of  selected cluster
logy <- log2(ysub+1)
myDF <- cbind(Row.Names = rownames(logy), logy)
head(myDF)
dfm <- melt(myDF,id.vars="Row.Names")
colors = brewer.pal(10, "Set3")

jpeg("cluster_4_boxplot.jpeg", units="in", family="Times New Roman",  width=7, height=5, res=300, pointsize = 9) #pointsize is font size| increase image size to see the key
ggplot(data=dfm) + 
geom_boxplot(aes(x=variable,y=value),fill = colors)+
#geom_boxplot(aes(x=variable,y=value))+
theme(legend.position = "bottom") +
guides(fill = guide_legend(
title.theme = element_text(size=12, face="bold", angle = 45)))+
 theme(legend.title=element_text(size=12, face="bold", angle = 45))+
  theme(axis.title.y = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10))+
theme(axis.text.x = element_text(face="plain", color="black", size=10, angle=45,, vjust = 1, hjust=1),
          axis.text.y = element_text(face="plain", color="black", size=10, angle=360))+
scale_y_continuous(expand=c(0,0), limits = c(0, 8)) + 
  xlab("Tissue") + 
ylab("Expression (log2 FPKM)")+ 
theme(plot.margin=unit(c(2,2,2,2),"cm"))
dev.off()


# Print out row labels in same order as shown in the heatmap for selected cluster.
test <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
write.table(test, "cluster_4.txt", sep="\t", append=F, quote=F, row.names=T)

# Print out all cluster and respective genes in a file
write.table(mycl, "All_cluster.txt", sep="\t", append=F, quote=F, row.names=T)






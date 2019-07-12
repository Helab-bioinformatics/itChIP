library(Matrix)
library(proxy)
library(gplots)

######filterering
df<- read.table("merge.countmatrix.sort.peak.txt",header = T,row.names = 1,sep = "\t")

###assign the selected genes to nearby peaks.
#marker <-read.table("marker.txt",header = F)
marker <-read.csv("select_feature_gene/germ_layers_ESC_refgene.csv",header = F)
#which(duplicated(marker$V1))
marker <- unique(as.character(marker$V1))
anno <- read.delim("peaks_anno.txt",header = T,sep = "\t")
anno.gene <- as.character(anno$Gene.Name)

anno.gene.use <- na.omit(anno.gene.use)
region.use <- paste(anno.gene.use$Chr,anno.gene.use$Start-1,anno.gene.use$End,sep = "_")
gene.region.match <- rbind(data.frame(region=region.use,anno.gene.use$Gene.Name))


#########plot heatmap of selected H3K27ac peak signal 
marker.all <-read.csv("select_feature_gene/germ_layers_ESC_refgene.csv",header = F)
mesoderm <- marker.all[1:149,]
endoderm <- marker.all[150:349,]
ectoderm <- marker.all[350:420,]
ESC <- marker.all[421:616,]

df.ESC <- df[as.character(gene.region.match$region)[na.omit(match(as.character(ESC$V1),gene.region.match$anno.gene.use.Gene.Name))],]#180
df.mesoderm <- df[as.character(gene.region.match$region)[na.omit(match(as.character(mesoderm$V1),gene.region.match$anno.gene.use.Gene.Name))],]#120
df.endoderm <- df[as.character(gene.region.match$region)[na.omit(match(as.character(endoderm$V1),gene.region.match$anno.gene.use.Gene.Name))],]#164
df.ectoderm <- df[as.character(gene.region.match$region)[na.omit(match(as.character(ectoderm$V1),gene.region.match$anno.gene.use.Gene.Name))],]#64


########agg gene(peak) signature
mesoderm.count <- apply(df.mesoderm,2,mean)
endoderm.count <- apply(df.endoderm,2,mean)
ectoderm.count <- apply(df.ectoderm,2,mean)
ESC.count <- apply(df.ESC,2,mean)
a<- rbind(ESC.count,mesoderm.count,endoderm.count,ectoderm.count)
a.nor<- t(t(a)/colSums(a))
m<-a.nor
scale_max <- 2
scale_min <- -2

m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
plot(density(as.matrix(m)))

m[is.nan(m)] = 0
m[m>scale_max] = scale_max
m[m<scale_min] = scale_min
plot(density(as.matrix(m)))

library(gplots)
#cols.ident<-c(rainbow(9))
hm <- heatmap.2(as.matrix(m),col="bluered", scale="none", dendrogram="col", labCol=F,Rowv=F, cexRow=1, Colv=T,density.info="none", trace="none")




#######
mesoderm.use <-data.frame(meso=m[2,match(mesoderm.index,colnames(m))],ESC=m[1,match(mesoderm.index,colnames(m))],celltype =celltype.mesoderm)
endoderm.use <-data.frame(endo=m[3,match(endoderm.index,colnames(m))],ESC=m[1,match(endoderm.index,colnames(m))],celltype =celltype.endoderm)
ectoderm.use <-data.frame(ecto=m[4,match(ectoderm.index,colnames(m))],ESC=m[1,match(ectoderm.index,colnames(m))],celltype =celltype.ectoderm) 

mesoderm.use <- mesoderm.use[order(mesoderm.use$ESC,decreasing = T),]
endoderm.use <- endoderm.use[order(endoderm.use$ESC,decreasing = T),]
ectoderm.use <- ectoderm.use[order(ectoderm.use$ESC,decreasing = T),]


plot(mesoderm.use$ESC,mesoderm.use$meso,col=rainbow(7)[factor(mesoderm.use$celltype)],pch=16,cex=0.6,xlim=c(-2,4),ylim=c(-2,4))
lines(lowess(mesoderm.use$meso,mesoderm.use$ESC))
legend("topright", 
       col=rainbow(7), 
       lty=c(1,1,1,1,1),
       cex=1,
       lwd=4) 


plot(endoderm.use$ESC,endoderm.use$endo,col=rainbow(6)[factor(endoderm.use$celltype)],pch=16,cex=0.6,xlim=c(-2,4),ylim=c(-2,4))
lines(lowess(endoderm.use$endo,endoderm.use$ESC))
legend("topright", 
       col=rainbow(6), 
       lty=c(1,1,1,1,1),
       cex=1,
       lwd=4) 


plot(ectoderm.use$ESC,ectoderm.use$ecto,col=rainbow(6)[factor(ectoderm.use$celltype)],pch=16,cex=0.6,xlim=c(-2,4),ylim=c(-2,4))
lines(lowess(ectoderm.use$ecto,ectoderm.use$ESC))
legend("topright", 
       col=rainbow(6), 
       lty=c(1,1,1,1,1,1),
       cex=1,
       lwd=4)


#projection the signal near feature gene
dim <-mesoderm.use[intersect(colnames(df.select.nor),rownames(mesoderm.use)),]
#dim <-ectoderm.use[intersect(colnames(df.select.nor),rownames(ectoderm.use)),]
#dim <-endoderm.use[intersect(colnames(df.select.nor),rownames(endoderm.use)),]

data.expr <-read.delim("get_germlayer_tss100kb_k27ac_signal/merge.matrix.sort.TSS100KB.txt",header = T,row.names = 1,sep="\t")
data.expr<- t(t(data.expr)/colSums(data.expr))
dim(data.expr)


library(ggplot2)
color.table <- as.data.frame(dim)
data.expr <- data.expr[,match(rownames(dim),colnames(data.expr))]
feature <- selected.region

for (i in 1:length(feature))
  
{ 
  data.change <- data.expr
  color.table[, 3] <- as.numeric(t(data.change[feature[i],]))
  color.table.sort <- color.table[order(color.table[, 3]), ]
  color.table.sort[, 4] <- rank(color.table.sort[,3])
  
  colnames(color.table.sort)[3:4] <- c("Expression","Rank")
  
  pdf(paste0("./feature_gene_projection/", feature[i], ".projection.pdf"),width = 5.05,height=4)
  p <- ggplot(color.table.sort, aes(x=ESC,y = meso,colour = Expression))+
    geom_point(size=1.5,alpha=0.8) + 
    scale_colour_gradientn(colours = c("grey","yellow","red"))+
    labs(x="dim1", y="dim2")+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))#make the title center
  print(p)  
  dev.off()
}


library(ggplot2)
library(gplots)
library(gridExtra)
library(car)
library(vegan)
library(reshape)
library(reshape2)
library(broom)
library(pvclust)
library(dendextend)
library(dplyr)

rm(list=ls())

Fldata<-read.csv("MOTEAcroThermalExp_2021.csv", header=TRUE)   ######data for figure 3
Blmetrics<-read.csv("MOTEbleaching_response_metrics_2021.csv", header=TRUE)	#####qPCR data for figure 3
dim(Fldata)
dim(Blmetrics)
rownames(Blmetrics)<-Blmetrics[,1]
Blmetrics<-Blmetrics[,-1]
rownames(Fldata)<-Fldata[,1]
Fldata<-t(Fldata[,-1])

rownames(Blmetrics)
colnames(Blmetrics)
rownames(Fldata)

FinDa<-cbind(Fldata,Blmetrics)

col<-c(0,1,2,3,4)
stage<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47)
metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ")
metrix<-c(1,2,3,4,5,6,7,8)
h<-1
metaphys<-list()
for( I in 1:5)
{
	for (TT in 1:47)
	{
		for (FF in 1:8)
		{
			uni<-paste(metric[FF], stage[TT], col[I], sep=".")
			metaphys[[h]]<-data.frame("uni"=uni, "metric"=metric[FF], "metrix"=metrix[FF], "stage"=stage[TT], "col"=(col[I]+1))
			h<-h+1
		}
	}
}
metafluor<-data.frame(do.call(rbind, metaphys))

###############
###############
###############
###############

library(igraph)

cor_mat <- as.matrix(as.dist(cor(t(t(FinDa)), method = c("spearman"), use="complete.obs")))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = TRUE)
cor_edge_list <- as_data_frame(cor_g, 'edges')
dim(cor_edge_list)
summary(cor_edge_list)

summary(cor_mat)

possibles<-colnames(Blmetrics)
#possibles<-rhar
finset<-subset(cor_edge_list, abs(weight)>0.5)
dim(finset)
summary(finset)

finer<-list()
gg<-1
for (t in 1:dim(finset)[1])
{	if(finset[t,1] %in% possibles)
	{	if(!(finset[t,2] %in% possibles))
		{
			finer[[gg]]<-cbind(finset[t,], "lm"=1, "phys"=2)
			gg<-gg+1
		}
	}
	if(finset[t,2] %in% possibles)
	{	if(!(finset[t,1] %in% possibles))
		{
			finer[[gg]]<-cbind(finset[t,], "phys"=2, "lm"=1)
			gg<-gg+1
		}
	}
}
resulter<-data.frame(do.call(rbind, finer))
#resulter<-finset

dim(resulter)
res1<-data.frame("vert"=resulter$from, "type"=resulter$phys, "weight"=resulter$weight)
res2<-data.frame("vert"=resulter$to, "type"=resulter$lm, "weight"=(resulter$weight)^3)
res3<-rbind(res1,res2)
res4<-res3[!duplicated(res3$vert),]

ModelSet<- resulter[order(abs(resulter$weight)),]
dim(ModelSet)

wavlen<-data.frame("uni"=possibles, "metric"="absorb", "metrix"=9,"stage"=35,"col"=6) 
aller<-rbind(metafluor, wavlen)
alls<-aller[match(res4$vert, aller$uni),]
res5<-data.frame(res4, alls)

new_g<-graph_from_data_frame(resulter[,1:3],F)
g6<-simplify(new_g, remove.multiple = T, remove.loops = T)

V(g6)$type<-res5$type
V(g6)$clust<-res5$col
V(g6)$shape<-res5$metrix
E(g6)$color<-as.factor(sign(E(g6)$weight))
E(g6)$weight<-(((abs(E(g6)$weight)-range(abs(E(g6)$weight))[1])/(range(abs(E(g6)$weight))[2]-range(abs(E(g6)$weight))[1]))+1)^2




ll <- layout_with_dh(g6)
#ll <- layout_with_fr(g6)

samp<-paste("N-spearman", ".pdf", sep="")
pdf(file = samp, width = 10, height = 10, bg="transparent")
plot(g6, layout=ll, vertex.shape=c("circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "square")[V(g6)$shape], vertex.size=c(15, 8)[V(g6)$type], vertex.color=c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e","grey")[V(g6)$shape], vertex.frame.color="black", vertex.frame.cex=2.5, vertex.label.color="black", vertex.label.cex=c(0.0008, 0.0001)[V(g6)$type], vertex.label.dist=0.01, edge.width=0.26, edge.curved=0.2, edge.color=c("black", "black")[E(g6)$color])
dev.off()

############
############
############
############

gg<-list()
for(I in 1:nrow(ModelSet))
{
	x<-as.matrix(FinDa[ModelSet$to[I]])
	y<-as.matrix(FinDa[ModelSet$from[I]])
	cor.test(x, y, method=c("spearman"))
	#plot(x, y, pch=19, cex=2.2, col=c("black", "black"), xlim=c(-100,50))
	#abline(lm(y ~ x), col="red", lwd=3)
	zero<-lm(y ~ x)[1][[1]][1]
	slope<-lm(y ~ x)[1][[1]][2]
	sl<-c("positive")
	if(lm(y ~ x)[1][[1]][2]<0)
	{sl<-c("inverse")}
	#Sys.sleep(0.01)
	gg[[I]]<-data.frame("metricFrom"=ModelSet$from[I],"metricto"=ModelSet$to[I],"upper"=(30*slope)+zero, "lower"=(-30*slope)+zero, "slope"=sl, "weight"=abs(ModelSet$weight[I]))

}
lter<-data.frame(do.call(rbind, gg))
lter[1:10,]

#######
#######
#######
#######

FinRanges<-list()
TotalRange<-list()
rangeRes<-list()
for(I in 1:nrow(FinDa))
{   
	genotype<-rownames(FinDa)[I]
	for(K in 1:nrow(lter))
	{
		val<-FinDa[I,lter$metricFrom[K]]
		if(lter$slope[K]=="positive")
		{
			if(between(val, lter$lower[K], lter$upper[K]))
			#if(val>lter$lower[K])
			{
				#rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=lter$weight[K])
				rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=1)
			}
			else
			{
				rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=0)
			}
		}
		if(lter$slope[K]=="inverse")
		{
			if(between(val, lter$upper[K], lter$lower[K]))
			#if(val<lter$upper[K])
			{
				#rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=lter$weight[K])
				rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=1)
			}
			else
			{
				rangeRes[[K]]<-data.frame("metric"=lter$metricFrom[K], "InRange"=0)
			}
		}
		
	}
	ResRange<-data.frame(do.call(rbind, rangeRes))	
	FinRanges[[I]]<-data.frame("genotype"=genotype, ResRange)
	TotalRange[[I]]<-data.frame("genotype"=genotype, "summer"=sum(ResRange$InRange))
}
FinalRanges<-data.frame(do.call(rbind, FinRanges))
TotalRanges<-data.frame(do.call(rbind, TotalRange))
tot<-TotalRanges[order(TotalRanges$summer),]

DD<-data.frame(tot[match(rownames(Blmetrics), tot$genotype),], Blmetrics)

modelRes<-DD[order(DD$summer),]
modelRes


cor.test(modelRes$summer, modelRes$X679.455, method=c("spearman"))
summary(lm(modelRes$summer ~ modelRes$X679.455))

samp<-paste("ModelCorrelation", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")
par(mfrow=c(1,1), oma = c(3, 3, 0.1, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

x<-modelRes$summer
y<-modelRes$X679.455
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=19, cex=2.2, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)])
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()

############
############
############
############
############
############
############
############
############
############

samp<-paste("correlation-link-Quant.1.2-PerTau2", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")
par(mfrow=c(1,1), oma = c(3, 1, 1, 4), mar = c(0.1, 0.1, 0.1, 0.1))

x<-FinDa$PerTau2
y<-FinDa$Quant.1.2
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=2.25, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-50,125), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.75)
axis(4, at=, col.axis="black", las=2, cex.axis=1.75)
box(col="black")
abline(lm(y ~ x), col="red", lwd=5)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-X420-Connect.46.2", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

x<-FinDa$X420.585
y<-FinDa$Connect.46.2
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()

FinDa[,1:5]

samp<-paste("correlation-link-X679-NPQ.22.0", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

x<-FinDa$X679.455
y<-FinDa$NPQ.22.0
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-PerTau2-qP.12.3", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

x<-FinDa$PerTau2
y<-FinDa$qP.12.3
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-50,125), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-Quant.11.2-X420", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

x<-FinDa$X420.585
y<-FinDa$Quant.11.2
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()


samp<-paste("correlation-link-qP.1.4-X530", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

x<-FinDa$X530.457
y<-FinDa$qP.1.4
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-Connect.46.2-X530.457", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

y<-FinDa$Connect.46.2
x<-FinDa$X530.457
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-Tau2.47.0-X460.559", ".pdf", sep="")
pdf(file = samp, width = 5, height = 5, bg="transparent")

y<-FinDa$Tau2.47.0
x<-FinDa$X460.559
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,0), axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=2, cex.axis=1)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()



samp<-paste("correlation-link-NPQ-X505", ".pdf", sep="")
pdf(file = samp, width = 5, height = 3.5, bg="transparent")

par(mfrow=c(1,1), oma = c(3, 3, 1, 1), mar = c(0.1, 0.1, 0.1, 0.1))

x<-FinDa$X505.602
y<-FinDa$NPQ.22.0
cor.test(x, y, method=c("spearman"))
plot(x, y, pch=19, cex=1.75, col=c("black", "grey")[c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)], xlim=c(-100,30), ylim=c(0.2,0.9),axes=FALSE, ylab="", xlab="")
axis(1, at=, col.axis="black", las=1, cex.axis=1.25)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
box(col="black")
abline(lm(y ~ x), col="red", lwd=3)
summary(lm(y ~ x))

dev.off()

############
############
############
############
############
############
############
############

col<-c(0,1,2,3,4)
stage<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47)
metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ")
metrix<-c(1,2,3,4,5,6,7,8)
h<-1
metaphys<-list()
for( I in 1:5)
{
	for (TT in 1:47)
	{
		for (FF in 1:8)
		{
			uni<-paste(metric[FF], stage[TT], col[I], sep=".")
			metaphys[[h]]<-data.frame("uni"=uni, "metric"=metric[FF], "metrix"=metrix[FF], "stage"=stage[TT], "col"=(col[I]+1))
			h<-h+1
		}
	}
}
metafluor<-data.frame(do.call(rbind, metaphys))

newset<-subset(metafluor, stage=="1" | stage=="9" | stage=="10" | stage=="11" | stage=="47" | stage=="46" | stage=="45" | stage=="44"| stage=="20" | stage=="21" | stage=="22" | stage=="23"| stage=="30" | stage=="29" | stage=="40" | stage=="39" | stage=="38")
dim(newset)
modelReduced<-na.omit(lter[match(newset$uni, lter$metricFrom),])
dim(modelReduced)

modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.1\\.", ".2.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.9\\.", ".3.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.10\\.", ".4.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.11\\.", ".5.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.47\\.", ".34.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.46\\.", ".33.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.45\\.", ".32.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.44\\.", ".31.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.20\\.", ".11.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.21\\.", ".12.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.22\\.", ".13.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.23\\.", ".14.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.29\\.", ".19.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.30\\.", ".20.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.38\\.", ".27.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.39\\.", ".28.",modelReduced$metricFrom)
modelReduced$metricFrom
modelReduced$metricFrom<-gsub("\\.40\\.", ".29.",modelReduced$metricFrom)


dim(modelReduced)
ssss<-read.csv("MOTEfulldata_2021.csv", header=TRUE)   ######data for figure 3
rownames(ssss)<-ssss[,1]
ssss<-ssss[,-1]
dim(ssss)

colnames(ssss)<-sub('\\.','-', sub('\\.','-', colnames(ssss)))


FinRanges<-list()
TotalRange<-list()
rangeRes<-list()
for(I in 1:nrow(t(ssss)))
{   
	genotype<-rownames(t(ssss))[I]
	for(K in 1:nrow(modelReduced))
	{
		val<-t(ssss)[I,modelReduced$metricFrom[K]]
		if(modelReduced$slope[K]=="positive")
		{
			if(between(val, modelReduced$lower[K], modelReduced$upper[K]))
			#if(val>lter$lower[K])
			{
				rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=modelReduced$weight[K])
				#rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=1)
			}
			else
			{
				rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=0)
			}
		}
		if(modelReduced$slope[K]=="inverse")
		{
			if(between(val, modelReduced$upper[K], modelReduced$lower[K]))
			#if(val<lter$upper[K])
			{
				rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=modelReduced$weight[K])
				#rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=1)
			}
			else
			{
				rangeRes[[K]]<-data.frame("metric"=modelReduced$metricFrom[K], "InRange"=0)
			}
		}
		
	}
	ResRange<-data.frame(do.call(rbind, rangeRes))	
	FinRanges[[I]]<-data.frame("genotype"=genotype, ResRange)
	TotalRange[[I]]<-data.frame("genotype"=genotype, "summer"=sum(ResRange$InRange))
}
FinalRanges<-data.frame(do.call(rbind, FinRanges))
TotalRanges<-data.frame(do.call(rbind, TotalRange))
tot<-TotalRanges[order(TotalRanges$summer),]

metadata<-read.csv("MOTEmeta2021.csv", header=TRUE)	#####qPCR data for figure 3

dod<-data.frame(tot,metadata[,7:9][match(tot$genotype, metadata$uniq),])
dod


#par(mfrow=c(1,1), oma = c(4, 2, 1, 1), mar = c(1, 1, 1, 0.1))
# Get the stacked barplot
barplot(dod$summer, col=c("#a6611a","#dfc27d","#80cdc1","#018571")[as.factor(dod$Symb)], border="white", space=0.02, axes=F, las=2, cex.names=0.005)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)







samp<-paste("Predictive-MOTESupervisedBarFrac", ".pdf", sep="")
pdf(file = samp, width = 12, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(4, 2, 1, 1), mar = c(1, 1, 1, 0.1))
# Get the stacked barplot
barplot(dod$summer, col=c("#a6611a","#dfc27d","#80cdc1","#018571")[as.factor(dod$Symb)], border="white", space=0.02, axes=F, las=2, cex.names=0.005)
axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
dev.off()


samp<-paste("Predictive-MOTEcoralFrac", ".pdf", sep="")
pdf(file = samp, width = 12, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(4, 2, 1, 1), mar = c(1, 1, 1, 0.1))
# Get the stacked barplot
barplot(dod$summer, col=c("#d73027","#fc8d59","#fee090","#525252","#d9d9d9","#91bfdb","#4575b4")[as.factor(dod$stat)], border="white", space=0.02, axes=F, las=2, cex.names=0.5, ylim=c(0,0.1))
#axis(2, at=, col.axis="black", las=2, cex.axis=1.5)
dev.off()


shapiro.test(sqrt(dod$summer))
summary(aov(sqrt(dod$summer) ~ dod$Symb))
hyh<-aov(sqrt(dod$summer) ~ dod$Symb)
TukeyHSD(hyh, conf.level=.95)


############
############
############
############
############
############


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
library(colorspace)


rm(list=ls())

#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -

ss<-read.csv("MOTEfulldata_2021.csv", header=TRUE)   ######data for figure 3
data<-read.csv("MOTEqPCR_2021.csv", header=TRUE)	#####qPCR data for figure 3
metadata<-read.csv("MOTEmeta2021.csv", header=TRUE)	#####qPCR data for figure 3

colnames(ss)<-sub('\\.','-', sub('\\.','-', colnames(ss)))


rownames(ss)<-ss[,1]
ss<-ss[,-1]
stat<-ss
dim(stat)
rs<-rownames(stat)
cls<-colnames(stat)
Boo<-dim(stat)[[1]]
datum <- list()
replicates<-metadata$Symb[match(colnames(stat), metadata$uniq)]
pvale<-0.01
for (i in 1:Boo)
{	H<-rownames(stat)[i]	
	if(is.na(summary(aov(t(stat[i,])~replicates))[[1]][,5][1])==FALSE)
	{	if (as.matrix(shapiro.test(resid(aov(t(stat[i,])~replicates)))[[2]])>pvale)
		{	if(summary(aov(t(stat[i,])~replicates))[[1]][,5][1]<pvale)
			{	#print(summary(aov(t(stat[i,])~replicates))[[1]][,5][1])
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
		else
		{	if(kruskal.test(t(stat[i,])~replicates)[[3]]<pvale)
			{	#print(kruskal.test(t(stat[i,])~replicates)[[3]])
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
	}
}
allerie2<-data.matrix(do.call(rbind, datum))
colnames(allerie2)<-cls
allerie2<-na.omit(allerie2)
summary(allerie2)
dim(allerie2)
head(allerie2)

ss2<-t(allerie2)
ss3<-ss2[match(data$X, rownames(ss2)),]
allerie2<-t(na.omit(ss3))
dim(allerie2)

cb<-pvclust(allerie2, method.dist="euclidean", method.hclust="ward.D", nboot=10)
ca<-hclust(dist(allerie2, method="canberra"), method="ward.D2")
caclusters<-cutree(ca, k=6)     ######### Chose how many sections you want k=....

par(mfrow=c(1,1), oma = c(4, 2, 5, 5), mar = c(1, 1, 5, 0.1))
plot(cb)
cb$hclust$order

dd<-flip_leaves(as.dendrogram(cb$hclust), c(1), c(14))
plot(dd)
order.dendrogram(dd)
de<-flip_leaves(dd, c(75, 76, 67, 74, 30, 72, 73, 68, 69), c(15, 14, 1, 4, 6, 5, 2, 3))
plot(de)
order.dendrogram(de)
df<-flip_leaves(de, c(71, 66, 64, 65), c(15, 14, 1, 4, 6, 5, 2, 3, 75, 76, 67, 74, 30, 72, 73, 68, 69))
plot(df)
order.dendrogram(df)
dr<-flip_leaves(df, c(7, 9, 12), c(53, 52, 57, 19, 21, 10, 16, 8, 13, 11, 20, 17, 18))
plot(dr)
order.dendrogram(dr)
dh<-flip_leaves(dr, c(89, 35, 33, 34, 38, 45, 87, 40, 41, 82, 70, 90, 31, 32, 22, 23, 25, 24, 26, 29, 39, 36, 37, 80, 85, 42, 27, 28, 77, 78, 43, 84, 79, 46, 86, 48, 81, 50, 44, 83, 88, 47, 49), c(55, 62, 63, 51, 59, 58, 60, 61, 54, 56, 53, 52, 57, 19, 21, 10, 16, 8, 13, 11, 20, 17, 18, 7, 9, 12, 15, 14, 1, 4, 6, 5, 2, 3, 75, 76, 67, 74, 30, 72, 73, 68, 69, 71, 66, 64, 65))
plot(dh)
order.dendrogram(dh)
dc<-flip_leaves(dh, c(9), c(12,10))
plot(dc)
order.dendrogram(dc)
dx<-flip_leaves(dc, c(7), c(12,10,9))
plot(dx)
order.dendrogram(dx)
dxx<-flip_leaves(dx, c(80, 79, 46, 86, 48, 81, 50, 47, 49, 84, 85, 43, 77, 78, 38, 35, 33, 34, 82, 32, 31, 90, 40, 41, 45, 87, 70, 88, 44, 83, 89, 27, 28, 42, 22, 24, 26, 29, 23, 25, 39, 36, 37), c(51, 59, 55, 62, 63, 58, 60, 61, 54, 56, 53, 52, 57, 12, 10, 9, 7, 19, 21, 16, 8, 11, 13, 20, 17, 18, 15, 14, 1, 5, 2, 3, 4, 6, 75, 76, 67, 74, 30, 72, 73, 68, 69, 71, 66, 64, 65))
plot(dxx)
order.dendrogram(dxx)

gd<-as.dendrogram(dh)
ca <- color_branches(ca, k = 6, col = c("green","red", "blue", "gray", "orange", "black")) # add color to the lines
gd <- color_branches(gd, k = 3, col = c("#bdbdbd", "#737373", "black")) # add color to the lines


par(mfrow=c(3,1), oma = c(4, 2, 5, 5), mar = c(1, 1, 5, 0.1))

plot(gd, print.pv=TRUE, hang=0.1, cex=0.5)

colsplot<-labels(gd)
colsplot2<-data[match(labels(gd), data$X),]
data1<-colsplot2

# Get the stacked barplot
barplot(t(data1[,2:5]), col=c("#a6611a","#dfc27d","#80cdc1","#018571"), border="white", space=0.02, axes=F, las=2, cex.names=0.005)
axis(4, at=, col.axis="black", las=2, cex.axis=1.5)

colsplot<-labels(gd)
colsplot3<-metadata$Species[match(labels(gd), metadata$uniq)]
da1<-colsplot3

das1<-rowSums(data1[,2:5])
barplot(das1, col=c("#d73027","#fc8d59","#fee090","#525252","#d9d9d9","#91bfdb","#4575b4")[as.factor(da1)], border="white", space=0.02, axes=F, las=2, cex.names=0.5)
axis(4, at=, col.axis="black", las=2, cex.axis=1.5)

samp<-paste("HeatMap-MOTESupervised", ".pdf", sep="")
pdf(file = samp, width = 12, height = 8, bg="transparent")

col_breaks = c(seq(-3.0,-0.5,length.out=115), seq(-0.4,0.4,length.out=200), seq(0.5,3.0,length.out=115))
my_palette <- colorRampPalette(c("#2166ac","white","#d73027"))(length(col_breaks)-1)

par(oma=c(4,0.1,0.1,0.1))
heatmap.2(allerie2, col=my_palette, breaks=col_breaks, margins = c(2,3), density.info="none", trace="none", dendrogram = "both", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F, na.color = "pink", Colv=as.dendrogram(gd), Rowv=as.dendrogram(ca), labRow=ca$labels, RowSideColors=c("red","black", "orange", "gray", "green", "blue")[caclusters], scale="none", key=T, key.xlab="Z-score", key.title=NA, cexRow=0.00001, cexCol=0.001, srtCol=90, lmat=rbind(c(5,6,4), c(3,1,2)), lhei=c(2, 8), lwid=c(0.25,0.1, 3))

dev.off()


colsplot<-labels(gd)
colsplot2<-data[match(labels(gd), data$X),]
data1<-colsplot2

samp<-paste("HeatMap-MOTESupervisedBarFrac", ".pdf", sep="")
pdf(file = samp, width = 12, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(4, 2, 5, 5), mar = c(1, 1, 5, 0.1))
# Get the stacked barplot
barplot(t(data1[,2:5]), col=c("#a6611a","#dfc27d","#80cdc1","#018571"), border="white", space=0.02, axes=F, las=2, cex.names=0.005)
axis(4, at=, col.axis="black", las=2, cex.axis=1.5)

dev.off()

colsplot<-labels(gd)
colsplot3<-metadata$Species[match(labels(gd), metadata$uniq)]
da1<-colsplot3

samp<-paste("HeatMap-coralBarFrac", ".pdf", sep="")
pdf(file = samp, width = 12, height = 4, bg="transparent")

das1<-rowSums(data1[,2:5])
barplot(das1, col=c("#d73027","#fc8d59","#fee090","#525252","#d9d9d9","#91bfdb","#4575b4")[as.factor(da1)], border="white", space=0.02, axes=F, las=2, cex.names=0.5)
axis(4, at=, col.axis="black", las=2, cex.axis=1.5)

dev.off()

##########################
##########################
##########################
##########################
##########################

caclust<-data.frame(caclusters)
ca4<-subset(caclust, caclusters=="4")
ca1<-subset(caclust, caclusters=="1")
ca3<-subset(caclust, caclusters=="3")
ca6<-subset(caclust, caclusters=="6")
ca2<-subset(caclust, caclusters=="2")
ca5<-subset(caclust, caclusters=="5")

dim(ca1)
dim(ca2)
dim(ca3)
dim(ca4)
dim(ca5)
dim(ca6)

caa<-ca3
ca1Rich<-metafluor[match(rownames(caa), metafluor$uni),]
ca1Rich1<-data.frame(aggregate(ca1Rich[,2], by=list(ca1Rich$metric), FUN=NROW))
ca1m<-data.frame(ca1Rich1, "percent"=ca1Rich1$x/dim(ca1Rich)[1])
ca1Rich2<-data.frame(aggregate(ca1Rich[,5], by=list(ca1Rich$col), FUN=NROW))
ca1c<-data.frame(ca1Rich2, "percent"=ca1Rich2$x/dim(ca1Rich)[1])
ca1Rich3<-data.frame(aggregate(ca1Rich[,4], by=list(ca1Rich$stage), FUN=NROW))
ca1t<-data.frame(ca1Rich3, "percent"=ca1Rich3$x/dim(ca1Rich)[1])
ca1m
ca1c
ca1t

##########################
##########################
##########################
##########################
cbclusters<-cutree(dh, k=3)     ######### Chose how many sections you want k=....
cbclust<-data.frame(cbclusters)

cb1<-subset(cbclust, cbclusters=="2") #red
cb2<-subset(cbclust, cbclusters=="1") #grey
cb3<-subset(cbclust, cbclusters=="3") #blue

allerie4<-t(ssss)
dim(allerie4)

phen1<-allerie4[match(rownames(cb1),rownames(allerie4)),]
dim(phen1)

phen2<-allerie4[match(rownames(cb2),rownames(allerie4)),]
dim(phen2)

phen3<-allerie4[match(rownames(cb3),rownames(allerie4)),]
dim(phen3)

##########Stats###########
##########Stats###########
##########Stats###########

library(lme4)
library(lmerTest)
library(multcomp)

bbbb<-data.frame(metafluor, t(phen1), t(phen2), t(phen3))
dim(bbbb)

factor<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
fact<-c("P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P2", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3", "P3")

Bquant1<-subset(bbbb, metric=="Quant" & col=="1")
Bquant2<-subset(bbbb, metric=="Quant" & col=="2")
Bquant3<-subset(bbbb, metric=="Quant" & col=="3")
Bquant4<-subset(bbbb, metric=="Quant" & col=="4")
Bquant5<-subset(bbbb, metric=="Quant" & col=="5")

quant1<-Bquant5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan), quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))




Bsigma1<-subset(bbbb, metric=="Sigma" & col=="1")
Bsigma2<-subset(bbbb, metric=="Sigma" & col=="2")
Bsigma3<-subset(bbbb, metric=="Sigma" & col=="3")
Bsigma4<-subset(bbbb, metric=="Sigma" & col=="4")
Bsigma5<-subset(bbbb, metric=="Sigma" & col=="5")

quant1<-Bsigma5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))




B1tau1<-subset(bbbb, metric=="Tau1" & col=="1")
B2tau1<-subset(bbbb, metric=="Tau1" & col=="2")
B3tau1<-subset(bbbb, metric=="Tau1" & col=="3")
B4tau1<-subset(bbbb, metric=="Tau1" & col=="4")
B5tau1<-subset(bbbb, metric=="Tau1" & col=="5")

quant1<-B5tau1[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))




B1tau2<-subset(bbbb, metric=="Tau2" & col=="1")
B2tau2<-subset(bbbb, metric=="Tau2" & col=="2")
B3tau2<-subset(bbbb, metric=="Tau2" & col=="3")
B4tau2<-subset(bbbb, metric=="Tau2" & col=="4")
B5tau2<-subset(bbbb, metric=="Tau2" & col=="5")

quant1<-B5tau2[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))

##########Stats###########
##########Stats###########
##########Stats###########

P1<-phen1
P2<-phen2
P3<-phen3


stupid<-list()
for(I in 1:ncol(P1))
{
	stupid[[I]]<-data.frame("phys"=colnames(P1)[I], "mean"=mean(as.numeric(P1[,I])), "sd"=sd(as.numeric(P1[,I]))/sqrt(nrow(P1)))
	#stupid[[I]]<-data.frame("phys"=colnames(P1)[I], "mean"=mean(as.numeric(P1[,I])), "sd"=sd(as.numeric(P1[,I])))
}
allerieP1<-data.frame(do.call(rbind, stupid))

stupid<-list()
for(I in 1:ncol(P2))
{
	stupid[[I]]<-data.frame("phys"=colnames(P2)[I], "mean"=mean(as.numeric(P2[,I])), "sd"=sd(as.numeric(P2[,I]))/sqrt(nrow(P2)))
	#stupid[[I]]<-data.frame("phys"=colnames(P2)[I], "mean"=mean(as.numeric(P2[,I])), "sd"=sd(as.numeric(P2[,I])))
}
allerieP2<-data.frame(do.call(rbind, stupid))

stupid<-list()
for(I in 1:ncol(P3))
{
	stupid[[I]]<-data.frame("phys"=colnames(P3)[I], "mean"=mean(as.numeric(P3[,I])), "sd"=sd(as.numeric(P3[,I]))/sqrt(nrow(P3)))
	#stupid[[I]]<-data.frame("phys"=colnames(P3)[I], "mean"=mean(as.numeric(P3[,I])), "sd"=sd(as.numeric(P3[,I])))

}
allerieP3<-data.frame(do.call(rbind, stupid))


pheno1<-data.frame(metafluor, allerieP1[,2:3])
pheno2<-data.frame(metafluor, allerieP2[,2:3])
pheno3<-data.frame(metafluor, allerieP3[,2:3])

metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ")
sets<-c(0.5,8,1,6000,4500000,0.75,1,0.5)
TT<-c(1,2,4,5)


samp<-paste("Profiles-BCD", ".pdf", sep="")
pdf(file = samp, width = 9.0, height = 7, bg="transparent")

par(mfrow=c(4,3), oma = c(4, 5.9, 0.1, 0.5), mar = c(1, 1, 0.1, 0.25))

for (G in 1:4)
{
	X<-TT[G]
	print(metric[X])
	ggg<-list(pheno1,pheno2,pheno3)
	for (I in 1:3)
	{	col1<-subset(ggg[[I]], metric==metric[X] & col=="1")
		col2<-subset(ggg[[I]], metric==metric[X] & col=="2")
		col3<-subset(ggg[[I]], metric==metric[X] & col=="3")
		col4<-subset(ggg[[I]], metric==metric[X] & col=="4")
		col5<-subset(ggg[[I]], metric==metric[X] & col=="5")
		
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col1$stage, col1$mean+col1$sd, type="l", lwd=0.5, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col1$stage, col1$mean-col1$sd, type="l", lwd=0.5, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean+col2$sd, type="l", lwd=0.5, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean-col2$sd, type="l", lwd=0.5, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean+col3$sd, type="l", lwd=0.5, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean-col3$sd, type="l", lwd=0.5, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean+col4$sd, type="l", lwd=0.5, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean-col4$sd, type="l", lwd=0.5, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean, type="l", lwd=1.25, col="green", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean+col5$sd, type="l", lwd=0.5, col="green", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean-col5$sd, type="l", lwd=0.5, col="green", ylim=c(0,sets[X]), axes=F)
		if(I==1)
		{	axis(2, at=, col.axis="black", las=2, cex.axis=2)}
		if(G==4)
		{	axis(1, at=, col.axis="black", las=1, cex.axis=2)}
		box(col="black")
		#grid(col="gray")
	}
}
dev.off()


##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########








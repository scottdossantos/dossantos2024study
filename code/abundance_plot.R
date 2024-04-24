m <- read.table("data/merged_readcounts_taxonomy_qiime_L6_16Sorder_color_sum.txt", header=T, row.names=1, sep="\t", comment.char="", check.names=F)

r <- read.table("data/td_OTU_tag_mapped_RDPlineage_blastcorrected_lactospecies.txt", header=T, row.names=1, sep="\t", comment.char="", check.names=F, skip=1)

r.taxon <- r$taxonomy
r.tax <- vector()
for(i in 1:length(r.taxon)){
 r.tax[i] <- as.character(strsplit(as.vector(r.taxon[i]), ";")[[1]][6])
}
r$taxonomy <- NULL

m.tax <- vector()
m$total <- NULL
m$color <- NULL

for(i in 1:nrow(m)){
 m.tax[i] <- as.character(strsplit(as.vector(rownames(m)[i]), ";")[[1]][6])
}

agg.m <- aggregate(m, by=list(m.tax), FUN=sum)
rownames(agg.m) <- agg.m$Group.1
agg.m$Group.1 <- NULL

agg.r <- aggregate(r, by=list(r.tax), FUN=sum)
rownames(agg.r) <- agg.r$Group.1
agg.r$Group.1 <- NULL

species <- c( "Lactobacillus_iners", "Lactobacillus_crispatus", "Megasphaera",  "Gardnerella", "Atopobium",  "Prevotella", "Sneathia", "Bifidobacterium" , "Dialister", "Streptococcus"   )
species.names <- c( "L_iners", "L_crispatus", "Megasphaera",  "Gardnerella", "Atopobium",  "Prevotella", "Sneathia", "Bifidobacterium" , "Dialister", "Streptococcus"   )

color=c("steelblue3", "skyblue1", "olivedrab3", "indianred1", "ivory2", "mediumpurple1", "pink", "mediumvioletred", "tan1", "#CC9933")

agg.m[agg.m<1] <- 0.01
agg.r.m[agg.r.m<1] <- 0.9
clr.m <- apply(agg.m, 2, function(x){log(x) - mean(log(x))})
clr.r <- apply(agg.r.m, 2, function(x){log(x) - mean(log(x))})

# http://www.surefoss.org/visualization%20/%20visualisierungen/visualizing-small-scale-paired-data-combining-boxplots-stripcharts-and-confidence-intervals-in-r/
pdf("figs/abundance.pdf", height=6, width=8)
par(mfrow=c(2,5))

for(i in 1:length(species)){
rRNA <- as.vector(as.matrix(clr.r[species[i],]))
mRNA <- as.vector(as.matrix(clr.m[species[i],]))
#Settinguptwoscreens
#par(mfrow=c(1,1))

#FirstGraph
s<-seq(length(rRNA))
#par(bty="l")
boxplot(rRNA,mRNA,main=species.names[i],xlab=NULL,ylab="Measure",names=c("16S","mRNA"),col=color[i], ylim=c(-4,10))
stripchart(list(rRNA,mRNA),vertical=T,pch=16,method="jitter",cex=0.5,add=T)

segments(rep(1.05,length(rRNA))[s],rRNA[s],rep(1.95,length(rRNA))[s],mRNA[s],col="grey",lwd=1, lty=2)

}
dev.off()

pdf("figs/ugly_16S_mRNA.pdf")
plot(as.matrix(clr.r[species[1],]), as.matrix(clr.m[species[1],1:24]), col=color[1], ylim=c(-4,9), xlim=c(-4,9), xlab="16S rRNA", ylab="mRNA", pch=19, cex=0.5)
lines(lowess(as.matrix(clr.r[species[1],]),as.matrix(clr.m[species[1],1:24])) , col=color[1])

for(i in 2:length(species)){
	points(as.matrix(clr.r[species[i],]), as.matrix(clr.m[species[i],1:24]), col=color[i], pch=19, cex=0.5)
	lines(lowess(as.matrix(clr.r[species[i],]),as.matrix(clr.m[species[i],1:24])), col=color[i])
}
abline(0,1, lty=2, col="grey")
legend(-4,9, species, fill=color, cex=0.7)
dev.off()

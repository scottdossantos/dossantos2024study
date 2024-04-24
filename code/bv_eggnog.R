# run code/setup.R
# run code/merge_table.R
# run code/bv_biplot.R
# eggnot biplot only
# uses data from the BV biplot dataset
#
#####
# NOTE: this was run once, figures loaded into .Rmd from fig directory
# data saved into Rdata/
# ALDEx2 stuff is kinda slow (minutes to less than an hour, but still ...
#####

# Very ugly code - leave for Greg


coglist <- read.table("~/Desktop/coglist.txt", header=T, stringsAsFactors=F, sep="\t",  quote="", comment.char="", row.names=2, fill=T)

# aggregate by eggNOG classification
egg.vec <- egg[rownames(bv.filt),6]
names(egg.vec) <- rownames(bv.filt)
egg.bv <- aggregate(bv.filt, by=list(egg.vec), FUN=sum)
rownames(egg.bv) <- egg.bv$Group.1
egg.bv$Group.1 <- NULL

# make a compositional biplot based on clr
egg.clr <- apply(egg.bv+0.5, 2, function(x) log(x) - mean(log(x)))
egg.pcx <- prcomp(t(egg.clr))
# to do make this a biplot with codaSeq and cog groups
# biplot(egg.pcx, var.axes=F, cex=c(1, 0.5), col=c("red", rgb(0,0,0,0.4)))


## what differs between BV from two sites
conds <- c(rep("L", 14), rep("E", 14))
#x <- aldex.clr(egg.bv, conds, denom="lvha")
#x.e <- aldex.effect(x)
#x.t <- aldex.ttest(x)
#x.all <- data.frame(x.e, x.t)
#save(x.all, file="Rdata/x.all.Rda")
load("Rdata/x.all.Rda")


### what differs between two types of BV
bv.conds <- rep("A", 28)
bv.conds[which(bv.pcx$x[,2] > 0)] <- "B"
#x.bv <- aldex.clr(egg.bv, bv.conds, denom="lvha")
#x.bv.e <- aldex.effect(x.bv)
#x.bv.t <- aldex.ttest(x.bv)
#x.bv.all <- data.frame(x.bv.e, x.bv.t)
#save(x.bv.all, file='Rdata/x.bv.all.Rda')
load('Rdata/x.bv.all.Rda')

### once more with a glm
### use denom that is the intersect of the two above denoms
### almost co-incedent anyway
#covariates <- data.frame("A" = conds,"B" = bv.conds)
#mm <- model.matrix(~ B + A, covariates)
#denom <- intersect(x.bv@denom, x@denom)
#x.g <- aldex.clr(egg.bv, mm, mc.samples=128, denom=denom)
#x.g.glm <- aldex.glm(x.g)
#save(x.g.glm, file="Rdata/x.g.glm.Rda")
load('Rdata/x.g.glm.Rda')


##### selbal
##### > 30 hours, 10 cores, no result
#library(selbal)
#
#conds <- c(rep("L", 12), rep("E", 14))
#x <- t(egg.bv)
#y <- factor(conds)
#
#bal.dic <- selbal.cv(x=x, y=y, n.fold=5, n.iter=10, logit.acc='AUC')



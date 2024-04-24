### this contains plotting code

plot(x.bv.all$diff.win, x.bv.all$diff.btw, pch=19, col=rgb(0,0,0,0.3))
points(x.bv.all$diff.win[x.g.glm[,14] < 0.1], x.bv.all$diff.btw[x.g.glm[,14] < 0.1], pch=19, col=rgb(1,0,0,0.7))
points(x.bv.all$diff.win[x.bv.all$effect > 2], x.bv.all$diff.btw[x.bv.all$effect > 2], pch=19, col=rgb(0,0,1,0.7))


plot(x.all$diff.win, x.all$diff.btw, pch=19, col=rgb(0,0,0,0.3))
points(x.all$diff.win[x.g.glm[,15] < 0.1], x.all$diff.btw[x.g.glm[,15] < 0.1], pch=19, col=rgb(1,0,0,0.7))
points(x.all$diff.win[x.all$effect > 2], x.all$diff.btw[x.all$effect > 2], pch=19, col=rgb(0,0,1,0.7))

coglist[rownames(x.bv.all)[x.bv.all$effect > 1.5],"Enz"]
coglist[rownames(x.g.glm)[x.g.glm[,14] < 0.1],"Enz"]


coglist[rownames(x.all)[x.all$effect > 2],"Enz"]

egg.lookup <- match(rownames(x.all), egg$V7)

egg.COG <- data.frame(egg$V5[egg.lookup])

rownames(egg.COG) <- rownames(x.all)

colnames(egg.COG) <- "group"

codaSeq.stripchart(aldex.out=x.all, group.table=egg.COG, group.label="group", sig.method="effect", x.axis="effect", sig.cutoff=2)



# get a list of enzymes that differ between BV types
coglist[rownames(x.bv.all)[x.g.glm[,14] < 0.1 & x.bv.all$effect < 1],"Enz"]
codaSeq.stripchart(aldex.out=x.bv.all, group.table=egg.COG, group.label="group", sig.method="we.eBH", x.axis="effect")




bv.ko.data <- data.frame(ko.both[,grep('v.', colnames(ko.both))])

bv.ko.filt <-  codaSeq.filter(bv.ko.data, min.occurrence=0.30,
    min.prop=0.00005, samples.by.row=F)
# need to use zCompositions (cmultRepl function)
# rather than add 0.5
bv.ko.clr <- apply(bv.ko.filt+0.5, 2, function(x) log(x) - mean(log(x)))

bv.ko.pcx <- prcomp(t(bv.ko.clr))
biplot(bv.ko.pcx, cex=c(0.7,0.1), col=c('black', rgb(1,0,0,0.5)), var.axes=F, scale=0)

# get % var explained
bv.ko.pcx$sd^2/sum(bv.ko.pcx$sd^2)

# ALDEx2 (dev version) 
conds <- c(rep('L',14), rep('E', 14))
devtools::load_all('~/Documents/0_git/ALDEx_bioc')

x <- aldex.clr(bv.ko.filt, conds)
x.e <- aldex.effect(x)
aldex.plot(x.e, test='effect')

rownames(x.e)[x.e$effect < -1]

# https://www.genome.jp/kegg-bin/search_pathway_object
# flagella and motility and quorum as reference
#K00575 K00798 K01698 K01932 K02390 K02392 K02396 K02397
#K02398 K02400 K02409 K02410 K02412 K02414 K02416 K02417
#K02422 K02556 K03406 K03407 K03408 K03410 K03411 K03412
#K03415 K06595 K06603 K09749 K09770 K11936 K13244 K14051



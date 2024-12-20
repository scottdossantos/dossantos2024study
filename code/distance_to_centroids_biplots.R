# calculating distance of each pathway to species centroids

############# steup
# load packages 
library(dplyr)
library(zCompositions)
library(CoDaSeq)
library(tidyr)
library(ggplot2)

# set path to github directory
user <- "sds/cc"
if(user=="sds/cc"){
  path.to.github<-"~/Documents/GitHub/dossantos2024study/"
}else if(user=="gg"){
  path.to.github<-"~/Documents/0_git/projects/dossantos2024study/"
}

# load feature table from github for combined dataset, change name and clean-up
load(paste(path.to.github,"Rdata/both.data.Rda",sep = ""))
data.all.both<-both.data
rm(both.data)

# split dataframe into bv only, europe bv only and europe bv only subsets
data.bv.london<-data.all.both[,c(9:22)]
data.bv.europe<-data.all.both[,c(23:36)]
data.bv.both<-cbind(data.all.both[,c(9:22)],data.all.both[,c(23:36)])

# filter data subsets for genes present in >30% samples comprising >0.005% RA
filt.bv.london<-codaSeq.filter(data.bv.london,
                               min.occurrence = 0.30,
                               min.prop = 0.00005,
                               samples.by.row = FALSE)

filt.bv.europe<-codaSeq.filter(data.bv.europe,
                               min.occurrence = 0.30,
                               min.prop = 0.00005,
                               samples.by.row = FALSE)

filt.bv.both<-codaSeq.filter(data.bv.both,
                             min.occurrence = 0.30,
                             min.prop = 0.00005,
                             samples.by.row = FALSE)

# load pre-computed taxa table and colours
tax.table <- read.table(paste(path.to.github,'1_VIRGO/1.taxon.tbl.txt', sep=""),
                        header=F, row.names=2)

tax.colors <- read.table(paste(path.to.github,"code/species_colors.txt", sep=""), sep="\t",
                         header=T, row.names=1, stringsAsFactors=F)

# load pre-computed pathway table and colours
pathway.table <- read.table(paste(path.to.github,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""),
                            sep="\t",header=T, row.names=1, fill=TRUE)

pathway.colors <- read.table(paste(path.to.github,"code/pathway_colors_table3.txt",sep = ""),
                             sep="\t", header=F, row.names = 2 ,stringsAsFactors=F)


############# biplot rotations

# zero replacement using count multiplicative replacement from zCompositions package
# note: input transposed due to calculation of abundance across rows
london.bv.clr.input<-cmultRepl(t(filt.bv.london), method = "CZM",label = 0)
europe.bv.clr.input<-cmultRepl(t(filt.bv.europe), method = "CZM",label = 0)
both.bv.clr.input<-cmultRepl(t(filt.bv.both), method = "CZM",label = 0)

# centre log-ratio transformation of zero-replaced data
# note: input transposed to reverse transposition introduced above
# note2: log(x/y) is equal to log(x) - log(y)
london.bv.clr<- as.data.frame(apply(t(london.bv.clr.input), 2, function(x) log(x) - mean(log(x))))
europe.bv.clr<- as.data.frame(apply(t(europe.bv.clr.input), 2, function(x) log(x) - mean(log(x))))
both.bv.clr<- as.data.frame(apply(t(both.bv.clr.input), 2, function(x) log(x) - mean(log(x))))

# perform pca on clr-transformed data
london.bv.pca<-prcomp(t(london.bv.clr))
europe.bv.pca<-prcomp(t(europe.bv.clr))
both.bv.pca<-prcomp(t(both.bv.clr))

# get rotation matrices from pca objects
london.bv.rotation<-as.data.frame(london.bv.pca$rotation[,c(1:3)])
europe.bv.rotation<-as.data.frame(europe.bv.pca$rotation[,c(1:3)])
both.bv.rotation<-as.data.frame(both.bv.pca$rotation[,c(1:3)])

############# london centroids

# make a vector of taxa and pathways for london dataset
london.bv.tax.vec <- tax.table[rownames(filt.bv.london),2]
london.bv.path.vec <- pathway.table[rownames(filt.bv.london),2]

# assign V numbers to vector elements
names(london.bv.tax.vec) <- rownames(filt.bv.london)
names(london.bv.path.vec) <- rownames(filt.bv.london)

# get number of genes corresponding to each taxon and pathway 
london.bv.no.genes.tax<-vector()
for (i in levels(factor(london.bv.tax.vec))) {
  london.bv.no.genes.tax[i]<-length(which(london.bv.tax.vec==i))
}

london.bv.no.genes.path<-vector()
for (j in levels(factor(london.bv.path.vec))) {
  london.bv.no.genes.path[j]<-length(which(london.bv.path.vec==j))
}

# pull all species represented by >100 genes and top 20 pathways (desc. order)
london.bv.no.genes.100<-names(which(london.bv.no.genes.tax >100))
london.bv.no.genes.20<-names(head(sort(london.bv.no.genes.path,decreasing = TRUE),20))

# create a list of named vectors containing indices of all v numbers corresponding
# to each taxa represented by >100 genes
london.bv.taxa.indices<-list()
for (spp in 1:length(london.bv.no.genes.100)) {
  london.bv.taxa.indices[[london.bv.no.genes.100[spp]]]<-which(london.bv.tax.vec==london.bv.no.genes.100[spp])
}

# create a list of named vectors containing indices of all v numbers corresponding
# to each of the top 20 pathways
london.bv.pathway.indices<-list()
for (path in 1:length(london.bv.no.genes.20)) {
  london.bv.pathway.indices[[london.bv.no.genes.20[path]]]<-which(london.bv.path.vec==london.bv.no.genes.20[path])
}

# finds the V number indices present in the london bv feature table that are not
# present in the list of species represented by >100 genes (i.e. all other genes)
london.bv.taxa.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.london)),
                                           unlist(london.bv.taxa.indices))

# finds the V number indices present in the london bv feature table that are not
# present in the list of top 20 pathways (i.e. all other pathways)
london.bv.pathway.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.london)),
                                              unlist(london.bv.pathway.indices))

# for each species, get the indices of all the V numbers with pathway = "Ribosome"
# in the london.bv.taxa vector (i.e. genes we want to use to calculate 'centroid')
london.bv.centroid.indices<-list()
for (tax in 1:(length(names(london.bv.taxa.indices))-1)) {
  vec.indices<-which(london.bv.taxa.indices[[tax]] %in% london.bv.pathway.indices[["Ribosome"]])
  vec.toAdd<-vector()
  for (index in 1:length(vec.indices)) {
    vec.toAdd<-append(vec.toAdd,london.bv.taxa.indices[[tax]][[vec.indices[index]]],after = length(vec.toAdd))
  }
  london.bv.centroid.indices[[names(london.bv.taxa.indices)[[tax]]]]<-vec.toAdd
}

# make a list of dataframes containing the rotation values for the corresponding
# species centroids (i.e. ribosomal genes for each species)
london.bv.taxa.centroids<-list()
for(tax in 1:length(names(london.bv.centroid.indices))){
  london.bv.taxa.centroids[[names(london.bv.centroid.indices)[[tax]]]]<-
    as.data.frame(london.bv.rotation[london.bv.centroid.indices[[tax]],])
}

# make a blank dataframe containing median PC1/PC2/PC3 values for each taxon
# (i.e. taxon centroids on the PCA plot)
london.bv.centroid.summary<-as.data.frame(cbind(med_PC1=c(rep(NA,length(london.bv.taxa.centroids))),
                                                med_PC2=c(rep(NA,length(london.bv.taxa.centroids))),
                                                med_PC3=c(rep(NA,length(london.bv.taxa.centroids)))),
                                          row.names=names(london.bv.centroid.indices))

# fill in summary dataframe with values from the list of rotation dataframes
# for each taxon
for (tax in 1:length(rownames(london.bv.centroid.summary))) {
  df.tax<-as.data.frame(london.bv.taxa.centroids[[tax]])
  med.PC1<-median(london.bv.taxa.centroids[[tax]]$PC1)
  med.PC2<-median(london.bv.taxa.centroids[[tax]]$PC2)
  med.PC3<-median(london.bv.taxa.centroids[[tax]]$PC3)
  london.bv.centroid.summary[tax,1]<-med.PC1
  london.bv.centroid.summary[tax,2]<-med.PC2
  london.bv.centroid.summary[tax,3]<-med.PC3
}

# make list of dataframes containing rotations for all V numbers in each pathway
london.bv.pathway.rotations<-list()
for(tax in 2:20){
  london.bv.pathway.rotations[[names(london.bv.pathway.indices)[[tax]]]]<-
    as.data.frame(london.bv.rotation[london.bv.pathway.indices[[tax]],])
}

# make list of dataframes (1 per taxon) containing values denoting the 
# difference between PC1 rotation of each V number and PC1 taxon centroid (i.e.
# median PC1 value for that taxon). Note: values raised to power of 2 and sq.
# rooted to get around +/- values
london.bv.cent.dif.PC1<-list()
for (t in 1:length(names(london.bv.pathway.rotations))) {
  df.rot<-as.data.frame(london.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(london.bv.centroid.summary$med_PC1)))
  for (u in 1:(length(london.bv.centroid.summary$med_PC1))) {
    minus<-london.bv.centroid.summary$med_PC1[u]
    minus.vec<-df.rot$PC1
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(london.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  london.bv.cent.dif.PC1[[names(london.bv.pathway.rotations)[t]]]<-df.new
}

# collapse list into large dataframe with new column to identify list element
# origin (i.e. which pathway each v number corresponds to)
london.bv.cent.dif.PC1.all<-bind_rows(london.bv.cent.dif.PC1,.id = "pathway")

# transform data into 'tidy' format for ggplot faceting & make boxplots of
# difference between PC1 taxon centroid and PC1 pathway centroids, for each 
# taxon
london.bv.plot.PC1<-london.bv.cent.dif.PC1.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC1 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(london.bv.cent.dif.PC1.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.035),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("london_bv_PC1_difference.png",units = "in",height = 6,width = 15,res = 400)
london.bv.plot.PC1
# dev.off()

# same for PC 2 
london.bv.cent.dif.PC2<-list()
for (t in 1:length(names(london.bv.pathway.rotations))) {
  df.rot<-as.data.frame(london.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(london.bv.centroid.summary$med_PC2)))
  for (u in 1:(length(london.bv.centroid.summary$med_PC2))) {
    minus<-london.bv.centroid.summary$med_PC2[u]
    minus.vec<-df.rot$PC2
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(london.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  london.bv.cent.dif.PC2[[names(london.bv.pathway.rotations)[t]]]<-df.new
}

london.bv.cent.dif.PC2.all<-bind_rows(london.bv.cent.dif.PC2,.id = "pathway")

london.bv.plot.PC2<-london.bv.cent.dif.PC2.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC2 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(london.bv.cent.dif.PC2.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.035),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("london_bv_PC2_difference.png",units = "in",height = 6,width = 15,res = 400)
london.bv.plot.PC2
# dev.off()

# same for PC3
london.bv.cent.dif.PC3<-list()
for (t in 1:length(names(london.bv.pathway.rotations))) {
  df.rot<-as.data.frame(london.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(london.bv.centroid.summary$med_PC3)))
  for (u in 1:(length(london.bv.centroid.summary$med_PC3))) {
    minus<-london.bv.centroid.summary$med_PC3[u]
    minus.vec<-df.rot$PC3
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(london.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  london.bv.cent.dif.PC3[[names(london.bv.pathway.rotations)[t]]]<-df.new
}

london.bv.cent.dif.PC3.all<-bind_rows(london.bv.cent.dif.PC3,.id = "pathway")

london.bv.plot.PC3<-london.bv.cent.dif.PC3.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC3 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(london.bv.cent.dif.PC3.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.035),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("london_bv_PC3_difference.png",units = "in",height = 6,width = 15,res = 400)
london.bv.plot.PC3
# dev.off()

############# europe centroids

# make a vector of taxa and pathways for europe dataset
europe.bv.tax.vec <- tax.table[rownames(filt.bv.europe),2]
europe.bv.path.vec <- pathway.table[rownames(filt.bv.europe),2]

# assign V numbers to vector elements
names(europe.bv.tax.vec) <- rownames(filt.bv.europe)
names(europe.bv.path.vec) <- rownames(filt.bv.europe)

# get number of genes corresponding to each taxon and pathway
europe.bv.no.genes.tax<-vector()
for (i in levels(factor(europe.bv.tax.vec))) {
  europe.bv.no.genes.tax[i]<-length(which(europe.bv.tax.vec==i))
}

europe.bv.no.genes.path<-vector()
for (j in levels(factor(europe.bv.path.vec))) {
  europe.bv.no.genes.path[j]<-length(which(europe.bv.path.vec==j))
}

# pull all species represented by >100 genes and top 20 pathways (desc. order)
europe.bv.no.genes.100<-names(which(europe.bv.no.genes.tax >100))
europe.bv.no.genes.20<-names(head(sort(europe.bv.no.genes.path,decreasing = TRUE),20))

# create a list of named vectors containing indices of all v numbers corresponding
# to each taxa represented by >100 genes
europe.bv.taxa.indices<-list()
for (spp in 1:length(europe.bv.no.genes.100)) {
  europe.bv.taxa.indices[[europe.bv.no.genes.100[spp]]]<-which(europe.bv.tax.vec==europe.bv.no.genes.100[spp])
}

# create a list of named vectors containing indices of all v numbers corresponding
# to each of the top 20 pathways
europe.bv.pathway.indices<-list()
for (path in 1:length(europe.bv.no.genes.20)) {
  europe.bv.pathway.indices[[europe.bv.no.genes.20[path]]]<-which(europe.bv.path.vec==europe.bv.no.genes.20[path])
}

# finds the V number indices present in the europe bv feature table that are not
# present in the list of species represented by >100 genes (i.e. all other genes)
europe.bv.taxa.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.europe)),
                                           unlist(europe.bv.taxa.indices))

# finds the V number indices present in the europe bv feature table that are not
# present in the list of top 20 pathways (i.e. all other pathways)
europe.bv.pathway.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.europe)),
                                              unlist(europe.bv.pathway.indices))

# for each species, get the indices of all the V numbers with pathway = "Ribosome"
# in the europe.bv.taxa vector (i.e. genes we want to use to calculate 'centroid')
europe.bv.centroid.indices<-list()
for (tax in 1:(length(names(europe.bv.taxa.indices))-1)) {
  vec.indices<-which(europe.bv.taxa.indices[[tax]] %in% europe.bv.pathway.indices[["Ribosome"]])
  vec.toAdd<-vector()
  for (index in 1:length(vec.indices)) {
    vec.toAdd<-append(vec.toAdd,europe.bv.taxa.indices[[tax]][[vec.indices[index]]],after = length(vec.toAdd))
  }
  europe.bv.centroid.indices[[names(europe.bv.taxa.indices)[[tax]]]]<-vec.toAdd
}

# make a list of dataframes containing the rotation values for the corresponding
# species centroids (i.e. ribosomal genes for each species)
europe.bv.taxa.centroids<-list()
for(tax in 1:length(names(europe.bv.centroid.indices))){
  europe.bv.taxa.centroids[[names(europe.bv.centroid.indices)[[tax]]]]<-
    as.data.frame(europe.bv.rotation[europe.bv.centroid.indices[[tax]],])
}

# make a blank dataframe containing median PC1/PC2/PC3 values for each taxon
# (i.e. taxon centroids on the PCA plot)
europe.bv.centroid.summary<-as.data.frame(cbind(med_PC1=c(rep(NA,length(europe.bv.taxa.centroids))),
                                                med_PC2=c(rep(NA,length(europe.bv.taxa.centroids))),
                                                med_PC3=c(rep(NA,length(europe.bv.taxa.centroids)))),
                                          row.names=names(europe.bv.centroid.indices))

# fill in summary dataframe with values from the list of rotation dataframes
# for each taxon
for (tax in 1:length(rownames(europe.bv.centroid.summary))) {
  df.tax<-as.data.frame(europe.bv.taxa.centroids[[tax]])
  med.PC1<-median(europe.bv.taxa.centroids[[tax]]$PC1)
  med.PC2<-median(europe.bv.taxa.centroids[[tax]]$PC2)
  med.PC3<-median(europe.bv.taxa.centroids[[tax]]$PC3)
  europe.bv.centroid.summary[tax,1]<-med.PC1
  europe.bv.centroid.summary[tax,2]<-med.PC2
  europe.bv.centroid.summary[tax,3]<-med.PC3
}

# make list of dataframes containing rotations for all V numbers in each pathway
europe.bv.pathway.rotations<-list()
for(tax in 2:20){
  europe.bv.pathway.rotations[[names(europe.bv.pathway.indices)[[tax]]]]<-
    as.data.frame(europe.bv.rotation[europe.bv.pathway.indices[[tax]],])
}

# make list of dataframes (1 per taxon) containing values denoting the 
# difference between PC1 rotation of each V number and PC1 taxon centroid (i.e.
# median PC1 value for that taxon). Note: values raised to power of 2 and sq.
# rooted to get around +/- values
europe.bv.cent.dif.PC1<-list()
for (t in 1:length(names(europe.bv.pathway.rotations))) {
  df.rot<-as.data.frame(europe.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(europe.bv.centroid.summary$med_PC1)))
  for (u in 1:(length(europe.bv.centroid.summary$med_PC1))) {
    minus<-europe.bv.centroid.summary$med_PC1[u]
    minus.vec<-df.rot$PC1
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(europe.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  europe.bv.cent.dif.PC1[[names(europe.bv.pathway.rotations)[t]]]<-df.new
}

# collapse list into large dataframe with new column to identify list element
# origin (i.e. which pathway each v number corresponds to)
europe.bv.cent.dif.PC1.all<-bind_rows(europe.bv.cent.dif.PC1,.id = "pathway")

# transform data into 'tidy' format for ggplot faceting & make boxplots of
# difference between PC1 taxon centroid and PC1 pathway centroids, for each 
# taxon
europe.bv.plot.PC1<-europe.bv.cent.dif.PC1.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC1 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",
                    values = pathway.colors[sort(unique(europe.bv.cent.dif.PC1.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.036),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("europe_bv_PC1_difference.png",units = "in",height = 6,width = 15,res = 400)
europe.bv.plot.PC1
# dev.off()

# same for PC 2 
europe.bv.cent.dif.PC2<-list()
for (t in 1:length(names(europe.bv.pathway.rotations))) {
  df.rot<-as.data.frame(europe.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(europe.bv.centroid.summary$med_PC2)))
  for (u in 1:(length(europe.bv.centroid.summary$med_PC2))) {
    minus<-europe.bv.centroid.summary$med_PC2[u]
    minus.vec<-df.rot$PC2
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(europe.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  europe.bv.cent.dif.PC2[[names(europe.bv.pathway.rotations)[t]]]<-df.new
}

europe.bv.cent.dif.PC2.all<-bind_rows(europe.bv.cent.dif.PC2,.id = "pathway")

europe.bv.plot.PC2<-europe.bv.cent.dif.PC2.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC2 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(europe.bv.cent.dif.PC2.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.036),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("europe_bv_PC2_difference.png",units = "in",height = 6,width = 15,res = 400)
europe.bv.plot.PC2
# dev.off()

# same for PC3
europe.bv.cent.dif.PC3<-list()
for (t in 1:length(names(europe.bv.pathway.rotations))) {
  df.rot<-as.data.frame(europe.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(europe.bv.centroid.summary$med_PC3)))
  for (u in 1:(length(europe.bv.centroid.summary$med_PC3))) {
    minus<-europe.bv.centroid.summary$med_PC3[u]
    minus.vec<-df.rot$PC3
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(europe.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  europe.bv.cent.dif.PC3[[names(europe.bv.pathway.rotations)[t]]]<-df.new
}

europe.bv.cent.dif.PC3.all<-bind_rows(europe.bv.cent.dif.PC3,.id = "pathway")

europe.bv.plot.PC3<-europe.bv.cent.dif.PC3.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC3 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 2,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(europe.bv.cent.dif.PC3.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.036),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("europe_bv_PC3_difference.png",units = "in",height = 6,width = 15,res = 400)
europe.bv.plot.PC3
# dev.off()

############# both centroids

# make a vector of taxa and pathways for both datasets
both.bv.tax.vec <- tax.table[rownames(filt.bv.both),2]
both.bv.path.vec <- pathway.table[rownames(filt.bv.both),2]

# assign V numbers to vector elements
names(both.bv.tax.vec) <- rownames(filt.bv.both)
names(both.bv.path.vec) <- rownames(filt.bv.both)

# get number of genes corresponding to each taxon and pathway 
both.bv.no.genes.tax<-vector()
for (i in levels(factor(both.bv.tax.vec))) {
  both.bv.no.genes.tax[i]<-length(which(both.bv.tax.vec==i))
}

both.bv.no.genes.path<-vector()
for (j in levels(factor(both.bv.path.vec))) {
  both.bv.no.genes.path[j]<-length(which(both.bv.path.vec==j))
}

# pull all species represented by >100 genes and top 20 pathways (desc. order)
both.bv.no.genes.100<-names(which(both.bv.no.genes.tax >100))
both.bv.no.genes.20<-names(head(sort(both.bv.no.genes.path,decreasing = TRUE),21))

# create a list of named vectors containing indices of all v numbers corresponding
# to each taxa represented by >100 genes
both.bv.taxa.indices<-list()
for (spp in 1:length(both.bv.no.genes.100)) {
  both.bv.taxa.indices[[both.bv.no.genes.100[spp]]]<-which(both.bv.tax.vec==both.bv.no.genes.100[spp])
}

# create a list of named vectors containing indices of all v numbers corresponding
# to each of the top 20 pathways
both.bv.pathway.indices<-list()
for (path in 1:length(both.bv.no.genes.20)) {
  both.bv.pathway.indices[[both.bv.no.genes.20[path]]]<-which(both.bv.path.vec==both.bv.no.genes.20[path])
}

# finds the V number indices present in the both bv feature table that are not
# present in the list of species represented by >100 genes (i.e. all other genes)
both.bv.taxa.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.both)),
                                         unlist(both.bv.taxa.indices))

# finds the V number indices present in the both bv feature table that are not
# present in the list of top 20 pathways (i.e. all other pathways)
both.bv.pathway.indices[["Other"]]<-setdiff(c(1:nrow(filt.bv.both)),
                                            unlist(both.bv.pathway.indices))

# for each species, get the indices of all the V numbers with pathway = "Ribosome"
# in the both.bv.taxa vector (i.e. genes we want to use to calculate 'centroid')
both.bv.centroid.indices<-list()
for (tax in 1:(length(names(both.bv.taxa.indices))-1)) {
  vec.indices<-which(both.bv.taxa.indices[[tax]] %in% both.bv.pathway.indices[["Ribosome"]])
  vec.toAdd<-vector()
  for (index in 1:length(vec.indices)) {
    vec.toAdd<-append(vec.toAdd,both.bv.taxa.indices[[tax]][[vec.indices[index]]],after = length(vec.toAdd))
  }
  both.bv.centroid.indices[[names(both.bv.taxa.indices)[[tax]]]]<-vec.toAdd
}

# make a list of dataframes containing the rotation values for the corresponding
# species centroids (i.e. ribosomal genes for each species)
both.bv.taxa.centroids<-list()
for(tax in 1:length(names(both.bv.centroid.indices))){
  both.bv.taxa.centroids[[names(both.bv.centroid.indices)[[tax]]]]<-
    as.data.frame(both.bv.rotation[both.bv.centroid.indices[[tax]],])
}

# make a blank dataframe containing median PC1/PC2/PC3 values for each taxon
# (i.e. taxon centroids on the PCA plot)
both.bv.centroid.summary<-as.data.frame(cbind(med_PC1=c(rep(NA,length(both.bv.taxa.centroids))),
                                              med_PC2=c(rep(NA,length(both.bv.taxa.centroids))),
                                              med_PC3=(rep(NA,length(both.bv.taxa.centroids)))),
                                        row.names=names(both.bv.centroid.indices))

# fill in summary dataframe with values from the list of rotation dataframes
# for each taxon
for (tax in 1:length(rownames(both.bv.centroid.summary))) {
  df.tax<-as.data.frame(both.bv.taxa.centroids[[tax]])
  med.PC1<-median(both.bv.taxa.centroids[[tax]]$PC1)
  med.PC2<-median(both.bv.taxa.centroids[[tax]]$PC2)
  med.PC3<-median(both.bv.taxa.centroids[[tax]]$PC3)
  both.bv.centroid.summary[tax,1]<-med.PC1
  both.bv.centroid.summary[tax,2]<-med.PC2
  both.bv.centroid.summary[tax,3]<-med.PC3
}

# make list of dataframes containing rotations for all V numbers in each pathway
both.bv.pathway.rotations<-list()
for(tax in 2:21){
  both.bv.pathway.rotations[[names(both.bv.pathway.indices)[[tax]]]]<-
    as.data.frame(both.bv.rotation[both.bv.pathway.indices[[tax]],])
}

# make list of dataframes (1 per taxon) containing values denoting the 
# difference between PC1 rotation of each V number and PC1 taxon centroid (i.e.
# median PC1 value for that taxon). Note: values raised to power of 2 and sq.
# rooted to get around +/- values
both.bv.cent.dif.PC1<-list()
for (t in 1:length(names(both.bv.pathway.rotations))) {
  df.rot<-as.data.frame(both.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(both.bv.centroid.summary$med_PC1)))
  for (u in 1:(length(both.bv.centroid.summary$med_PC1))) {
    minus<-both.bv.centroid.summary$med_PC1[u]
    minus.vec<-df.rot$PC1
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(both.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  both.bv.cent.dif.PC1[[names(both.bv.pathway.rotations)[t]]]<-df.new
}

# collapse list into large dataframe with new column to identify list element
# origin (i.e. which pathway each v number corresponds to)
both.bv.cent.dif.PC1.all<-bind_rows(both.bv.cent.dif.PC1,.id = "pathway")

# transform data into 'tidy' format for ggplot faceting & make boxplots of
# difference between PC1 taxon centroid and PC1 pathway centroids, for each 
# taxon
both.bv.plot.PC1<-both.bv.cent.dif.PC1.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC1 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 3,scales = "fixed")+
  scale_fill_manual(name="Pathway",
                    values = pathway.colors[sort(unique(both.bv.cent.dif.PC1.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.032),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("both_bv_PC1_difference.png",units = "in",height = 6,width = 15,res = 400)
both.bv.plot.PC1
# dev.off()

# same for PC 2 
both.bv.cent.dif.PC2<-list()
for (t in 1:length(names(both.bv.pathway.rotations))) {
  df.rot<-as.data.frame(both.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(both.bv.centroid.summary$med_PC2)))
  for (u in 1:(length(both.bv.centroid.summary$med_PC2))) {
    minus<-both.bv.centroid.summary$med_PC2[u]
    minus.vec<-df.rot$PC2
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(both.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  both.bv.cent.dif.PC2[[names(both.bv.pathway.rotations)[t]]]<-df.new
}

both.bv.cent.dif.PC2.all<-bind_rows(both.bv.cent.dif.PC2,.id = "pathway")

both.bv.plot.PC2<-both.bv.cent.dif.PC2.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC2 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 3,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(both.bv.cent.dif.PC2.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.032),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("both_bv_PC2_difference.png",units = "in",height = 6,width = 15,res = 400)
both.bv.plot.PC2
# dev.off()

# same for PC3
both.bv.cent.dif.PC3<-list()
for (t in 1:length(names(both.bv.pathway.rotations))) {
  df.rot<-as.data.frame(both.bv.pathway.rotations[[t]])
  df.new<-as.data.frame(matrix(data = NA,
                               nrow = nrow(df.rot),
                               ncol = length(both.bv.centroid.summary$med_PC3)))
  for (u in 1:(length(both.bv.centroid.summary$med_PC3))) {
    minus<-both.bv.centroid.summary$med_PC3[u]
    minus.vec<-df.rot$PC3
    names(minus.vec)<-rownames(t)
    minus.result<-minus.vec-minus
    minus.trans<-sqrt((minus.result^2))
    df.new[,u]<-minus.trans
  }
  colnames(df.new)<-rownames(both.bv.centroid.summary)
  rownames(df.new)<-rownames(df.rot)
  both.bv.cent.dif.PC3[[names(both.bv.pathway.rotations)[t]]]<-df.new
}

both.bv.cent.dif.PC3.all<-bind_rows(both.bv.cent.dif.PC3,.id = "pathway")

both.bv.plot.PC3<-both.bv.cent.dif.PC3.all%>%
  pivot_longer(cols = -pathway,
               names_to = "species",
               values_to = "dist")%>%
  ggplot(aes(x=pathway,y=dist))+
  geom_point(col=rgb(0.5,0.5,0.5,0.1),size=1)+
  geom_boxplot(aes(fill=pathway),outlier.shape = NA)+
  labs(x= "Pathway",y="PC3 difference (taxon centroid - pathway centroid)")+
  facet_wrap(~species,nrow = 3,scales = "fixed")+
  scale_fill_manual(name="Pathway",values = pathway.colors[sort(unique(both.bv.cent.dif.PC3.all$pathway)),])+
  scale_y_continuous(limits = c(0,0.032),expand = c(0,0))+
  theme_bw()+theme(axis.text.x = element_blank())+theme(legend.text = element_text(size = 10))

# png("both_bv_PC3_difference.png",units = "in",height = 6,width = 15,res = 400)
both.bv.plot.PC3
# dev.off()

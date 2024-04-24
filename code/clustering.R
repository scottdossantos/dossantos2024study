# https://cran.r-project.org/web/packages/umap/index.html
install.packages('umap')
install.packages('factoextra')
install.packages('cluster')
install.packages('zCompositions')
install.packages('ppclust')
install.packages('fclust')

load('Rdata/ko.both.Rda')

# heirarchical clustering
# dendrogram, based on distances of clr transformed data
# should be congruent with fuzzy clustering or k-means clustering
clusters <- hclust(dist(t(d.clr)))
plot(clusters)

# reading 
# https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giz107/5572529


##########
# kmeans clustering
# semi-unsupervised clustering
# we decide how many groups we believe
##########

library(factoextra)
library(cluster)

# do a quick clustering
# finds centroids that minimize within group SS (multivariate euclidian distance)
# in the multivariate PCA plot
km <- kmeans(t(d.clr), centers=4, nstart=25)

# scree plot showing win SS / total SS 
# implies 2 or 4 clusters
fviz_nbclust(t(d.clr), kmeans, method = "wss")

# show the plots with samples in each group
fviz_cluster(km, data=t(d.clr))

# try to use the 'gap statistic' to determine clusters
# look for local max
# this takes F.O.R.E.V.E.R
gap_stat <- clusGap(t(d.clr),
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)

#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

# we see that once we get to 3-4 clusters there is no more info to be gained

#############
# fuzzy clustering
# shows the probability of group membership with kmeans cluster
#############
library(ppclust)
# determine fuzzy cluster membership
# sqeuclidian is variance in distance
res.fcm <- fcm(t(d.clr), dmetric='sqeuclidean', centers=3) #as.data.frame(res.fcm$u) #summary(res.fcm)
res.fcm2 <- ppclust2(res.fcm, "kmeans")
tail(res.fcm$u) # examine wild-type split columns

factoextra::fviz_cluster(res.fcm2, data = t(d.clr), axes=c(1,2),
  ellipse.type = "norm", labelsize=10,  palette = "jco",
  repel = TRUE)

# probabilies for cluster membership
res.fcm$u

# plot group membership probabilities
# strong support for WT split into 2 groups
plot(1:86, res.fcm$u[,1]/res.fcm$u[,3], pch=19, cex=5, log='y')
# no support for observed SNF2 split 
points(1:86, res.fcm$u[,2]/res.fcm$u[,4], col='red' , pch=19, cex=5,)
abline(h=1, col='grey')


#######
# UMAP
#######

# https://cran.r-project.org/web/packages/umap/index.html
# install.packages('umap')

# orients data in PCA space to maximize distance between groups, but not the 
# original orientation
library(umap)

cols <- gsub(1, 'red', res.fcm$cluster)
cols <- gsub(2, 'blue', cols)
cols <- gsub(3, 'orange', cols)

# make the umap output
set.seed(18)
y.umap <- umap(t(d.clr), n_components=4)

plot(y.umap$layout[,1],y.umap$layout[,2], col=cols,pch=19, cex=2, xlab='dim 1', ylab='dim 2')

clusters <- hclust(dist(t(d.clr)) )

plot(clusters)



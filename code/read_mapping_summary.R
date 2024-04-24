# read mapping summary upstream of VIRGO
# intended for supplementary data
# shows that 6 samples from ENA had low total yields
#     because of a lot of human mapping
# shows that 1 samples from LDN had low total yields
#     because of human mapping or rRNA mapping
# LDN samples had more rRNA - perhaps the ENA was
#     depleted before deposition

# source("code/setup.R")

# paste the paste function from setup with the new location

d.map <- read.table(paste(locn,"code/species_colors.txt", sep=""), header=T, row.names=2, sep="\t", stringsAsFactors=F)

d.num <- apply(d.map[,2:6], 2, as.numeric)
d.num <- rbind(d.num, rep(1e5, 5), rep(2e8, 5))

#parcoord(d.num, var.label=T, log="y", col=c(rep("red", 24),rep("blue",22)) )

# stripchart
stripchart(list(d.num[1:24,1],d.num[1:24,5]), col=rgb(1,0,0,0.7),
    method="jitter", vertical=T, pch=19, log="y",
    ylim=c(min(d.num[1:46,5]), max(d.num[1:46,1])),
    group.names=c("input", "mapped"))


stripchart(list(d.num[25:46,1],d.num[25:46,5]), col=rgb(0,0,1,0.7), method="jitter", vertical=T, pch=19, add=T)

 legend(1.4, 6e7, legend=c("LDN", "ENA"), pch=19, col=c("red", "blue"))




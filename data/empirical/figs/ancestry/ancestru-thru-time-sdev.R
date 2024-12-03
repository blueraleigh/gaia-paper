library(sf)

Z = readRDS("../../analysis/results/ancestry-thru-time-subset-avg-chr18p.rds")
landgrid = st_read("../../data/landgrid.gpkg")
land = st_transform(st_read("../../data/land.gpkg"), crs=st_crs(landgrid))

sites = st_centroid(landgrid$geom)

wrld = st_read("../../data/TM_WORLD_BORDERS-0.1.gpkg")
wrld = st_make_valid(wrld)

sites = st_transform(sites, st_crs(wrld))

SUBREGION = wrld$SUBREGION

site_map = wrld$SUBREGION[st_nearest_feature(sites, wrld)]
# map Taiwan to Eastern Asia
site_map[which(site_map == 0)] = 30L
# separate out Siberia east of Ural Mountains into a Northern Asia
site_map[site_map == 151L & (sapply(sites, "[", 1) > 59 | 
    sapply(sites, "[", 1) < -170)] = max(site_map) + 1L

sites = st_centroid(landgrid$geom)

site_coords = st_coordinates(sites)

# UN M49 regions
UNM49 = c(
`18`="#fe9929", # 5 southern africa
`17`="#ec7014", # 4 middle africa
`14`="#cc4c02", # 2 eastern africa
`11`="#993404", # 1 western africa
`15`="#662506", # 3 north africa

`39`="#4eb3d3",  # 9 southern europe
`151`="#2b8cbe", # 12 eastern europe
`154`="#0868ac", # 13 northern europe
`155`="#084081",  # 14 western europe

`30`="#99d8c9",  # 6 eastern asia
`34`="#66c2a4",  # 7 southern asia
`35`="#41ae76",  # 8 southeastern asia
`143`="#238b45", # 10 central asia
`145`="#006d2c", # 11 western asia
`156`="#00441b"  # 15 northern asia
)


ME = 145
EU = c(39,151,154,155)
AS = c(30,34,35,143,156)
AF = c(11,14,15,17,18)

ME_I = which(site_map %in% ME)
EU_I = which(site_map %in% EU)
AS_I = which(site_map %in% AS)
AF_I = which(site_map %in% AF)


ord = c(AF_I,ME_I,AS_I,EU_I)

cols = c("#a6611a","#dfc27d","#80cdc1","#018571")
leg = c("Africa", "Europe", "Asia", "Middle East")
colv = character(177)
colv[AF_I] = cols[1]
colv[EU_I] = cols[2]
colv[AS_I] = cols[3]
colv[ME_I] = cols[4]


kya = character(201)
kya[c(5,41,81,201)] = c("10","100","200","500")
par(mfcol=c(4,4), mar=c(1,3,1,1),oma=c(2,2,2,1))

rg = 1
for (tm in c(5,41,81,201))
{
u = Z[,rg,tm]+attr(Z,"std.error")[,rg,tm]*sqrt(100)
l = Z[,rg,tm]-attr(Z,"std.error")[,rg,tm]*sqrt(100)
plot(Z[,rg,tm][ord], pch=19, cex=0.6, las=1, bty='l', ylab="Sample ancestry", 
    ylim=c(0,max(u)), xaxt='n', xlab='',col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, u[ord], col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, l[ord], col=colv[ord])
if (tm == 5)
{
    title("Africa", adj=0)
}
if (tm == 81)
{
mtext("Sample ancestry", 2, line=3, at=0.2)
}
legend("topright", legend=sprintf("%s kya", kya[tm]), bty="n", cex=0.7)
}

rg = 2
for (tm in c(5,41,81,201))
{
u = Z[,rg,tm]+attr(Z,"std.error")[,rg,tm]*sqrt(100)
l = Z[,rg,tm]-attr(Z,"std.error")[,rg,tm]*sqrt(100)
plot(Z[,rg,tm][ord], pch=19, cex=0.6, las=1, bty='l', ylab="Sample ancestry", 
    ylim=c(0,max(u)), xaxt='n', xlab='',col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, u[ord], col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, l[ord], col=colv[ord])
if (tm == 5)
{
    title("Europe", adj=0)
}
legend("topright", legend=sprintf("%s kya", kya[tm]), bty="n", cex=0.7)
}

rg = 3
for (tm in c(5,41,81,201))
{
u = Z[,rg,tm]+attr(Z,"std.error")[,rg,tm]*sqrt(100)
l = Z[,rg,tm]-attr(Z,"std.error")[,rg,tm]*sqrt(100)
plot(Z[,rg,tm][ord], pch=19, cex=0.6, las=1, bty='l', ylab="Sample ancestry", 
    ylim=c(0,max(u)), xaxt='n', xlab='',col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, u[ord], col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, l[ord], col=colv[ord])
if (tm==201)
{
legend(0, -0.02, legend=leg, pch=19, col=cols, lty=1, ncol=4, bty='n',
    inset=-.2,xpd=NA, xjust=0.5)
}
if (tm == 5)
{
    title("Asia", adj=0)
}
legend("topright", legend=sprintf("%s kya", kya[tm]), bty="n", cex=0.7)
}


rg = 4
for (tm in c(5,41,81,201))
{
u = Z[,rg,tm]+attr(Z,"std.error")[,rg,tm]*sqrt(100)
l = Z[,rg,tm]-attr(Z,"std.error")[,rg,tm]*sqrt(100)
plot(Z[,rg,tm][ord], pch=19, cex=0.6, las=1, bty='l', ylab="Sample ancestry", 
    ylim=c(0,max(u)), xaxt='n', xlab='',col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, u[ord], col=colv[ord])
segments(1:177, Z[,rg,tm][ord], 1:177, l[ord], col=colv[ord])
if (tm == 5)
{
    title("Middle East", adj=0)
}
legend("topright", legend=sprintf("%s kya", kya[tm]), bty="n", cex=0.7)
}

dev.print(pdf, file="ancestry-thru-time-sdev.pdf")

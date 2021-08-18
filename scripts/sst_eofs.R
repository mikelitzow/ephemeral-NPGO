# global sst EOFs
# exploratory analysis

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(FactoMineR)
library(nlme)
library(MuMIn)
library(oce)
library(ggpubr)

# set palette
new.col <- oceColorsPalette(64)

# set theme
theme_set(theme_bw())

# load and process SST data
nc <- nc_open("./data/nceiErsstv5_d29b_802d_76af.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 20-68 deg. N, 120-250 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

x; y


SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


m <- months(d)  # Extracts months from the date vector
yr <- years(d)

# limit to 1949 - 2021 (to match SLP data)
m <- m[yr %in% 1949:2021]
SST <- SST[yr %in% 1949:2021,]
d <- d[yr %in% 1949:2021]
yr <- yr[yr %in% 1949:2021]

# restrict to NE Pacific
SST <- SST[,lat %in% 20:64 & lon >= 150]

# clean up x, y, lat, lon
y <- y[y %in% 20:64]
x <- x[x >= 150]

# reset lat/lon
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   


# and plot 
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col)
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# looks good!

# identify columns in SST matrix corresponding to land
land <- is.na(colMeans(SST)) 

# For analysis, we only use the columns of the matrix with non-missing values:
X <- SST[,!land] 

# remove seasonal means
f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 

xtra <- 12*((length(d)/12)-floor(length(d)/12)) # add trailing months

mu <- rbind(mu, mu[1:xtra,])
anom <- X - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}

# get a vector of weights (square root of the cosine of latitude)
lat.weights <- lat[!land]
weight <- sqrt(cos(lat.weights*pi/180))

# limit to FMA 
# (peak PDO months from Newman et al. 2016, will match with SLP for NDJ)

d.FMA <- d[m %in% c("Feb", "Mar", "Apr")]
yr.FMA <- yr[m %in% c("Feb", "Mar", "Apr")]
anom.detr.FMA <- anom.detr[m %in% c("Feb", "Mar", "Apr"),]

# these are the four eras we want to calculate...
length(2009:2021)
length(1989:2008)
length(1969:1988)
length(1949:1968)

# EOF by era 
# weighting the columns
EOF.all <- svd.triplet(cov(anom.detr.FMA), col.w=weight)
EOF.1 <- svd.triplet(cov(anom.detr.FMA[yr.FMA %in% 1949:1968,]), col.w=weight)
EOF.2 <- svd.triplet(cov(anom.detr.FMA[yr.FMA %in% 1969:1988,]), col.w=weight)
EOF.3 <- svd.triplet(cov(anom.detr.FMA[yr.FMA %in% 1989:2008,]), col.w=weight)
EOF.4 <- svd.triplet(cov(anom.detr.FMA[yr.FMA %in% 2009:2021,]), col.w=weight)

# get loadings for EOF1-2 by era
eig.1.all <- EOF.all$U[,1]
eig.2.all <- EOF.all$U[,2]

eig.1.1 <- -EOF.1$U[,1] # reversing sign to match full time series pattern!
eig.2.1 <- EOF.1$U[,2]

eig.1.2 <- -EOF.2$U[,1] # reversing sign to match full time series pattern!
eig.2.2 <- EOF.2$U[,2]

eig.1.3 <- EOF.3$U[,1] 
eig.2.3 <- EOF.3$U[,2]

eig.1.4 <- EOF.4$U[,1] 
eig.2.4 <- EOF.4$U[,2]


# get % variance explained by era
var.all <- 100*round(prop.table(EOF.all$vs),3)
var.1 <- 100*round(prop.table(EOF.1$vs),3)
var.2 <- 100*round(prop.table(EOF.2$vs),3)
var.3 <- 100*round(prop.table(EOF.3$vs),3)
var.4 <- 100*round(prop.table(EOF.4$vs),3)

# set the limit for plotting 
lim.1 <- range(eig.1.all, eig.1.1, eig.1.2, eig.1.3, eig.1.4)
lim.2 <- range(eig.2.all, eig.2.1, eig.2.2, eig.2.3, eig.2.4)

lim.1
lim.2

# and plot

png("./figs/sst_EOF.png", 4, 10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# Full time series!
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1949-2021 (", var.all[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- eig.2.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1949-2021 (", var.all[2], "%)", sep=""), cex=0.8)

#################
# Era 1!
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1949-1968 (", var.1[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- eig.2.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1949-1969 (", var.1[2], "%)", sep=""), cex=0.8)

#############
# Era 2
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1969-1988 (", var.2[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- eig.2.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1969-1988 (", var.2[2], "%)", sep=""), cex=0.8)


#######
# Era 3
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1989-2008 (", var.3[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- eig.2.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1989-2008 (", var.3[2], "%)", sep=""), cex=0.8)

#############
# Era 4
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 2009-2021 (", var.4[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- eig.2.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 2009-2021 (", var.4[2], "%)", sep=""), cex=0.8)

dev.off()


# plot distributions of changes in loadings from full time series for each era
diff.1.1 <- eig.1.all - eig.1.1
diff.1.2 <- eig.1.all - eig.1.2
diff.1.3 <- eig.1.all - eig.1.3
diff.1.4 <- eig.1.all - eig.1.4

diff.2.1 <- eig.2.all - eig.2.1
diff.2.2 <- eig.2.all - eig.2.2
diff.2.3 <- eig.2.all - eig.2.3
diff.2.4 <- eig.2.all - eig.2.4

# combine into one df
plot.diff <- data.frame(mode = rep(c("Mode 1", "Mode 2"), each = 4*length(diff.1.1)),
                        era = rep(rep(c("1949-1968", "1969-1988", "1989-2008", "2009-2021"), each = length(diff.1.1)), 2),
                        difference = c(diff.1.1, diff.1.2, diff.1.3, diff.1.4, diff.2.1, diff.2.2, diff.2.3, diff.2.4))



m1 <- ggplot(filter(plot.diff, mode == "Mode 1"), aes(difference)) +
  geom_density(fill = "grey", lty = 0) +
  facet_wrap(~era, ncol = 1) +
  ggtitle("Mode 1")

m2 <- ggplot(filter(plot.diff, mode == "Mode 2"), aes(difference)) +
  geom_density(fill = "grey", lty = 0) +
  facet_wrap(~era, ncol = 1) +
  ggtitle("Mode 2")

png("./figs/sst_EOF_loading_differences_by_era.png", width=4, height=5.5, units='in', res=300)
ggarrange(m1, m2, ncol=2)
dev.off()

# examine SD in differences 
sd1 <- data.frame(mode = "Mode 1",
                  era = c("1949-1968", "1969-1988", "1989-2008", "2009-2021"),
                  SD = c(sd(diff.1.1), sd(diff.1.2), sd(diff.1.3), sd(diff.1.4)))

sd2 <- data.frame(mode = "Mode 2",
                  era = c("1949-1968", "1969-1988", "1989-2008", "2009-2021"),
                  SD = c(sd(diff.2.1), sd(diff.2.2), sd(diff.2.3), sd(diff.2.4)))

plot.sd <- rbind(sd1, sd2)

ggplot(plot.sd, aes(era, SD, fill = mode)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab("Standard deviation")

ggsave("./figs/SD_SST_loading_differences_by_era.png", width = 4.5, height = 5, units = 'in')


# some thoughts:
# mode 2 variance about twice as great in some areas, similar to mode 1 in other eras
# would be good to examine correlations between PC1 and PC2 (defined for era 1) across different eras
# also need SDE to relate PC1/2 of SLP/SST
# and...look at time-dependent correlations between mode 1/2 and upwelling TS as in Di Lorenzo et al. 2008


# next steps:
# 1) get full time series PC scores and correlate these with EOFs for each era

# calculate PC scores
pc1 <- anom.detr[yr %in% 1914:2013,] %*% EOF.all$U[,1]
pc2 <- anom.detr[yr %in% 1914:2013,] %*% EOF.all$U[,2]
pc3 <- anom.detr[yr %in% 1914:2013,] %*% EOF.all$U[,3]
pc4 <- anom.detr[yr %in% 1914:2013,] %*% EOF.all$U[,4]

# and get a shorter year index for subsetting pcs
short.yr <- yr[yr %in% 1914:2013]
names(pc1) <- names(pc2) <- names(pc3) <- names(pc4) <- short.yr

# get correlations with era-specific eigenvectors...
cor1.a <- cor2.a <- cor3.a <- cor4.a <- NA

for(i in 1:ncol(anom.detr)){
  
  cor1.a[i] <- cor(pc1, anom.detr[yr %in% 1914:2013, i])
  cor2.a[i] <- cor(pc2, anom.detr[yr %in% 1914:2013, i]) 
  cor3.a[i] <- cor(pc3, anom.detr[yr %in% 1914:2013, i])
  cor4.a[i] <- cor(pc4, anom.detr[yr %in% 1914:2013, i])
}

cor1.1 <- cor2.1 <- cor3.1 <- cor4.1 <- NA

for(i in 1:ncol(anom.detr)){
  
  cor1.1[i] <- cor(pc1[short.yr %in% 1914:1938], anom.detr[yr %in% 1914:1938, i])
  cor2.1[i] <- cor(pc2[ short.yr %in% 1914:1938], anom.detr[yr %in% 1914:1938, i]) 
  cor3.1[i] <- cor(pc3[short.yr %in% 1914:1938], anom.detr[yr %in% 1914:1938, i])
  cor4.1[i] <- cor(pc4[short.yr %in% 1914:1938], anom.detr[yr %in% 1914:1938, i])
}

cor1.2 <- cor2.2 <- cor3.2 <- cor4.2 <- NA

for(i in 1:ncol(anom.detr)){
  
  cor1.2[i] <- cor(pc1[short.yr %in% 1939:1963], anom.detr[yr %in% 1939:1963, i])
  cor2.2[i] <- cor(pc2[short.yr %in% 1939:1963], anom.detr[yr %in% 1939:1963, i]) 
  cor3.2[i] <- cor(pc3[short.yr %in% 1939:1963], anom.detr[yr %in% 1939:1963, i])
  cor4.2[i] <- cor(pc4[short.yr %in% 1939:1963], anom.detr[yr %in% 1939:1963, i])
}

cor1.3 <- cor2.3 <- cor3.3 <- cor4.3 <- NA

for(i in 1:ncol(anom.detr)){
  
  cor1.3[i] <- cor(pc1[short.yr %in% 1964:1988], anom.detr[yr %in% 1964:1988, i])
  cor2.3[i] <- cor(pc2[short.yr %in% 1964:1988], anom.detr[yr %in% 1964:1988, i]) 
  cor3.3[i] <- cor(pc3[short.yr %in% 1964:1988], anom.detr[yr %in% 1964:1988, i])
  cor4.3[i] <- cor(pc4[short.yr %in% 1964:1988], anom.detr[yr %in% 1964:1988, i])
}

cor1.4 <- cor2.4 <- cor3.4 <- cor4.4 <- NA

for(i in 1:ncol(anom.detr)){
 # i <- 1
  cor1.4[i] <- cor(pc1[short.yr %in% 1989:2013], anom.detr[yr %in% 1989:2013, i])
  cor2.4[i] <- cor(pc2[short.yr %in% 1989:2013], anom.detr[yr %in% 1989:2013, i]) 
  cor3.4[i] <- cor(pc3[short.yr %in% 1989:2013], anom.detr[yr %in% 1989:2013, i])
  cor4.4[i] <- cor(pc4[short.yr %in% 1989:2013], anom.detr[yr %in% 1989:2013, i])
}



# reset the limits for plotting 
lim.1 <- range(cor1.a, cor1.1, cor1.2, cor1.3, cor1.4)
lim.2 <- range(cor2.a, cor2.1, cor2.2, cor2.3, cor2.4)
lim.3 <- range(cor3.a, cor3.1, cor3.2, cor3.3, cor3.4)
lim.4 <- range(cor4.a, cor4.1, cor4.2, cor4.3, cor4.4)

lim.1
lim.2
lim.3
lim.4

lims <- range(lim.1,
              lim.2,
              lim.3,
              lim.4)

png("correlations with PCs by era and full time series.png", 11, 8, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfcol=c(4,5), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# Full time series!
z <- rep(NA, ncol(SST))
z[!land] <- cor1.a
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC1 correlation 1914-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor2.a
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC2 correlation 1914-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor3.a
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC3 correlation 1914-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor4.a
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC4 correlation 1914-2013", cex=0.8)

# Era 1!
z <- rep(NA, ncol(SST))
z[!land] <- cor1.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC1 correlation 1914-1938", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor2.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC2 correlation 1914-1938", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor3.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC3 correlation 1914-1938", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor4.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC4 correlation 1914-1938", cex=0.8)

#############
# Era 2
z <- rep(NA, ncol(SST))
z[!land] <- cor1.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC1 correlation 1939-1963", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor2.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC2 correlation 1939-1963", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor3.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC3 correlation 1939-1963", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor4.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC4 correlation 1939-1963", cex=0.8)

#######
# Era 3
z <- rep(NA, ncol(SST))
z[!land] <- cor1.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC1 correlation 1964-1988", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor2.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC2 correlation 1964-1988", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor3.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC3 correlation 1964-1988", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor4.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC4 correlation 1964-1988", cex=0.8)

#############
# Era 4
z <- rep(NA, ncol(SST))
z[!land] <- cor1.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC1 correlation 1989-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor2.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC2 correlation 1989-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor3.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC3 correlation 1989-2013", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- cor4.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lims[2], lims[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, add=T, lwd=1)
mtext("PC4 correlation 1989-2013", cex=0.8)

dev.off()


# next steps:
# 2) regress SST anomalies onto MEI / AMO / PDO / NPGO for each era!

# calculate AMO following Trenberth and Shea 2006
# Mean N. Atlantic SST, 0º-60ºN, 0º-80ºW 
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1914-01-01):1:(2018-11-01T00:00:00Z)][(0.0):1:(0.0)][(0):1:(60)][(280):1:(358)]", "~amo.sst")
# load and process SST data
nc.amo <- nc_open("~amo.sst")

# extract dates

ncvar_get(nc.amo, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc.amo, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 54-62 deg. N, 200-226 deg. E
x.amo <- ncvar_get(nc.amo, "longitude")
y.amo <- ncvar_get(nc.amo, "latitude")

SST.amo <- ncvar_get(nc.amo, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST.amo <- aperm(SST.amo, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST.amo <- matrix(SST.amo, nrow=dim(SST.amo)[1], ncol=prod(dim(SST.amo)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat.amo <- rep(y.amo, length(x.amo))   
lon.amo <- rep(x.amo, each = length(y.amo))   
dimnames(SST.amo) <- list(as.character(d), paste("N", lat.amo, "E", lon.amo, sep=""))

# plot mean temperature pattern to check
png("check.png", 6,8, units="in", res=300)
SST.mean.amo <- colMeans(SST.amo)
z <- t(matrix(SST.mean.amo,length(y.amo)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x.amo,y.amo,z, col=tim.colors(64))
contour(x.amo, y.amo, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)
dev.off()

# need to black out some outlying areas
hudson <- lat.amo > 50 & lon.amo < 286 
SST.amo[,hudson] <- NA  
  
great <- lat.amo >40 & lon.amo <=284 
SST.amo[,great] <- NA

pac <- lat.amo < 10 & lon.amo <= 282
SST.amo[,pac] <- NA

med <- lon.amo > 354 & lat.amo == 36 
SST.amo[,med] <- NA


# get raw average AMO

# get a vector of weights (square root of the cosine of latitude)
amo.weight <- sqrt(cos(lat.amo*pi/180))

ff <- function(x) weighted.mean(x, amo.weight, na.rm=T)
amo.raw <- apply(SST.amo, 1, ff)

# now get monthly mean and remove
amo.months <- months(d)

mu <- tapply(amo.raw, amo.months, mean)  

mu <- mu[rep(1:12, floor(length(d)/12))] 

xtra <- 12*((length(d)/12)-floor(length(d)/12)) # add trailing months

mu <- c(mu, mu[1:xtra])
amo <- amo.raw - mu   # compute matrix of anomalies

# get global SST mean
ff <- function(x) weighted.mean(x, weight, na.rm=T)

glob.mean <- apply(X[yr >= 1914,], 1, ff)

glob.mean <- glob.mean-mean(glob.mean)

amo <- amo-glob.mean

plot(1:length(amo), amo, type="l", col="grey")
abline(h=0)

amo.sm <- rollmean(amo, 121, mean, fill=NA)

# set up data frame for analysis of 1964:2013 correlations

cor.dat <- data.frame(year=yr[yr %in% 1964:2013], month=1:12, mei=NA, amo=amo[years(names(amo)) %in% 1964:2013], pdo=NA, npgo=NA, pc1=NA, pc2=NA, pc3=NA, pc4=NA)

mei <- read.csv("mei.csv")
mei <- mei %>%
  select(5:16)

# rename months as numbers
names(mei) <- 1:12
mei <- gather(mei)
mei$key <- as.numeric(mei$key)

mei$year <- as.numeric(rep(1950:2017))

mei <- arrange(mei, year)

cor.dat$mei <- mei$value[mei$year %in% 1964:2013]

# pdo
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

cor.dat$pdo <- pdo$value[pdo$YEAR %in% 1964:2013]

# and npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=818, fill=T, col.names = c("Year", "month", "value"))

cor.dat$npgo <- npgo$value[npgo$Year %in% 1964:2013]

# and add era-specific PC1-4
cor.dat$pc1[cor.dat$year %in% 1964:1988] <- anom.detr[yr %in% 1964:1988,] %*% EOF.3$U[,1]
cor.dat$pc2[cor.dat$year %in% 1964:1988] <- anom.detr[yr %in% 1964:1988,] %*% EOF.3$U[,2]
cor.dat$pc3[cor.dat$year %in% 1964:1988] <- anom.detr[yr %in% 1964:1988,] %*% EOF.3$U[,3]
cor.dat$pc4[cor.dat$year %in% 1964:1988] <- anom.detr[yr %in% 1964:1988,] %*% EOF.3$U[,4]

cor.dat$pc1[cor.dat$year %in% 1989:2013] <- anom.detr[yr %in% 1989:2013,] %*% EOF.4$U[,1]
cor.dat$pc2[cor.dat$year %in% 1989:2013] <- anom.detr[yr %in% 1989:2013,] %*% EOF.4$U[,2]
cor.dat$pc3[cor.dat$year %in% 1989:2013] <- anom.detr[yr %in% 1989:2013,] %*% EOF.4$U[,3]
cor.dat$pc4[cor.dat$year %in% 1989:2013] <- anom.detr[yr %in% 1989:2013,] %*% EOF.4$U[,4]

cor.dat$era <- 1
cor.dat$era[cor.dat$year > 1988] <- 2
cor.dat$era <- as.factor(cor.dat$era)
# and save so we don't have to re-run the EOFs!
write.csv(cor.dat, "global eof correlation data.csv")


plot.pc1 <- cor.dat %>%
  select(3:7,11) %>%
  gather(index, value, -pc1, -era)

plot.pc2 <- cor.dat %>%
  select(3:6,8,11) %>%
  gather(index, value, -pc2, -era)

plot.pc3 <- cor.dat %>%
  select(3:6,9,11) %>%
  gather(index, value, -pc3, -era)

plot.pc4 <- cor.dat %>%
  select(3:6,10,11) %>%
  gather(index, value, -pc4, -era)

quartz()
ggplot(plot.pc1, aes(value, pc1, color=era)) + geom_point() + facet_wrap(~index, scales="free") + 
  geom_smooth(method="lm", se=F)

quartz()
ggplot(plot.pc2, aes(value, pc2, color=era)) + geom_point() + facet_wrap(~index, scales="free") + 
  geom_smooth(method="lm", se=F)

quartz()
ggplot(plot.pc3, aes(value, pc3, color=era)) + geom_point() + facet_wrap(~index, scales="free") + 
  geom_smooth(method="lm", se=F)

quartz()
ggplot(plot.pc4, aes(value, pc4, color=era)) + geom_point() + facet_wrap(~index, scales="free") + 
  geom_smooth(method="lm", se=F)

# test for changing relationships with gls model
era <- cor.dat$era
for(xx in 3:6){
 # xx <- 3
  covar <- cor.dat[,xx]
  
  for(yy in 7:10){
   # yy <- 7
   resp <- cor.dat[,yy]
  
   mod.nonst <- gls(resp ~ covar*era, correlation = corAR1())
   mod.st <- gls(resp ~ covar, correlation = corAR1())
   
   print(paste(colnames(cor.dat)[yy], " = ", colnames(cor.dat)[xx], sep=""))
   print(AICc(mod.nonst, mod.st))
   print(summary(mod.nonst))
  }
}

# plot bar plots
cor1 <- cor(cor.dat[cor.dat$era==1, 3:10])
cor2 <- cor(cor.dat[cor.dat$era==2, 3:10])

pc.correlations <- data.frame(correlation=c(cor1[5:8,1:4], cor2[5:8,1:4]), index=rep(c("MEI", "AMO", "PDO", "NPGO"), each=4), pc=c("PC1", "PC2", "PC3", "PC4"), 
                              era=rep(c("1964-1988", "1989-2013"), each=16))

ggplot(pc.correlations, aes(index, correlation, fill=era)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~pc) + geom_hline(yintercept = 0, color="dark grey")


diag(cor1) <- diag(cor2) <- 0

index.correlations <- data.frame(correlation=c(cor1[1:4,1:4], cor2[1:4,1:4]), x=rep(c("MEI", "AMO", "PDO", "NPGO"), each=4), y=c("MEI", "AMO", "PDO", "NPGO"), 
                                                 era=rep(c("1964-1988", "1989-2013"), each=16))

ggplot(index.correlations, aes(x, correlation, fill=era)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~y) + geom_hline(yintercept = 0, color="dark grey")


# running 10-yr correlations between AMO and other indices
cor.dat$amo.mei <- cor.dat$amo.pdo <- cor.dat$amo.npgo <- NA
for(i in 120:600){
  cor.dat$amo.mei[i] <- cor(cor.dat$amo[(i-119):i], cor.dat$mei[(i-119):i])
  cor.dat$amo.pdo[i] <- cor(cor.dat$amo[(i-119):i], cor.dat$pdo[(i-119):i]) 
  cor.dat$amo.npgo[i] <- cor(cor.dat$amo[(i-119):i], cor.dat$npgo[(i-119):i])
}

dec.yr <- as.numeric(as.character(cor.dat$year)) + (cor.dat$month-0.5)/12

plot.dat <- cor.dat %>%
  select(12:14) %>%
  gather()
plot.dat$year <- dec.yr

ggplot(plot.dat, aes(year, value, color=key)) + geom_line() + geom_hline(yintercept = 0, col="dark grey")

plot(dec.yr, cor.dat$amo, type="l", col="grey")
lines(dec.yr, rollmean(cor.dat$amo, 120, fill=NA), col="red", lwd=1.3)
abline(h=0)

indices <- cor.dat %>%
  select(3:6,11) %>%
  gather(index, value, -era)

pcs <- cor.dat %>%
  select(7:10,11) %>%
  gather(axis, pc.score, -era)

plot.pcs <- left_join(indices, pcs)

png("global pc vs. index plots.png", 10,8, units="in", res=300)
ggplot(plot.pcs, aes(value, pc.score, color=era)) + geom_point() + facet_grid(axis~index, scales="free") + 
  geom_smooth(method="lm", se=F)
dev.off()

# plot regression of SST onto the four indices for each era
mei.1 <- mei.2 <- amo.1 <- amo.2 <- pdo.1 <- pdo.2 <- npgo.1 <- npgo.2 <- NA

for(i in 1:ncol(anom.detr)){
  #i <- 1
  temp <- anom.detr[yr %in% 1964:1988,i]
  mei.1[i] <- summary(lm(temp ~ cor.dat$mei[cor.dat$year %in% 1964:1988]))$coefficients[2,1]
  amo.1[i] <- summary(lm(temp ~ cor.dat$amo[cor.dat$year %in% 1964:1988]))$coefficients[2,1]
  pdo.1[i] <- summary(lm(temp ~ cor.dat$pdo[cor.dat$year %in% 1964:1988]))$coefficients[2,1]
  npgo.1[i] <- summary(lm(temp ~ cor.dat$npgo[cor.dat$year %in% 1964:1988]))$coefficients[2,1]
  
  temp <- anom.detr[yr %in% 1989:2013,i]
  mei.2[i] <- summary(lm(temp ~ cor.dat$mei[cor.dat$year %in% 1989:2013]))$coefficients[2,1]
  amo.2[i] <- summary(lm(temp ~ cor.dat$amo[cor.dat$year %in% 1989:2013]))$coefficients[2,1]
  pdo.2[i] <- summary(lm(temp ~ cor.dat$pdo[cor.dat$year %in% 1989:2013]))$coefficients[2,1]
  npgo.2[i] <- summary(lm(temp ~ cor.dat$npgo[cor.dat$year %in% 1989:2013]))$coefficients[2,1]
  
}

# and plot
mei.diff <- mei.2-mei.1
amo.diff <- amo.2-amo.1
pdo.diff <- pdo.2-pdo.1
npgo.diff <- npgo.2-npgo.1

dlim <- range(mei.diff, amo.diff, pdo.diff, npgo.diff)

png("global sst regressed on indices.png", 11, 8, units="in", res=300)
par(mar=c(1.5,2.5,1,0.5),  mfrow=c(4,3), oma=c(0,0,0,0.2))

zlim <- range(mei.1, mei.2)

z <- rep(NA, ncol(anom.detr))
z[!land] <- mei.1 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-MEI 1964-1988", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- mei.2 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-MEI 1989-2013", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- mei.diff
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-dlim[2],dlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("(1989-2013)-(1964-2013)", cex=0.8)

zlim <- range(amo.1, amo.2)
z <- rep(NA, ncol(anom.detr))
z[!land] <- amo.1 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-AMO 1964-1988", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- amo.2 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-AMO 1989-2013", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- amo.diff
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-dlim[2],dlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("(1989-2013)-(1964-2013)", cex=0.8)

zlim <- range(pdo.1, pdo.2)
z <- rep(NA, ncol(anom.detr))
z[!land] <- pdo.1 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-PDO 1964-1988", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- pdo.2 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-PDO 1989-2013", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- pdo.diff
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-dlim[2],dlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("(1989-2013)-(1964-2013)", cex=0.8)

zlim <- range(npgo.1, npgo.2)
z <- rep(NA, ncol(anom.detr))
z[!land] <- npgo.1 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-NPGO 1964-1988", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- npgo.2 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("SST-NPGO 1989-2013", cex=0.8)

z <- rep(NA, ncol(anom.detr))
z[!land] <- npgo.diff
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, zlim=c(-dlim[2],dlim[2]), ylab="", xlab="", yaxt="n", xaxt="n")
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F,add=T, lwd=1)
mtext("(1989-2013)-(1964-2013)", cex=0.8)

dev.off()

###
# next step..."blending" of EOFs...ie, if we calculate PCs based on 1964-1988 definition, are they now correlated?

# pc1 <- anom.detr[yr %in% 1914:2013,] %*% EOF.all$U[,1]
pc1 <- pc2 <- pc3 <- pc4 <- NA

use <- anom.detr[yr %in% 1964:2018,]
for(i in 1:nrow(use)){
 # i <- 1
 temp <- use[i,]
 pc1[i] <- 8763*weighted.mean(temp*EOF.3$U[,1], weight)
 pc2[i] <- 8763*weighted.mean(temp*EOF.3$U[,2], weight) 
 pc3[i] <- 8763*weighted.mean(temp*EOF.3$U[,3], weight)
 pc4[i] <- 8763*weighted.mean(temp*EOF.3$U[,4], weight)
}

pc.ts <- data.frame(year=as.numeric(as.character(yr[yr %in% 1964:2018])), month=c(rep(1:12, length.out=12*54), 1:11), pc1=pc1, pc2=pc2, pc3=pc3, pc4=pc4, 
                    pc1.pc2=NA, pc1.pc3=NA, pc1.pc4=NA, pc2.pc3=NA, pc2.pc4=NA, pc3.pc4=NA)

pc.ts$dec.yr <- pc.ts$year+(pc.ts$month-0.5)/12

for(i in 120:659){
  pc.ts$pc1.pc2[i] <- cor(pc.ts$pc1[(i-119):i], pc.ts$pc2[(i-119):i])
  pc.ts$pc1.pc3[i] <- cor(pc.ts$pc1[(i-119):i], pc.ts$pc3[(i-119):i])
  pc.ts$pc1.pc4[i] <- cor(pc.ts$pc1[(i-119):i], pc.ts$pc4[(i-119):i])
  
  pc.ts$pc2.pc3[i] <- cor(pc.ts$pc2[(i-119):i], pc.ts$pc3[(i-119):i])
  pc.ts$pc2.pc4[i] <- cor(pc.ts$pc2[(i-119):i], pc.ts$pc4[(i-119):i])
  
  pc.ts$pc3.pc4[i] <- cor(pc.ts$pc3[(i-119):i], pc.ts$pc4[(i-119):i])
}

plot.pc1 <- pc.ts %>%
  select(7:9) %>%
  gather()
plot.pc1$year <- pc.ts$dec.yr

ggplot(plot.pc1, aes(year, value, color=key)) + geom_line() + geom_hline(yintercept = 0, col="dark grey")

plot.pc2 <- pc.ts %>%
  select(7,10,11) %>%
  gather()
plot.pc2$year <- pc.ts$dec.yr

ggplot(plot.pc2, aes(year, value, color=key)) + geom_line() + geom_hline(yintercept = 0, col="dark grey")

plot.pc3 <- pc.ts %>%
  select(8,10,12) %>%
  gather()
plot.pc3$year <- pc.ts$dec.yr

ggplot(plot.pc3, aes(year, value, color=key)) + geom_line() + geom_hline(yintercept = 0, col="dark grey")


plot.pc4 <- pc.ts %>%
  select(9,11,12) %>%
  gather()
plot.pc4$year <- pc.ts$dec.yr

ggplot(plot.pc4, aes(year, value, color=key)) + geom_line() + geom_hline(yintercept = 0, col="dark grey")

# and...identify correlations with e.g., precip or streamflow and AMO/NPGO from the literature and show how they are no cross-correlated with each other?
# and...identify changing SLP regressions on either indices or sst EOFs
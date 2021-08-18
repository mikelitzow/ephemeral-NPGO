# slp EOFs
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

# load and process SLP data
nc <- nc_open("./data/hawaii_soest_f19d_3925_d70b_8736_2da9_5f29.nc")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 20-70 deg. N, 120-250 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

x; y


SLP <- ncvar_get(nc, "slp", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


m <- months(d)  # Extracts months from the date vector
yr <- years(d)

# restrict to NE Pacific
SLP <- SLP[,lat %in% 20:64 & lon >= 150]

# clean up x, y, lat, lon
y <- y[y %in% 20:64]
x <- x[x >= 150]

# reset lat/lon
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   


# and plot 
SLP.mean <- colMeans(SLP)
z <- t(matrix(SLP.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col)
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# looks good!

# remove seasonal means
f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(SLP, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 

xtra <- 12*((length(d)/12)-floor(length(d)/12)) # add trailing months

mu <- rbind(mu, mu[1:xtra,])
anom <- SLP - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}

# get a vector of weights (square root of the cosine of latitude)
weight <- sqrt(cos(lat*pi/180))

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

eig.1.1 <- EOF.1$U[,1] # reversing sign to match full time series pattern!
eig.2.1 <- EOF.1$U[,2]

eig.1.2 <- EOF.2$U[,1] # reversing sign to match full time series pattern!
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

png("./figs/slp_EOF.png", 4, 10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# Full time series!
z  <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1949-2021 (", var.all[1], "%)", sep=""), cex=0.8)

z  <- eig.2.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1949-2021 (", var.all[2], "%)", sep=""), cex=0.8)

#################
# Era 1!
z  <- eig.1.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1949-1968 (", var.1[1], "%)", sep=""), cex=0.8)

z  <- eig.2.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1949-1968 (", var.1[2], "%)", sep=""), cex=0.8)

#############
# Era 2
z  <- eig.1.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1969-1988 (", var.2[1], "%)", sep=""), cex=0.8)

z  <- eig.2.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1969-1988 (", var.2[2], "%)", sep=""), cex=0.8)

#######
# Era 3
z  <- eig.1.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1989-2008 (", var.3[1], "%)", sep=""), cex=0.8)

z  <- eig.2.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1989-2008 (", var.3[2], "%)", sep=""), cex=0.8)

#############
# Era 4
z  <- eig.1.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 2009-2021 (", var.4[1], "%)", sep=""), cex=0.8)

z  <- eig.2.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
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
  facet_wrap(~era, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-3.5,3.5)) +
  ggtitle("Mode 1")

m2 <- ggplot(filter(plot.diff, mode == "Mode 2"), aes(difference)) +
  geom_density(fill = "grey", lty = 0) +
  facet_wrap(~era, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-3.5,3.5)) +
  ggtitle("Mode 2")

png("./figs/slp_EOF_loading_differences_by_era.png", width=4, height=5.5, units='in', res=300)
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

ggsave("./figs/SD_SLP_loading_differences_by_era.png", width = 4.5, height = 5, units = 'in')


# sst EOFs
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

m <- as.numeric(m)
m.FMA <- m[m %in% 2:4]

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


## set up dataframe to correlate with PDO and NPGO ------------------------
pdo <- read.csv("./data/PDO.csv", skip = 1)

str(pdo)

pdo <- pdo %>%
  rename(pdo = Value) %>%
  mutate(year = floor(Date/100)) %>%
  mutate(m = Date-year*100) %>%
  select(-Date)

npgo <- read.csv("./data/NPGO.csv")

str(npgo)

npgo <- npgo %>% 
  rename(npgo = NPGO, m = MONTH, year = YEAR) 


pcs <- data.frame(year = as.numeric(as.character(yr.FMA)),
                  m = m.FMA,
                  pc1.all = scale(anom.detr.FMA %*% EOF.all$U[,1]),
                  pc2.all = scale(anom.detr.FMA %*% EOF.all$U[,2]),
                  pc1.era = c(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% -EOF.1$U[,1]),
                              scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% -EOF.2$U[,1]),
                              scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% EOF.3$U[,1]),
                              scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% EOF.4$U[,1])),
                  pc2.era = c(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% -EOF.1$U[,2]),
                              scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% -EOF.2$U[,2]),
                              scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% -EOF.3$U[,2]),
                              scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% -EOF.4$U[,2])))

pcs <- left_join(pcs, pdo)
pcs <- left_join(pcs, npgo)

## and plot------------------------------------------------

# set the limit for plotting 
lim.1 <- range(eig.1.all, eig.1.1, eig.1.2, eig.1.3, eig.1.4)
lim.2 <- range(-eig.2.all, -eig.2.1, -eig.2.2, -eig.2.3, -eig.2.4)

lim.1
lim.2

png("./figs/sst_EOF.png", 4, 10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2
cex = 0.7

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# Full time series!
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1949-2021 (", var.all[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.all, pcs$pdo), 2), sep = ""), cex = cex)

z <- rep(NA, ncol(SST))
z[!land] <- -eig.2.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1949-2021 (", var.all[2], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (NPGO) = ", round(cor(pcs$pc2.all, pcs$npgo, use = "p"), 2), sep = ""), cex = cex)

#################
# Era 1!
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1949-1968 (", var.1[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.era[pcs$year %in% 1949:1968], pcs$pdo[pcs$year %in% 1949:1968]), 2), sep = ""), cex = cex)

z <- rep(NA, ncol(SST))
z[!land] <- -eig.2.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1949-1969 (", var.1[2], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (NPGO) = ", round(cor(pcs$pc2.era[pcs$year %in% 1949:1968], pcs$npgo[pcs$year %in% 1949:1968], use = "p"), 2), sep = ""), cex = cex)


#############
# Era 2
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1969-1988 (", var.2[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.era[pcs$year %in% 1969:1988], pcs$pdo[pcs$year %in% 1969:1988]), 2), sep = ""), cex = cex)

z <- rep(NA, ncol(SST))
z[!land] <- -eig.2.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1969-1988 (", var.2[2], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (NPGO) = ", round(cor(pcs$pc2.era[pcs$year %in% 1969:1988], pcs$npgo[pcs$year %in% 1969:1988]), 2), sep = ""), cex = cex)

#######
# Era 3
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1989-2008 (", var.3[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.era[pcs$year %in% 1989:2008], pcs$pdo[pcs$year %in% 1989:2008]), 2), sep = ""), cex = cex)

z <- rep(NA, ncol(SST))
z[!land] <- -eig.2.3
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 1989-2008 (", var.3[2], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (NPGO) = ", round(cor(pcs$pc2.era[pcs$year %in% 1989:2008], pcs$npgo[pcs$year %in% 1989:2008]), 2), sep = ""), cex = cex)

#############
# Era 4
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 2009-2021 (", var.4[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.era[pcs$year %in% 2009:2021], pcs$pdo[pcs$year %in% 2009:2021]), 2), sep = ""), cex = cex)

z <- rep(NA, ncol(SST))
z[!land] <- -eig.2.4
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.2[2], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF2 2009-2021 (", var.4[2], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (NPGO) = ", round(cor(pcs$pc2.era[pcs$year %in% 2009:2021], pcs$npgo[pcs$year %in% 2009:2021]), 2), sep = ""), cex = cex)

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

## calculate PC1/2 for the full period, project across all eras, then calculate correlations with actual PC1/2 for each era -----------------------------------

pc1.all <- scale(anom.detr.FMA %*% EOF.all$U[,1])
pc2.all <- scale(anom.detr.FMA %*% EOF.all$U[,2])

# combine save for comparison with PDO/NPGO
pcs <- data.frame(decimal.year <- as.numeric(as.character(yr.FMA)) + (m.FMA - 0.5)/12,
                  pc1 = pc1.all,
                  pc2 = pc2.all)
write.csv(pcs, "./output/winter.sst.pcs.csv", row.names = F)

# calculate PC1-2 by projecting EOF for first era - call these "stationary" PCs
# the nomenclature is: stationary (st) / nonstationary (nst).mode.era.variable
st.pc1.1.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% EOF.all$U[,1]))
st.pc2.1.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% EOF.all$U[,2]))

st.pc1.2.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% EOF.all$U[,1]))
st.pc2.2.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% EOF.all$U[,2]))

st.pc1.3.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% EOF.all$U[,1]))
st.pc2.3.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% EOF.all$U[,2]))

st.pc1.4.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% EOF.all$U[,1]))
st.pc2.4.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% EOF.all$U[,2]))

# now the "nonstationary" PCs - i.e., those fit separately to each era
nst.pc1.1.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% -EOF.1$U[,1]))
nst.pc2.1.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1949:1968,] %*% EOF.1$U[,2]))

nst.pc1.2.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% -EOF.2$U[,1]))
nst.pc2.2.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1969:1988,] %*% EOF.2$U[,2]))

nst.pc1.3.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% EOF.3$U[,1]))
nst.pc2.3.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 1989:2008,] %*% EOF.3$U[,2]))

nst.pc1.4.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% EOF.4$U[,1]))
nst.pc2.4.sst <- as.vector(scale(anom.detr.FMA[yr.FMA %in% 2009:2021,] %*% EOF.4$U[,2]))

# combine into a data frame
# get # of total rows
length <- 3*length(nst.pc1.2.sst) + length(nst.pc1.4.sst)
st.nst.PCs <- data.frame(era = rep(c("1949-1968", "1969-1988", "1989-2008", "2009-2021"), times = c(length(nst.pc1.1.sst), length(nst.pc1.2.sst), length(nst.pc1.3.sst), length(nst.pc1.4.sst))),
                         mode = rep(c("PC1", "PC2"), each = length),
                         st.pc.score = c(st.pc1.1.sst, st.pc1.2.sst, st.pc1.3.sst, st.pc1.4.sst, st.pc2.1.sst, st.pc2.2.sst, st.pc2.3.sst, st.pc2.4.sst),
                         nst.pc.score = c(nst.pc1.1.sst, nst.pc1.2.sst, nst.pc1.3.sst, nst.pc1.4.sst, nst.pc2.1.sst, nst.pc2.2.sst, nst.pc2.3.sst, nst.pc2.4.sst))

cors <- data.frame(era = rep(c("1949-1968", "1969-1988", "1989-2008", "2009-2021"), each = 2),
                   mode = rep(c("PC1", "PC2"), times = 4),
                   cor = c(cor(st.pc1.1.sst, nst.pc1.1.sst),
                           cor(st.pc2.1.sst, nst.pc2.1.sst),
                           cor(st.pc1.2.sst, nst.pc1.2.sst),
                           cor(st.pc2.2.sst, nst.pc2.2.sst),
                           cor(st.pc1.3.sst, nst.pc1.3.sst),
                           cor(st.pc2.3.sst, nst.pc2.3.sst),
                           cor(st.pc1.4.sst, nst.pc1.4.sst),
                           cor(st.pc2.4.sst, nst.pc2.4.sst)))

cors <- cors %>%
  arrange(mode)
cors

ggplot(st.nst.PCs, aes(st.pc.score, nst.pc.score)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_grid(era~mode)

cc <- data.frame(era = filter(st.nst.PCs, mode == "PC1")$era,
                 PC1 = filter(st.nst.PCs, mode == "PC1")$st.pc.score,
                 PC2 = filter(st.nst.PCs, mode == "PC2")$st.pc.score)

ggplot(cc, aes(PC1, PC2)) +
  geom_point() +
  facet_wrap(~era)


cc1 <- cc %>%
  filter(era == "1949-1968")

cor1 <- round(cor(cc1$PC1, cc1$PC2), 2) # -0.24

cc1 <- cc %>%
  filter(era == "1969-1988")

cor2 <- round(cor(cc1$PC1, cc1$PC2), 2) # -0.42

cc1 <- cc %>%
  filter(era == "1989-2008")

cor3 <- round(cor(cc1$PC1, cc1$PC2), 2) # 0.54

cc1 <- cc %>%
  filter(era == "2009-2021")

cor4 <-  round(cor(cc1$PC1, cc1$PC2), 2) # 0.13

# annotate
anno <- data.frame(x1 = -2,
                   y1 = 2,
                   labs = c(paste("r = ", cor1, sep = ""),
                            paste("r = ", cor2, sep = ""),
                            paste("r = ", cor3, sep = ""),
                            paste("r = ", cor4, sep = "")),
                   era = c("1949-1968", "1969-1988", "1989-2008", "2009-2021"))

ggplot(cc, aes(PC1, PC2)) +
  geom_point() +
  facet_wrap(~era) +
  geom_text(data = anno, aes(x = x1, y = y1, label = labs))

### now correlate with PDO/NPGO -----------------------------------------
pdo <- read.csv("./data/PDO.csv", skip = 1)

str(pdo)

pdo <- pdo %>%
  rename(pdo = Value) %>%
  mutate(year = floor(Date/100)) %>%
  mutate(m = Date-year*100) %>%
  select(-Date)

npgo <- read.csv("./data/NPGO.csv")

str(npgo)

npgo <- npgo %>% 
  rename(npgo = NPGO, m = MONTH, year = YEAR) 

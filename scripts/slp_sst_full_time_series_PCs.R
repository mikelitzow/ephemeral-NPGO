# create PC1-2 time series for SLPa and SSTa for integration
 
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

### load and process SLP data -------------------
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

# EOF
# weighting the columns
EOF.slp <- svd.triplet(cov(anom.detr), col.w=weight)

# calculate PC1-2
pc1.slp <- anom.detr %*% EOF.slp$U[,1]
pc2.slp <- anom.detr %*% EOF.all$U[,2]

# give month and year slp-specific names for later...
m.slp <- as.numeric(m)
yr.slp <- as.numeric(as.character(yr))

### load and process SST -----------------------
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

# limit to 1948 - 2021 (to match SLP data)
m <- m[yr %in% 1948:2021]
SST <- SST[yr %in% 1948:2021,]
d <- d[yr %in% 1948:2021]
yr <- yr[yr %in% 1948:2021]

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

# EOF
# weighting the columns
EOF.sst <- svd.triplet(cov(anom.detr), col.w=weight)

# calculate PC1-2
pc1.sst <- anom.detr %*% EOF.sst$U[,1]
pc2.sst <- anom.detr %*% EOF.sst$U[,2]

# check that dates are identical for slp & sst
identical(m.slp, as.numeric(m))
identical(yr.slp, as.numeric(as.character(yr)))

# they are

# combine into a data frame;
# scale all PCs to remove units
sst.pcs <- data.frame(year = yr.slp, 
                      m = m.slp,
                      month = m,
                      variable = "sst",
                      pc1 = scale(pc1.sst),
                      pc2 = scale(pc2.sst))

slp.pcs <- data.frame(year = yr.slp, 
                      m = m.slp,
                      month = m,
                      variable = "slp",
                      pc1 = scale(pc1.slp),
                      pc2 = scale(pc2.slp))

all.pcs <- rbind(sst.pcs, slp.pcs)

# plot to check
plot.dat <- all.pcs %>%
  pivot_longer(cols = c(-m, -month, -year, -variable)) %>%
  mutate(decimal.year = year + (m-0.5)/12)

ggplot(plot.dat, aes(decimal.year, value)) +
  geom_line() +
  facet_grid(name~variable)

# looks good!

# save
write.csv(all.pcs, "./output/slp_sst_PCs_1948-2021.csv", row.names = F)
  
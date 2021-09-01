# test for era-dependent changes in correlations between SST EOFs and upwelling indices

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

# EOF by era
# weighting the columns
EOF.1 <- svd.triplet(cov(anom.detr[yr %in% 1949:1968,]), col.w=weight)
EOF.2 <- svd.triplet(cov(anom.detr[yr %in% 1969:1988,]), col.w=weight)
EOF.3 <- svd.triplet(cov(anom.detr[yr %in% 1989:2008,]), col.w=weight)
EOF.4 <- svd.triplet(cov(anom.detr[yr %in% 2009:2021,]), col.w=weight)

# calculate PC1-2
pc1.1 <- anom.detr[yr %in% 1949:1968,] %*% EOF.1$U[,1]
pc2.1 <- anom.detr[yr %in% 1949:1968,] %*% EOF.1$U[,2]

pc1.2 <- anom.detr[yr %in% 1969:1988,] %*% EOF.2$U[,1]
pc2.2 <- anom.detr[yr %in% 1969:1988,] %*% EOF.2$U[,2]

pc1.3 <- anom.detr[yr %in% 1989:2008,] %*% EOF.1$U[,1]
pc2.3 <- anom.detr[yr %in% 1989:2008,] %*% EOF.1$U[,2]

pc1.4 <- anom.detr[yr %in% 2009:2021,] %*% EOF.2$U[,1]
pc2.4 <- anom.detr[yr %in% 2009:2021,] %*% EOF.2$U[,2]

# combine into dfs
decimal.yr <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12

era1 <- data.frame(decimal.yr = decimal.yr[yr >= 1949 & yr < 1969],
                   pc1 = pc1.1,
                   pc2 = pc2.1)

era2 <- data.frame(decimal.yr = decimal.yr[yr >= 1969 & yr < 1989],
                   pc1 = pc1.2,
                   pc2 = pc2.2)

era3 <- data.frame(decimal.yr = decimal.yr[yr >= 1989 & yr < 2009],
                   pc1 = pc1.3,
                   pc2 = pc2.3)

era4 <- data.frame(decimal.yr = decimal.yr[yr >= 2009],
                   pc1 = pc1.4,
                   pc2 = pc2.4)


pcs <- rbind(era1, era2, era3, era4)


## correlations with Bakun indices--------------
bakun <- read.csv("./data/bakun.csv")
head(bakun)
tail(bakun)

# organize
unique(bakun$LAT)

# limit to sites associated with PDO/NPGO

bakun$mode <- if_else(bakun$LAT %in% c("30N", "31N", "33N", "36N"), "NPGO",
                      if_else(bakun$LAT %in% c("39N", "42N", "45N", "48N"), "PDO", "NA"))

# change months to #s for ease of processing!
names(bakun)[4:15] <- 1:12

bakun <- bakun %>%
  filter(mode %in% c("PDO", "NPGO")) %>%
  pivot_longer(cols = c(-LAT, -LONG, -YEAR, -mode)) %>%
  mutate(name = as.integer(name)) %>%
  rename(year = YEAR, month = name) %>%
  mutate(decimal.yr = year + (month-0.5)/12) %>%
  group_by(mode, decimal.yr) %>%
  summarise(anom = mean(value))

# plot

ggplot(bakun, aes(decimal.yr, anom)) +
  geom_line() + 
  facet_wrap(~mode)

bakun <- left_join(bakun, pcs) %>%
  pivot_longer(cols = c(pc1, pc2))

bakun$era <- case_when(
  bakun$decimal.yr < 1969 ~ "1949-1968",
  bakun$decimal.yr > 1969 & bakun$decimal.yr < 1989 ~ "1969-1988",
  bakun$decimal.yr > 1989 & bakun$decimal.yr < 2009 ~ "1989-2008",
  bakun$decimal.yr > 2009 ~ "2009-2021"
)

ggplot(bakun, aes(value, anom)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_grid(era~mode)

cors <- bakun %>%
  group_by(mode, name, era) %>%
  summarise(cor = cor(anom, value, use = "p")
    
  )

# try again with individual Bakun time series
bakun <- read.csv("./data/bakun.csv")
head(bakun)
tail(bakun)

# organize
unique(bakun$LAT)

# limit to sites associated with PDO/NPGO
bakun <- bakun %>%
  filter(!is.na(LAT))

# change months to #s for ease of processing!
names(bakun)[4:15] <- 1:12

bakun <- bakun %>%
  pivot_longer(cols = c(-LAT, -LONG, -YEAR)) %>%
  mutate(name = as.integer(name)) %>%
  rename(year = YEAR, month = name) %>%
  mutate(decimal.yr = year + (month-0.5)/12) 

bakun <- left_join(bakun, pcs)

bakun$era <- case_when(
  bakun$decimal.yr < 1969 ~ "1949-1968",
  bakun$decimal.yr > 1969 & bakun$decimal.yr < 1989 ~ "1969-1988",
  bakun$decimal.yr > 1989 & bakun$decimal.yr < 2009 ~ "1989-2008",
  bakun$decimal.yr > 2009 ~ "2009-2021"
)

cors <- bakun %>%
  group_by(LAT, era) %>%
  summarise(corPC1 = cor(pc1, value, use = "p"),
            corPC2 = cor(pc2, value, use = "p"))



# restrict to the area relevant to each mode from Di Lorenzo et al. 2008

cors <- cors %>%
  filter(LAT != "31N") 

cors <- cors %>%
  filter(LAT != "60N") 

cors <- cors %>%
  filter(LAT != "999") 

cor_long <- cors %>%
  pivot_longer(cols = c(corPC1, corPC2))

ggplot(cor_long, aes(era, value)) +
  geom_bar(stat = "identity") +
  facet_grid(name~LAT)

sd <- cors %>%
  group_by(LAT, era) %>%
  summarise(sdPC1 = sd(corPC1                                                                                        ),
            sdPC2 = sd(corPC2))

sd  

NPGOsd <- NPGO %>%
  group_by(LAT) %>%
  summarise(sd = sd(corPC2))
NPGOsd  

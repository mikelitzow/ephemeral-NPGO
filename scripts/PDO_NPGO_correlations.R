# calculate correlations between PDO/NPGO and PC1/2 from our calculations

library(tidyverse)

# load data

pcs <- read.csv("./output/winter.sst.pcs.csv")

str(pcs)

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

dat <- left_join(pcs, pdo)
dat <- left_join(dat, npgo)

ggplot(dat, aes(pc1, pdo)) +
  geom_point()


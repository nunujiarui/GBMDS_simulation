# This script contains codes for real data application
# This example is for section 5.2 Geographical data

library(vegan)
library(tidyverse)
library(tibble)
library(mvtnorm)
library(MCMCpack)
library(MASS)
library(truncnorm)
library(reshape2)
library(parallel)
library(gridExtra)
library(grid)
library(sn)
library(fGarch)
library(ggrepel)
library(gridExtra)
library(philentropy)

set.seed(10)

## source ASMC models
# helper functions
source(file = "R/ASMC_helper_fun.R")
# ASMC function
source(file = "R/ASMC_fun.R")
# truncated Normal
source(file = "R/ASMC_truncatedN.R")
# truncated T
source(file = "R/ASMC_truncatedT.R")
# truncated skewed Normal
source(file = "R/ASMC_truncatedSkewedN.R")
# function to plot ASMC result
source(file = "R/ASMC_plot.R")

dist.metric <- c("manhattan", "euclidean", "chebyshev")
# change this index to try different distance metrics
dist.metric.index <- 3

## read in data
data.all <- read.csv(file = "data/us-cities-demographics.csv", sep = ";")

data <- data.all %>%
  distinct(City, State, .keep_all = TRUE) %>%
  dplyr::select(City, Median.Age, Number.of.Veterans,
                Foreign.born, Average.Household.Size)

names(data)[1] <- "city"
data$city <- as.character(data$city)

us.cities.all <- read.csv(file = "data/uscities.csv")

us.cities.all.distinct <- distinct(us.cities.all, city, state_id,
                                   .keep_all = TRUE)

us.cities.all.distinct$city <- as.character(us.cities.all.distinct$city)

us.cities.org <- us.cities.all.distinct %>%
  top_n(15, population) %>%
  dplyr::select(city, lat, lng, population, density) %>%
  left_join(data, by = "city") %>%
  tibble::column_to_rownames(var = "city")
# us.cities contains information about
# "city", "lat", "lng", "population", "density", "Median.Age", "Number.of.Veterans", "Foreign.born", "Average.Household.Size"

us.cities.loc <- us.cities.org %>%
  dplyr::select(lng, lat)

## add some noises in the data
n.lat.noise <- 5
n.lng.noise <- 5

lat.sd <- sd(us.cities.org$lat)
lng.sd <- sd(us.cities.org$lng)

lat.noise <- replicate(n.lat.noise, rnorm(n = nrow(us.cities.org), mean = 0, sd = lat.sd))
colnames(lat.noise) <- paste0("lat.noise.", seq(1, n.lat.noise))

lng.noise <- replicate(n.lng.noise, rnorm(n = nrow(us.cities.org), mean = 0, sd = lng.sd))
colnames(lng.noise) <- paste0("lng.noise.", seq(1, n.lng.noise))

#data <- cbind(us.cities.org, lat.noise, lng.noise)

data <- cbind(us.cities.org[,c(1,2)], lat.noise, lng.noise)

## calculate observed distance from data with different set of weights
data <- scale(data)

# set weights
#aa <- 1/6   # this corresponds to the case when R = 1  (equal weight)
aa <- 4/9   # this corresponds to the case when R = 4  (low signal-to-noise ratio)
#aa <- 10/15 # this corresponds to the case when R = 10 (high signal-to-noise ratio)
weight <- c(rep(aa/2, 2), rep((1-aa)/10, 10))

data.weight <- t(apply(data, 1, function(x) sqrt(weight) * x))
dis <- philentropy::distance(data.weight, method = dist.metric[dist.metric.index])
colnames(dis) =  rownames(dis) <- rownames(data.weight)

dis.df <- data.frame(dis = dis[upper.tri(dis)])
ggplot(dis.df, aes(x = dis)) +
  geom_histogram(bins = 40) +
  labs(x = paste0(dist.metric[dist.metric.index], " dissimilarity")) +
  theme_bw()

## general settings
p <- 2
n <- nrow(dis)

cmds.result <- cmdscale(d = dis, k = p,
                        eig = TRUE, add = FALSE, x.ret = FALSE)
class(cmds.result) <- append(class(cmds.result), "CMDS")


## set hyperparameters
sim.a.initial <- 5
SSR.initial <- SSRFun(d.mat = dis, delta.mat = as.matrix(dist(cmds.result$points)))
sim.m <- n * (n - 1)/2
sim.b.initial <- SSR.initial/sim.m

sim.alpha.initial <- 1/2
sample.cov <- cov(cmds.result$points)
sim.beta.initial <- (1/2)*diag(sample.cov)

df.initial <- 5
c.initial <- -2
d.initial <- 2

constant.multiple <- 2.38^2/5

hyperparList <- list(a = sim.a.initial, b = sim.b.initial,
                     alpha = sim.alpha.initial, beta = sim.beta.initial,
                     df = df.initial,
                     c = c.initial, d = d.initial, constant.multiple = constant.multiple)

reference.x.sd <- diag(rep(0.01, p))

tuningparList <- list(K = 200, phi = 0.8, eps = 0.5)
n.core <- detectCores() - 1

## run GBMDS-ASMC
#model <- truncatedN(hyperparList, p, reference.x.sd)
#model <- truncatedT(hyperparList, p, reference.x.sd)
model <- truncatedSkewedN(hyperparList, p, reference.x.sd)

start.time <- Sys.time()
asmc.result <- ASMC(model = model,
                    dist.mat = dis,
                    tuningparList, n.core, cmds.result = cmds.result$points,
                    metric = dist.metric[dist.metric.index])
end.time <- Sys.time()
end.time - start.time


# log marginal likelihood estimate
asmc.result$logZ

# posterior inference
index.asmc <- which.min(asmc.result$SSR.output)
asmc.res <- asmc.result$xi.output[[index.asmc]]

## compute STRESS values
cmds.stress <- stressFun(d.mat = dis,
                         delta.mat = philentropy::distance(cmds.result$points, method = "euclidean"))
asmc.stress <- stressFun(d.mat = dis,
                         delta.mat = philentropy::distance(asmc.res, method = "euclidean"))

cat("STRESS value for CMDS is", round(cmds.stress, 4))
cat("STRESS value for BMDS with ASMC is", round(asmc.stress, 4))

# scatter plot in 2D
grid.arrange(plot(cmds.result, data.label = rownames(data),
                  plot.title = "cmds",
                  rotate = T, rotate.angle = pi/3,
                  flip = F, flip.h = T, flip.v = F),
             plot(asmc.result, data.label = rownames(data), ci.level = 0.95,
                  plot.title = "asmc",
                  procruste.transform = TRUE,
                  rotate = T, rotate.angle = pi/3,
                  flip = F, flip.h = T, flip.v = F), ncol=2)



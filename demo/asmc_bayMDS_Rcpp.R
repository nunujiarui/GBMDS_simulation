library(bayMDS)


data(cityDIST) # 30 * 30

start.time <- Sys.time()
out <- bmds(cityDIST,min_p=2,max_p=2) 
end.time <- Sys.time()
end.time - start.time # 4 sec
#if(interactive()){bayMDSApp(out)}


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
library(R.matlab)
library(Rcpp)

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
# source(file = "R/ASMC_truncatedSkewedN.R")
source(file = "R/ASMC_truncatedSkewedN_Rcpp.R")
# function to plot ASMC result
source(file = "R/ASMC_plot.R")

source(file = "Rcpp/cppFun.R")
Rcpp::sourceCpp(file = "Rcpp/mvnorm.cpp")

## choose dissimilarity metric
dist.metric <- c("manhattan", "euclidean", "chebyshev")
# change this index to try different distance metrics
dist.metric.index <- 2



## create containers to store the results
tn.logZ = tn.stress = tsn.logZ = tsn.stress = tt.logZ = tt.stress <- numeric()

i = 1

dis <- cityDIST
  
## MDS with annealed SMC
p = 2 # dims in lower space

# use cmds to get initial values
cmds.result <- cmdscale(d = dis, k = p,
                        eig = TRUE, add = FALSE, x.ret = FALSE)
class(cmds.result) <- append(class(cmds.result), "CMDS")

n <- nrow(dis)

# hyperparameters
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

tuningparList <- list(K = 200, phi = 0.80, eps = 0.5)
n.core <- detectCores() - 2

## 1. annealed SMC with truncated Normal
# model <- truncatedN(hyperparList, p, reference.x.sd)
# start.time <- Sys.time()
# #Rprof("/Users/erin/Desktop/ASMC.out")
# asmc.result1 <- ASMC(model = model,
#                      dist.mat = dis,
#                      tuningparList, n.core, cmds.result = cmds.result$points,
#                      metric = dist.metric[dist.metric.index])
# #Rprof(NULL)
# end.time <- Sys.time()
# end.time - start.time 
# 
# #summaryRprof(filename = "/Users/erin/Desktop/ASMC.out")
# 
# 
# tn.logZ[i] = asmc.result1$logZ
# 
# # posterior inference
# tn.index.asmc <- which.min(asmc.result1$SSR.output)
# tn.asmc.res <- asmc.result1$xi.output[[tn.index.asmc]]
# 
# tn.stress[i] <- round(stressFun(d.mat = dis,
#                                 delta.mat = as.matrix(dist(tn.asmc.res, method = "euclidean"))), 4)


## 2. anneal SMC with truncated skewed Normal
model <- truncatedSkewedN(hyperparList, p, reference.x.sd)

# run annealed SMC
set.seed(1055)
start.time <- Sys.time()
asmc.result2 <- ASMC(model = model,
                     dist.mat = dis,
                     tuningparList, n.core, cmds.result = cmds.result$points,
                     metric = dist.metric[dist.metric.index])
end.time <- Sys.time()
end.time - start.time  

tsn.logZ[i] = asmc.result2$logZ

# posterior inference
tsn.index.asmc <- which.min(asmc.result2$SSR.output)
tsn.asmc.res <- asmc.result2$xi.output[[tsn.index.asmc]]

tsn.stress[i] <- round(stressFun(d.mat = dis,
                                 delta.mat = as.matrix(dist(tsn.asmc.res, method = "euclidean"))), 4)

# ## 3. anneal SMC with truncated T
# model <- truncatedT(hyperparList, p, reference.x.sd)
# 
# start.time <- Sys.time()
# asmc.result3 <- ASMC(model = model,
#                      dist.mat = dis,
#                      tuningparList, n.core, cmds.result = cmds.result$points,
#                      metric = dist.metric[dist.metric.index])
# end.time <- Sys.time()
# end.time - start.time
# 
# tt.logZ[i] = asmc.result3$logZ
# 
# # posterior inference
# tt.index.asmc <- which.min(asmc.result3$SSR.output)
# tt.asmc.res <- asmc.result3$xi.output[[tt.index.asmc]]
# 
# tt.stress[i] <- round(stressFun(d.mat = dis,
#                                 delta.mat = as.matrix(dist(tt.asmc.res, method = "euclidean"))), 4)

##### Scotch Whisky Dataset #####

# https://outreach.mathstat.strath.ac.uk/outreach/nessie/nessie_whisky.html
# https://www.kaggle.com/datasets/koki25ando/scotch-whisky-dataset/data

library(R.matlab)

data <- readMat("data/whisky.mat")

dis <- data$whiskydist # 86 * 86

start.time <- Sys.time()
out <- bmds(dis, min_p=2, max_p=2, niter = 5000) 
end.time <- Sys.time()
end.time - start.time # 14 sec


## MDS with annealed SMC
p = 2 # dims in lower space

# use cmds to get initial values
cmds.result <- cmdscale(d = dis, k = p,
                        eig = TRUE, add = FALSE, x.ret = FALSE)
class(cmds.result) <- append(class(cmds.result), "CMDS")

n <- nrow(dis)

# hyperparameters
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

tuningparList <- list(K = 200, phi = 0.80, eps = 0.5)
n.core <- detectCores() - 0

## 2. anneal SMC with truncated skewed Normal
model <- truncatedSkewedN(hyperparList, p, reference.x.sd)

# run annealed SMC
start.time <- Sys.time()
asmc.result2 <- ASMC(model = model,
                     dist.mat = dis,
                     tuningparList, n.core, cmds.result = cmds.result$points,
                     metric = dist.metric[dist.metric.index])
end.time <- Sys.time()
end.time - start.time # 41 sec

tsn.logZ = asmc.result2$logZ

# posterior inference
tsn.index.asmc <- which.min(asmc.result2$SSR.output)
tsn.asmc.res <- asmc.result2$xi.output[[tsn.index.asmc]]

tsn.stress <- round(stressFun(d.mat = dis,
                              delta.mat = as.matrix(dist(tsn.asmc.res, method = "euclidean"))), 10)





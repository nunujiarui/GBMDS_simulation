# This script contains codes for simulation study
# This example is for section 5.1 NIPS text data with incremental dimensions

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
library(text2vec)
library(smacof)

## source ASMC models
# helper functions
source(file = "R/ASMC_helper_fun.R")
# ASMC function
source(file = "R/ASMC_fun.R")
source(file = "R/ASMC_incr_fun.R")
# truncated Normal
source(file = "R/ASMC_truncatedN.R")
source(file = "R/ASMC_truncatedN_incr.R")
# truncated T
source(file = "R/ASMC_truncatedT.R")
source(file = "R/ASMC_truncatedT_incr.R")
# truncated skewed Normal
source(file = "R/ASMC_truncatedSkewedN.R")
source(file = "R/ASMC_truncatedSkewedN_incr.R")

## choose dissimilarity metric
dist.metric <- c("manhattan", "euclidean", "chebyshev", "cosine")
# change this index to try different distance metrics
dist.metric.index <- 4

## read in data
data <- readRDS(file = "data/NIPS_data_small.rds")

# calculate observed distances
if (dist.metric.index == 4){
  dis <- 1-philentropy::distance(t(data), method = dist.metric[dist.metric.index], mute.message = TRUE)
} else{
  dis <- philentropy::distance(t(data), method = dist.metric[dist.metric.index], mute.message = TRUE)
}


hist(dis[upper.tri(dis)],
     xlab = "dissimilarity",
     main = paste0("Histogram of ", dist.metric[dist.metric.index], " dissimilarity"))

## set incremental size
n.t1 <- 50      # correspond to n
incr.step <- 5  # correspond to m

## create containers to store the results
time1 <- c(rep(NA, 20))
stress1 <- c(rep(NA, 20))
time2 <- c(rep(NA, 20))
stress2 <- c(rep(NA, 20))
time12 <- c(rep(NA, 20))
stress12 <- c(rep(NA, 20))

for (i in 1:20){

  ## select initial dis and incremental dis
  dis.t1 <- dis[1:n.t1, 1:n.t1]
  dis.t2 <- dis[1:(n.t1+incr.step), 1:(n.t1+incr.step)]

  ## MDS with annealed SMC
  p = 2 # dims in lower space

  # use cmds to get initial values
  t1.cmds.res <- cmdscale(d = dis.t1, k = p,
                          eig = FALSE, add = FALSE, x.ret = FALSE)

  t2.cmds.res <- cmdscale(d = dis.t2, k = p,
                          eig = FALSE, add = FALSE, x.ret = FALSE)

  ## set hyperparameters
  constant.multiple <- 2.38^2

  n.t1 <- nrow(dis.t1)

  sim.a.initial.t1 <- 5
  SSR.initial.t1 <- SSRFun(d.mat = dis.t1, delta.mat = as.matrix(dist(t1.cmds.res)))
  sim.m.t1 <- n.t1 * (n.t1 - 1)/2
  sim.b.initial.t1 <- SSR.initial.t1/sim.m.t1

  sim.alpha.initial.t1 <- 1/2
  sample.cov.t1 <- cov(t1.cmds.res)
  sim.beta.initial.t1 <- (1/2)*diag(sample.cov.t1)

  df.initial.t1 <- 5
  c.initial.t1 <- -2
  d.initial.t1 <- 2

  hyperpars.t1 <- list(a = sim.a.initial.t1, b = sim.b.initial.t1,
                       alpha = sim.alpha.initial.t1, beta = sim.beta.initial.t1, df = df.initial.t1,
                       c = c.initial.t1, d = d.initial.t1, constant.multiple = constant.multiple)

  n.t2 <- nrow(dis.t2)

  sim.a.initial.t2 <- 5
  SSR.initial.t2 <- SSRFun(d.mat = dis.t2, delta.mat = as.matrix(dist(t2.cmds.res)))
  sim.m.t2 <- n.t2 * (n.t2 - 1)/2
  sim.b.initial.t2 <- SSR.initial.t2/sim.m.t2

  sim.alpha.initial.t2 <- 1/2
  sample.cov.t2 <- cov(t2.cmds.res)
  sim.beta.initial.t2 <- (1/2)*diag(sample.cov.t2)

  df.initial.t2 <- 5
  c.initial.t2 <- -2
  d.initial.t2 <- 2

  hyperpars.t2 <- list(a = sim.a.initial.t2, b = sim.b.initial.t2,
                       alpha = sim.alpha.initial.t2, beta = sim.beta.initial.t2, df = df.initial.t2,
                       c = c.initial.t2, d = d.initial.t2, constant.multiple = constant.multiple)

  reference.x.sd <- diag(rep(0.01, p))


  tuningparList1 <- list(K = 100, phi = 0.80, eps = 0.5)
  n.core <- detectCores() - 2

  ## set 1 (original part) only
  model <- truncatedT(hyperparList = hyperpars.t1, p, reference.x.sd)

  n <- nrow(dis.t1)
  start.time <- Sys.time()
  t1.asmc.result <- ASMC(model = model,
                         dist.mat = dis.t1,
                         tuningparList = tuningparList1, n.core, cmds.result = t1.cmds.res,
                         metric = dist.metric[dist.metric.index])
  end.time <- Sys.time()
  time1[i] <- as.numeric(end.time - start.time, units = "secs")

  # posterior inference
  index.asmc.t1 <- which.min(t1.asmc.result$SSR.output)
  t1.asmc.res <- t1.asmc.result$xi.output[[index.asmc.t1]]

  stress1[i] <- stressFun(d.mat = dis.t1,
                          delta.mat = 1-philentropy::distance(t1.asmc.res,
                                                            method = dist.metric[dist.metric.index], mute.message = TRUE))

  # store results
  t1.asmc.res.all <- list(xi = t1.asmc.res, sigma2 = t1.asmc.result$sigma2.output[[index.asmc.t1]],
                          lambda = t1.asmc.result$lambda.output[[index.asmc.t1]],
                          psi = t1.asmc.result$psi.output[[index.asmc.t1]],
                          g = t1.asmc.result$g.output[[index.asmc.t1]])

  ## set 2 (original + incremental part) only
  n <- nrow(dis.t2)
  start.time <- Sys.time()
  t2.asmc.result <- ASMC(model = model,
                         dist.mat = dis.t2,
                         tuningparList = tuningparList1, n.core, cmds.result = t2.cmds.res,
                         metric = dist.metric[dist.metric.index])
  end.time <- Sys.time()
  time2[i] <- as.numeric(end.time - start.time, units = "secs")

  # posterior inference
  index.asmc.t2 <- which.min(t2.asmc.result$SSR.output)
  t2.asmc.res <- t2.asmc.result$xi.output[[index.asmc.t2]]

  stress2[i] <- stressFun(d.mat = dis.t2,
                          delta.mat = 1-philentropy::distance(t2.asmc.res,
                                                              method = dist.metric[dist.metric.index], mute.message = TRUE))

  ## set 1 and 2 combine
  reference.x.sd <- cov(t1.asmc.res.all$xi)
  model <- truncatedT_incr(hyperparList = hyperpars.t2, p, reference.x.sd)

  n.incr <- nrow(dis.t2) - nrow(dis.t1)
  tuningparList2 <- list(K = 50, # reduce the particle size to half
                         phi = 0.80, eps = 0.5)

  start.time <- Sys.time()
  t12.asmc.result <- ASMC_incr(model = model,
                               dist.mat = dis.t2,
                               tuningparList =  tuningparList2, n.core, prev.result = t1.asmc.res.all,
                               metric = dist.metric[dist.metric.index])
  end.time <- Sys.time()
  time12[i] <- as.numeric(end.time - start.time, units = "secs")

  # posterior inference
  index.asmc.t12 <- which.min(t12.asmc.result$SSR.output)
  t12.asmc.res <- t12.asmc.result$xi.output[[index.asmc.t12]]

  stress12[i] <- stressFun(d.mat = dis.t2,
                           delta.mat = 1-philentropy::distance(t12.asmc.res,
                                                               method = dist.metric[dist.metric.index], mute.message = TRUE))

}

# computation time
summary(time2)  # without adaptive inference
summary(time12) # with adaptive inference
# STRESS value
summary(stress2)  # without adaptive inference
summary(stress12) # with adaptive inference

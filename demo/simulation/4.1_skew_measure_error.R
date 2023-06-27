# This script contains codes for simulation study
# This example is for section 4.1 Experiment 1: Data with skewed errors

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

## create containers to store the results
tn.logZ = tn.stress = tsn.logZ = psi = tsn.stress = tt.logZ = tt.stress <- numeric()

## run simulations 20 times
for (i in 1:20){

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

  ## choose dissimilarity metric
  dist.metric <- c("manhattan", "euclidean", "chebyshev")
  # change this index to try different distance metrics
  dist.metric.index <- 2

  ## simulate data
  n.col = 20
  n.row = 50

  true.obs <- matrix(data = c(rnorm(n.col * n.row/2, mean = 0, sd = 1),
                              rnorm(n.col * n.row/4, 100, 10), rnorm(n.col * n.row/4, -10, 1)), nrow = n.row)
  true.dis <- as.matrix(dist(true.obs, method = dist.metric[dist.metric.index]))
  true.dis.df <- data.frame(dis = true.dis[upper.tri(true.dis)])

  ggplot(true.dis.df, aes(x = dis)) +
    geom_histogram(bins = 20) +
    labs(x = paste0(dist.metric[dist.metric.index], " dissimilarities of \n true measurements"),
         y = "Frequency") +
    theme_bw()

  summary(as.vector(true.obs))

  n.total <- n.col * n.row
  pct.outlier1 <- 0.2
  n.outlier1 <- round(n.total * pct.outlier2, 0)
  pct.outlier2 <- 0.02
  n.outlier2 <- round(n.total * pct.outlier2, 0)

  # add some minor measurement errors to all observations
  obs.obs <- true.obs + rnorm(n.total)
  # add relatively large measurement errors to some observations
  index.outlier1 <- sample(1:n.total, size = n.outlier1)
  obs.obs[index.outlier1] <- true.obs[index.outlier1] + rnorm(n.outlier1, mean = 10, sd = 1)
  # add extremely large measurement errors to few observations (considered as outliers)
  index.outlier2 <- sample(1:n.total, size = n.outlier2)
  obs.obs[index.outlier2] <- true.obs[index.outlier2] + rnorm(n.outlier2, mean = 20, sd = 1)


  # step <- 3
  # row.ind <- seq(1, 50, step)
  # col.ind <- seq(1, 20, step)
  # obs.obs[row.ind,] <- true.obs[row.ind,] + rnorm(length(row.ind), mean = 10, sd = 1)
  # obs.obs[,col.ind] <- true.obs[,col.ind] + rnorm(length(col.ind), mean = 10, sd = 1)
  # # add extremely large measurement errors to few observations (considered as outliers)
  # index.outlier <- sample(1:n.total, size = n.outlier)
  # obs.obs[index.outlier] <- true.obs[index.outlier] + rnorm(n.outlier, mean = 20, sd = 1)

  ## calculate the dissimilarities with simulated data
  obs.dis <- as.matrix(dist(obs.obs, method = dist.metric[dist.metric.index]))
  obs.dis.df <- data.frame(dis = obs.dis[upper.tri(obs.dis)])

  ggplot(obs.dis.df, aes(x = dis)) +
    geom_histogram(bins = 20) +
    labs(x = paste0(dist.metric[dist.metric.index], " dissimilarities of \n observed measurements"),
         y = "Frequency") +
    theme_bw()

  # get observed measurement errors with Gaussian + skewed error
  measure.error <- obs.dis - true.dis
  measure.error.df <- data.frame(dis = measure.error[upper.tri(measure.error)])

  ggplot(measure.error.df, aes(x = dis)) +
    geom_histogram(bins = 20) +
    labs(x = "measurement error", y = "Frequency") +
    theme_bw()

  ## MDS with annealed SMC
  p = 2 # dims in lower space

  # observed dissimilarity
  dis <- obs.dis

  # true dissimilarity
  dis.true <- true.dis

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

  tuningparList <- list(K = 200, # number of particles
                        phi = 0.80, eps = 0.5)
  n.core <- detectCores() - 2

  ## 1. annealed SMC with truncated Normal
  model <- truncatedN(hyperparList, p, reference.x.sd)
  start.time <- Sys.time()
  asmc.result1 <- ASMC(model = model,
                       dist.mat = dis,
                       tuningparList, n.core, cmds.result = cmds.result$points,
                       metric = dist.metric[dist.metric.index])
  end.time <- Sys.time()
  end.time - start.time

  tn.logZ[i] = asmc.result1$logZ

  # posterior inference
  tn.index.asmc <- which.min(asmc.result1$SSR.output)
  tn.asmc.res <- asmc.result1$xi.output[[tn.index.asmc]]

  tn.stress[i] <- round(stressFun(d.mat = dis.true,
                                  delta.mat = as.matrix(dist(tn.asmc.res, method = "euclidean"))), 4)


  ## 2. anneal SMC with truncated skewed Normal
  model <- truncatedSkewedN(hyperparList, p, reference.x.sd)

  # run annealed SMC
  start.time <- Sys.time()
  asmc.result2 <- ASMC(model = model,
                       dist.mat = dis,
                       tuningparList, n.core, cmds.result = cmds.result$points,
                       metric = dist.metric[dist.metric.index])
  end.time <- Sys.time()
  end.time - start.time

  tsn.logZ[i] = asmc.result2$logZ
  psi[i] = mean(asmc.result2$psi.output)

  # posterior inference
  tsn.index.asmc <- which.min(asmc.result2$SSR.output)
  tsn.asmc.res <- asmc.result2$xi.output[[tsn.index.asmc]]

  tsn.stress[i] <- round(stressFun(d.mat = dis.true,
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
  # tt.stress[i] <- round(stressFun(d.mat = dis.true,
  #                                 delta.mat = as.matrix(dist(tt.asmc.res, method = "euclidean"))), 4)

}

# # make plot
# data.frame(TN = tn.logZ, TSN = tsn.logZ) %>%
#   pivot_longer(everything()) %>%
#   ggplot(aes(x = name, y = value, fill = name)) +
#   geom_boxplot() +
#   labs(y = "log marginal likelihood", x = "model") +
#   scale_fill_discrete(name = "Model",
#                       labels = c("TN", "TSN")) +
#   theme_bw() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())






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
#library(smacof)
library(doMC)

## source ASMC models
# helper functions
source(file = "R/ASMC_helper_fun.R")
# ASMC function
source(file = "R/ASMC_fun_Rcpp.R")
source(file = "R/ASMC_incr_fun_Rcpp.R")
# truncated skewed Normal
source(file = "R/ASMC_truncatedSkewedN.R")
source(file = "R/ASMC_truncatedSkewedN_incr.R")
# truncated T
source(file = "R/ASMC_truncatedT.R")
source(file = "R/ASMC_truncatedT_incr.R")

Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_SN_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_SN_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_SN_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_SN_incr_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_SN_incr_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_SN_incr_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_T_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_T_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_T_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_T_incr_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_T_incr_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_T_incr_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/bisectionFun_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/dmvnrm_arma_fast_rcpp.cpp")

## define dissimilarity metric
dist.metric <-"cosine"

## set incremental size
n.t1 <- 50      # correspond to n0
incr.step <- 5  # correspond to n1

## create containers to store the results
n.run <- 30
time1 <- c(rep(NA, n.run))
stress1 <- c(rep(NA, n.run))
time2 <- c(rep(NA, n.run))
stress2 <- c(rep(NA, n.run))
time12 <- c(rep(NA, n.run))
stress12 <- c(rep(NA, n.run))

n1 <- 2
n2 <- round((n.t1+incr.step)/2)
n20 <- round((n.t1)/2)

for (i in 1:1){
  
  # permute data columnwise
  data <- readRDS(file = "data/NIPS_data_small.rds")
  data <- data[,sample(ncol(data))]
  
  # calculate Cosine distances for text
  dis <- 1-philentropy::distance(t(data), method = "cosine", mute.message = TRUE)
  
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
  
  reference.x.sd <- diag(rep(0.01, p))
  
  hyperpars.t1 <- list(a = sim.a.initial.t1, b = sim.b.initial.t1,
                       alpha = sim.alpha.initial.t1, beta = sim.beta.initial.t1, df = df.initial.t1,
                       c = c.initial.t1, d = d.initial.t1, constant_multiple = constant.multiple,
                       reference_x_sd = reference.x.sd)
  
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
                       c = c.initial.t2, d = d.initial.t2, constant_multiple = constant.multiple,
                       reference_x_sd = reference.x.sd)

  tuningparList1 <- list(K = 100, phi = 0.90, eps = 0.5)
  n.core <- detectCores() - 1

  ## set 1 (original part) only
  #model <- truncatedSkewedN(hyperparList = hyperpars.t1, p, reference.x.sd)
  model <- truncatedT(hyperparList = hyperpars.t1, p, reference.x.sd)
  
  n <- nrow(dis.t1)
  start.time <- Sys.time()
  t1.asmc.result <- ASMC_Rcpp(model = model,
                              dist.mat = dis.t1,
                              tuningparList = tuningparList1, 
                              n.core, cmds.result = t1.cmds.res,
                              metric = dist.metric,
                              upper_bound = 1, n.update = n1, 
                              n.update.x = n20)
  end.time <- Sys.time()
  time1[i] <- as.numeric(end.time - start.time, units = "secs")
  time1[i]
  
  # t1.asmc.result$accept_rate %>% mean()
  # plot(t1.asmc.result$accept_rate, type = "l",
  #      xlab = "ASMC iteration", ylab = "acceptce rate",
  #      main = paste0("n.update = ", n1, ", n.update.x = ", n2))
  # t1.asmc.result$rESS %>% mean()
  # plot(t1.asmc.result$rESS, type = "l",
  #      xlab = "ASMC iteration", ylab = "rESS",
  #      main = paste0("n.update = ", n1, ", n.update.x = ", n2))
  
  # posterior inference
  index.asmc.t1 <- which.min(t1.asmc.result$SSR.output)
  t1.asmc.res <- t1.asmc.result$xi.output[[index.asmc.t1]]

  stress1[i] <- stressFun(d.mat = dis.t1,
                          delta.mat = 1-philentropy::distance(t1.asmc.res,
                                                            method = dist.metric, mute.message = TRUE))
  stress1[i]
  #plot(t1.asmc.res[,1], t1.asmc.res[,2])
  
  # store results
  t1.asmc.res.all <- list(xi = t1.asmc.res, sigma2 = t1.asmc.result$sigma2.output[[index.asmc.t1]],
                          lambda = t1.asmc.result$lambda.output[[index.asmc.t1]],
                          psi = t1.asmc.result$psi.output[[index.asmc.t1]],
                          g = t1.asmc.result$g.output[[index.asmc.t1]])

  ## set 2 (original + incremental part) only
  
  n <- nrow(dis.t2)
  start.time <- Sys.time()
  t2.asmc.result <- ASMC_Rcpp(model = model,
                              dist.mat = dis.t2,
                              tuningparList = tuningparList1, 
                              n.core, cmds.result = t2.cmds.res,
                              metric = dist.metric,
                              upper_bound = 1, n.update = n1, 
                              n.update.x = n2)
  end.time <- Sys.time()
  time2[i] <- as.numeric(end.time - start.time, units = "secs")

  # posterior inference
  index.asmc.t2 <- which.min(t2.asmc.result$SSR.output)
  t2.asmc.res <- t2.asmc.result$xi.output[[index.asmc.t2]]

  stress2[i] <- stressFun(d.mat = dis.t2,
                          delta.mat = 1-philentropy::distance(t2.asmc.res,
                                                              method = dist.metric, mute.message = TRUE))
  ## set 1 and 2 combine
  #reference.x.sd <- cov(t1.asmc.res.all$xi)
  hyperpars.t12 <- list(a = sim.a.initial.t2, b = sim.b.initial.t2,
                        alpha = sim.alpha.initial.t2, beta = sim.beta.initial.t2, df = df.initial.t2,
                        c = c.initial.t2, d = d.initial.t2, constant_multiple = constant.multiple,
                        reference_x_sd = cov(t1.asmc.res.all$xi))
  
  #model <- truncatedSkewedN_incr(hyperparList = hyperpars.t12, p, reference.x.sd)
  model <- truncatedT_incr(hyperparList = hyperpars.t12, p,
                           reference.x.sd = cov(t1.asmc.res.all$xi))

  n.incr <- nrow(dis.t2) - nrow(dis.t1)
  tuningparList2 <- list(K = 100, # reduce the particle size to half
                         phi = 0.90, eps = 0.5)

  start.time <- Sys.time()
  t12.asmc.result <- ASMC_incr_Rcpp(model = model,
                                    dist.mat = dis.t2,
                                    tuningparList =  tuningparList2, n.core, 
                                    prev.result = t1.asmc.res.all$xi,
                                    metric = dist.metric,
                                    upper_bound = 1, n.update = n1, 
                                    n.update.x = n2)
  end.time <- Sys.time()
  time12[i] <- as.numeric(end.time - start.time, units = "secs")

  # posterior inference
  index.asmc.t12 <- which.min(t12.asmc.result$SSR.output)
  t12.asmc.res <- t12.asmc.result$xi.output[[index.asmc.t12]]

  stress12[i] <- stressFun(d.mat = dis.t2,
                           delta.mat = 1-philentropy::distance(t12.asmc.res,
                                                               method = dist.metric, 
                                                               mute.message = TRUE))
  
}

# computation time
a <- 1
b <- n.run
summary(time2[a:b])  # without adaptive inference
summary(time12[a:b]) # with adaptive inference


sd(time2[a:b])
sd(time12[a:b])


# STRESS value
summary(stress2[a:b])  # without adaptive inference
summary(stress12[a:b]) # with adaptive inference




# This script contains codes for simulation study
# This example is for section 4.1 Experiment 3: Data with known trye underlying dimension

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
library(Rcpp)
library(doMC)
library(bayMDS)

## create containers to store the results
n.run <- 20
cmds.stress <- matrix(data = NA, nrow = 7, ncol = n.run)
tn.logZ = tn.stress = tsn.logZ = psi = tsn.stress = 
  tt.logZ = tt.stress <- matrix(data = NA, nrow = 7, ncol = n.run)
tn.time = tsn.time = tt.time <- matrix(data = NA, nrow = 7, ncol = n.run)
bmds.stress1 = bmds.stress2 = bmds.time <- matrix(data = NA, nrow = 7, ncol = n.run)

## source ASMC models
# helper functions
source(file = "R/ASMC_helper_fun.R")
# ASMC function with rcpp implementation
source(file = "R/ASMC_fun.R")
source(file = "R/ASMC_fun_Rcpp.R")
# truncated Normal
source(file = "R/ASMC_truncatedN.R")
# truncated T
source(file = "R/ASMC_truncatedT.R")
# truncated skewed Normal
source(file = "R/ASMC_truncatedSkewedN.R")

# helper rcpp functions
Rcpp::sourceCpp(file = "helper_Rcpp/bisectionFun_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/dmvnrm_arma_fast_rcpp.cpp")

# model truncated Normal rcpp functions
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_cpp.cpp")

# model truncated skewed Normal rcpp functions
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_SN_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_SN_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_SN_cpp.cpp")

# model truncated T rcpp functions
Rcpp::sourceCpp(file = "helper_Rcpp/initialFun_T_cpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/likelihoodFun_T_rcpp.cpp")
Rcpp::sourceCpp(file = "helper_Rcpp/proposalFun_T_cpp.cpp")

n1 = 4 # number of parameters to propose in each ASMC iteration
n2 = 100 # number of observations (X) to propose in each ASMC iteration

## run simulations 20 times
for (i in 1:n.run){
  
  ## choose dissimilarity metric
  dist.metric <- "euclidean"
  
  ## simulate 100 random samples from a 5-dim multivariate Gaussian distribution
  ## with mean 0 and variance I, the identity matrix
  X <- rmvnorm(n = 100, mean = rep(0, 5), sigma=diag(5))
  
  ## get true pairwise Euclidean dissimilarities
  true.dis <- as.matrix(dist(X, method = dist.metric))
  
  ## Given true pairwise Euclidean dissimilarities, 
  ## generate the observed dissimilarities from some normal distribution
  f1 <- function(x){
    res <- 0
    while (res <= 0){
      res <- rnorm(1, mean = x, sd = 1)
    }
    return(res)
  }
  obs.dis <- apply(true.dis, 1:2, f1)
  diag(obs.dis) <- 0
  
  # observed dissimilarity
  dis <- obs.dis
  
  # true dissimilarity
  dis.true <- true.dis
  
  for (p in 2:8){
   
    ## 0. MDS with annealed SMC
    # use cmds to get initial values
    cmds.result <- cmdscale(d = dis, k = p,
                            eig = TRUE, add = FALSE, x.ret = FALSE)
    class(cmds.result) <- append(class(cmds.result), "CMDS")
    
    cmds.stress[(p-1),i] <- round(stressFun(d.mat = dis.true,
                                      delta.mat = as.matrix(dist(cmds.result$points, method = "euclidean"))), 4)
    
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

    reference.x.sd <- diag(rep(0.01, p))
    
    hyperparList <- list(a = sim.a.initial, b = sim.b.initial,
                         alpha = sim.alpha.initial, beta = sim.beta.initial,
                         df = df.initial,
                         c = c.initial, d = d.initial, constant_multiple = constant.multiple,
                         reference_x_sd = reference.x.sd)
    
    tuningparList <- list(K = 200, # number of particles
                          phi = 0.80, eps = 0.5)
    n.core <- detectCores() - 1
    
    ## 1. annealed SMC with truncated Normal
    model <- truncatedN(hyperparList, p, reference.x.sd)
    start.time <- Sys.time()
    asmc.result1 <- ASMC_Rcpp(model = model,
                              dist.mat = dis,
                              tuningparList, n.core, cmds.result = cmds.result$points,
                              metric = dist.metric,
                              upper_bound = 1e5, n.update = n1,
                              n.update.x = n2)
    end.time <- Sys.time()
    tn.time[(p-1),i] <- difftime(end.time, start.time, units = "secs")
    
    tn.logZ[(p-1),i] = asmc.result1$logZ
    
    # posterior inference
    tn.index.asmc <- which.min(asmc.result1$SSR.output)
    tn.asmc.res <- asmc.result1$xi.output[[tn.index.asmc]]
    
    (tn.stress[(p-1),i] <- round(stressFun(d.mat = dis.true,
                                    delta.mat = as.matrix(dist(tn.asmc.res, method = "euclidean"))), 4))
    
    tn.iteration <- asmc.result1$iteration
    
    #set.seed(549)
    ## 2. anneal SMC with truncated skewed Normal
    model <- truncatedSkewedN(hyperparList, p, reference.x.sd)
    
    # run annealed SMC
    start.time <- Sys.time()
    asmc.result2 <- ASMC_Rcpp(model = model,
                              dist.mat = dis,
                              tuningparList, n.core, cmds.result = cmds.result$points,
                              metric = dist.metric,
                              upper_bound = 1e5, n.update = n1, 
                              n.update.x = n2)
    end.time <- Sys.time()
    tsn.time[(p-1),i] <- difftime(end.time, start.time, units = "secs")
    
    tsn.logZ[(p-1),i] = asmc.result2$logZ
    psi[(p-1),i] = mean(asmc.result2$psi.output)
    
    # posterior inference
    tsn.index.asmc <- which.min(asmc.result2$SSR.output)
    tsn.asmc.res <- asmc.result2$xi.output[[tsn.index.asmc]]
    
    tsn.stress[(p-1),i] <- round(stressFun(d.mat = dis.true,
                                     delta.mat = as.matrix(dist(tsn.asmc.res, method = "euclidean"))), 4)
    
    tsn.iteration <- asmc.result2$iteration
    
    ## 3. anneal SMC with truncated T
    model <- truncatedT(hyperparList, p, reference.x.sd)
    
    start.time <- Sys.time()
    asmc.result3 <- ASMC_Rcpp(model = model,
                              dist.mat = dis,
                              tuningparList, n.core, cmds.result = cmds.result$points,
                              metric = dist.metric,
                              upper_bound = 1e5, n.update = n1, 
                              n.update.x = n2)
    end.time <- Sys.time()
    tt.time[(p-1),i] <- end.time - start.time
    
    tt.logZ[(p-1),i] = asmc.result3$logZ
    
    # posterior inference
    tt.index.asmc <- which.min(asmc.result3$SSR.output)
    tt.asmc.res <- asmc.result3$xi.output[[tt.index.asmc]]
    
    tt.stress[(p-1),i] <- round(stressFun(d.mat = dis.true,
                                    delta.mat = as.matrix(dist(tt.asmc.res, method = "euclidean"))), 4)
    
  }
  
  ## 4. bayMDS R package
  start.time <- Sys.time()
  out <- bmds(dis,min_p=2,max_p=8) 
  end.time <- Sys.time()
  bmds.time[(p-1),i] <- difftime(end.time, start.time, units = "secs")
  bmds.stress1[,i] <- unlist(out$stress.L)
  for (p in 2:8){
    bmds.stress2[(p-1),i] <- round(stressFun(d.mat = dis.true,
                                             delta.mat = as.matrix(dist(out$x_bmds[[p]], 
                                                                        method = "euclidean"))), 4)
    
  }
  rm(out)
  gc()
  
}





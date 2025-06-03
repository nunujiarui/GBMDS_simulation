# This script is for Experiment 3: Data with outliers

# devtools::install_github("https://github.com/SFU-Stat-ML/GBMDS")

library(parallel)
library(doMC)
library(fGarch)
library(GBMDS)

## choose dissimilarity metric
dist.metric <- c("euclidean")
# change this index to try different distance metrics
dist.metric.index <- 1

## specify percentage of outliers to add
outlier.pct <- 0.05

## create containers to store the results
tsn.logZ = tsn.stress = tt.logZ = tt.stress <- numeric()
tsn.time = tt.time <- numeric()
bmds.stress1 = bmds.stress2 = bmds.stress3 = bmds.stress4 = bmds.time1 = bmds.time2 <- numeric()


## run simulations 20 times
for (i in 1:20){
  
  ## simulate 100 random samples from a 10-dim multivariate Gaussian distribution
  ## with mean 0 and variance I, the identity matrix
  data <- rmvnorm(n = 100, mean = rep(0, 10), sigma=diag(10))
  true.dis <- as.matrix(dist(data, method = dist.metric[dist.metric.index], p = 2))
  
  ## Given true pairwise Euclidean dissimilarities, 
  ## generate the observed dissimilarities from some normal distribution
  f1 <- function(x){
    res <- 0
    while (res <= 0){
      res <- rnorm(1, mean = x, sd = 0.5)
    }
    return(res)
  }
  obs.dis <- apply(true.dis, 1:2, f1)
  diag(obs.dis) <- 0
  
  n.dis <- nrow(obs.dis) * (nrow(obs.dis) -1) / 2
  
  n.outlier <- round(n.dis * outlier.pct, 0)
  
  outlier.index <- sample(n.dis, n.outlier)
  
  dis <- matrix(data = 0, nrow = nrow(obs.dis), ncol = ncol(obs.dis))
  
  dis[upper.tri(dis)] <- obs.dis[upper.tri(obs.dis)]
  dis[upper.tri(dis)][outlier.index] <- obs.dis[upper.tri(obs.dis)][outlier.index] * 2
  dis <- dis + t(dis)
  
  dis.df <- data.frame(dis = dis[upper.tri(dis)])
  
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
  reference.x.sd <- diag(rep(0.01, p))
  
  hyperparList <- list(a = sim.a.initial, b = sim.b.initial,
                       alpha = sim.alpha.initial, beta = sim.beta.initial,
                       df = df.initial,
                       c = c.initial, d = d.initial, constant_multiple = constant.multiple,
                       reference_x_sd = reference.x.sd)
  
  tuningparList <- list(K = 200, phi = 0.80, eps = 0.5)
  n.core <- detectCores()-1
  
  ## 2. anneal SMC with truncated skewed Normal
  model <- truncatedSkewedN(hyperparList, p, reference.x.sd)
  
  # run annealed SMC
  start.time <- Sys.time()
  asmc.result2 <- ASMC_Rcpp(model = model,
                            dist.mat = dis,
                            tuningparList, n.core, cmds.result = cmds.result$points,
                            metric = dist.metric[dist.metric.index],
                            upper_bound = 1e5, n.update = 4, 
                            n.update.x = nrow(dis))
    
  end.time <- Sys.time()
  tsn.time[i] <- difftime(end.time, start.time, units = "secs")
  
  tsn.logZ[i] = asmc.result2$logZ
  
  # posterior inference
  tsn.index.asmc <- which.min(asmc.result2$SSR.output)
  tsn.asmc.res <- asmc.result2$xi.output[[tsn.index.asmc]]
  
  tsn.stress[i] <- round(stressFun(d.mat = dis,
                                   delta.mat = as.matrix(dist(tsn.asmc.res, method = "euclidean"))), 4)
  
  tsn.iteration <- asmc.result2$iteration
  
  ## 3. anneal SMC with truncated T
  model <- truncatedT(hyperparList, p, reference.x.sd)
  
  start.time <- Sys.time()
  asmc.result3 <- ASMC_Rcpp(model = model,
                            dist.mat = dis,
                            tuningparList, n.core, cmds.result = cmds.result$points,
                            metric = dist.metric[dist.metric.index],
                            upper_bound = 1e5, n.update = 4, 
                            n.update.x = nrow(dis))
  end.time <- Sys.time()
  tt.time[i] <- difftime(end.time, start.time, units = "secs")
  
  tt.logZ[i] = asmc.result3$logZ
  
  # posterior inference
  tt.index.asmc <- which.min(asmc.result3$SSR.output)
  tt.asmc.res <- asmc.result3$xi.output[[tt.index.asmc]]
  
  tt.stress[i] <- round(stressFun(d.mat = dis,
                                  delta.mat = as.matrix(dist(tt.asmc.res, method = "euclidean"))), 4)
  
  tt.iteration <- asmc.result3$iteration
  
  ## 4. bayMDS R package
  start.time <- Sys.time()
  out1 <- bmds(dis, min_p=2, max_p=2, niter = tn.iteration * 200) 
  end.time <- Sys.time()
  bmds.time1[i] <- difftime(end.time, start.time, units = "secs")
  bmds.stress1[i] <- out1$stress.L[[2]]
  bmds.stress2[i] <- round(stressFun(d.mat = dis,
                                     delta.mat = as.matrix(dist(out1$x_bmds[[2]], method = "euclidean"))), 4)
  rm(out1)
  gc()
  
  start.time <- Sys.time()
  out2 <- bmds(dis, min_p=2, max_p=2, niter = tt.iteration * 200) 
  end.time <- Sys.time()
  bmds.time2[i] <- difftime(end.time, start.time, units = "secs")
  bmds.stress3[i] <- out2$stress.L[[2]]
  bmds.stress4[i] <- round(stressFun(d.mat = dis,
                                     delta.mat = as.matrix(dist(out2$x_bmds[[2]], method = "euclidean"))), 4)
  rm(out2)
  gc()
  
}

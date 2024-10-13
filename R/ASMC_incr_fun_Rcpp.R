# This script contains the main function to perform annealed SMC algorithm with adaptive inference.
#'
#' @param model likelihood model, including " " ...
#' @param dist.mat      distance matrix
#' @param tuningparList  SMC tuning parameters
#' @param n.core         the number of cores
#' @param prev.result    results from the previous ASMC
#' @param metric         distance metric used in the model
#' @return results of weighted particles, marginal likelihood estimates
#' @examples
#' print(" ")


ASMC_incr_Rcpp <- function(model, dist.mat, tuningparList, n.core, 
                           prev.result, metric, upper_bound){

  n <- nrow(dist.mat)
  K <- tuningparList$K
  
  hyperparList <- model$hyperparList

  ### setting

  x.mean <- list(rep(0, p))
  x.mean.list <- rep(x.mean, n)

  ### initialization
  # initialize SMC iteration index
  r <- 1

  # initialize annealing parameter
  tau <- numeric()
  tau[1] <- 0
  tauDiff <- numeric()
  tauDiff[1] <- 0

  # initialize marginal log likelihood estimate
  logZ <- 0

  # initialize the relative effective sample size
  rESS <- numeric()
  rESS[1] <- 0

  # initialize particles with independent samples
  theta <- list()
  xi <- list()
  sigma2 <- numeric()
  sigma2.list <- list()
  lambda <- list()
  psi <- numeric()

  if (any(class(model) %in% c("truncatedT_incr"))){
    g <- list()
  }
  
  initialK <- function(k){
    #res <- initialFun(model, cmds.result, dist.mat, metric)
    if (any(class(model) %in% c("truncatedN_incr"))){
      res <- initialFun_incr_cpp(prev.result, dist.mat, metric, hyperparList)
    } else if (any(class(model) %in% c("truncatedT_incr"))){
      res <- initialFun_T_incr_cpp(prev.result, dist.mat, metric, hyperparList)
    } else if (any(class(model) %in% c("truncatedSkewedN_incr"))){
      res <- initialFun_SN_incr_cpp(prev.result, dist.mat, metric, hyperparList)
    }
    return(res)
  }
  
  if (n.core > 1){
    theta <- foreach(i = 1:K, .combine = 'c') %dopar% {
      temp <- initialK(i)
      list(temp)
    }
    xi <- lapply(1:K, function(k){theta[[k]]$x})
    sigma2 <- unlist(lapply(1:K, function(k){theta[[k]]$sigma2}))
    psi <- unlist(lapply(1:K, function(k){theta[[k]]$psi}))
    lambda <- lapply(1:K, function(k){theta[[k]]$lambda})
    if (any(class(model) %in% c("truncatedT_incr"))){
      g <- lapply(1:K, function(k){theta[[k]]$g})
    }
  } else{
    if (any(class(model) %in% c("truncatedN_incr"))){
      theta <- lapply(1:K, function(i){initialFun_cpp(cmds.result, dist.mat, metric, hyperparList)})
    } else if (any(class(model) %in% c("truncatedT_incr"))){
      theta <- lapply(1:K, function(i){initialFun_T_cpp(cmds.result, dist.mat, metric, hyperparList)})
    } else if (any(class(model) %in% c("truncatedSkewedN_incr"))){
      theta <- lapply(1:K, function(i){initialFun_SN_cpp(cmds.result, dist.mat, metric, hyperparList)})
    }
    #theta <- lapply(1:K, function(i){initialFun_cpp(cmds.result, dist.mat, metric, hyperparList)})
    xi <- lapply(1:K, function(k){theta[[k]]$x})
    sigma2 <- unlist(lapply(1:K, function(k){theta[[k]]$sigma2}))
    psi <- unlist(lapply(1:K, function(k){theta[[k]]$psi}))
    lambda <- lapply(1:K, function(k){theta[[k]]$lambda})
    if (any(class(model) %in% c("truncatedT_incr"))){
      g <- lapply(1:K, function(k){theta[[k]]$g})
    }
  }

  sigma2.list[[r]] <- sigma2

  n.incr <- nrow(dist.mat) - nrow(prev.result)
  reference.mean <- rep(list(rbind(prev.result,
                                   matrix(data = 0, nrow = n.incr, ncol = model$p))), K)
  #reference.mean <- rep(list(prev.result$xi), K)


  # set initialize weights to unity
  W <- rep(1/K, K)
  logW <- log(W)

  ### Annealed SMC
  while (tau[r] < 1) {

    cat("iteration:",r,"\n")
    r <- r+1

    # evaluate the log-likelihood
    logL <- rep(NA, K)
    SSR <- rep(NA, K)

    logLFun <- function(k){
      n.incr <- nrow(dist.mat) - nrow(prev.result)
      if (any(class(model) %in% c("truncatedT_incr"))){
        previousVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                            lambda = lambda[[k]], SSR = SSR[k], g = g[[k]])
      } else{
        previousVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                            lambda = lambda[[k]], SSR = SSR[k])
      }
      
      if (any(class(model) %in% c("truncatedN_incr"))){
        lik.result <- likelihoodFun_incr_cpp(dist.mat, upper_bound, previousVal, metric, hyperparList, n.incr)
      } else if (any(class(model) %in% c("truncatedT_incr"))){
        lik.result <- likelihoodFun_T_incr_cpp(dist.mat, upper_bound, previousVal, metric, hyperparList, n.incr)
      } else if (any(class(model) %in% c("truncatedSkewedN_incr"))){
        lik.result <- likelihoodFun_SN_incr_cpp(dist.mat, upper_bound, previousVal, metric, hyperparList, n.incr)
      }
      
      #lik.result <- likelihoodFun(model,  dist.mat, proposal.result = previousVal, metric, n.incr)
      xi.list <- lapply(seq_len(nrow(xi[[k]])), function(i) xi[[k]][i,])
      xi.prev.list <- xi.list[1:(n-n.incr)]
      #lambda.list <- rep(list(lambda[[k]]), n)
      #prior.d.x <- mapply(mvtnorm::dmvnorm, xi.list, x.mean.list, lambda.list, log = TRUE)
      reference.mean.k <- lapply(seq_len(nrow(reference.mean[[k]])), function(i) reference.mean[[k]][i,])[1:(n-n.incr)]
      ref.d.x <- mapply(mvtnorm::dmvnorm, xi.prev.list, reference.mean.k, rep(list(reference.x.sd), n-n.incr), log = TRUE)
      logPriorRef <- sum(ref.d.x)
      logL <- lik.result$loglikelihood + sum(lik.result$logprior) - logPriorRef
      SSR <- lik.result$SSR
      
      result <- list(logL = logL, SSR = SSR)
      
      return(result)
    }
    
    if (n.core > 1){
      result <- foreach(i = 1:K, .combine = 'cbind') %dopar% {
        logLFun(i)
      }
      logL <- unlist(result[1, ])
      SSR <- unlist(result[2, ])
    } else{
      result <- sapply(1:K, logLFun)
      logL <- unlist(result[1, ])
      SSR <- unlist(result[2, ])
    }

    ## determine next annealing parameter
    tauDiff[r] <- bisectionFun_cpp(low = 0, high = 1, W, logL, phi = tuningparList$phi)
    tau[r] <- tau[r-1] + tauDiff[r]
    cat("annealing parameter:", tau[r],"\n")

    # if tau is set greater than 1, fix by setting to 1
    if(tau[r] > 1){
      tau[r] <- 1
      tauDiff[r] <- 1-tau[r-1]
    }

    ## compute pre-sampling unnormalized weights
    logw <- logW + tauDiff[r]*logL
    # normalize the weights
    logmax <- max(logw)
    W <- exp(logw-logmax)/sum(exp(logw-logmax))
    logW <- log(W)

    ## sample particles using invariant Metropolis-Hastings kernel
    MCMCmove <- function(k){

      n.incr <- nrow(dist.mat) - nrow(prev.result)
      if (any(class(model) %in% c("truncatedT_incr"))){
        currentVal <- list(x = xi[[k]], sigma2 = sigma2[k],
                           lambda = lambda[[k]], SSR = SSR[k], g = g[[k]])
      } else if (any(class(model) %in% c("truncatedSkewedN_incr"))){
        currentVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                           lambda = lambda[[k]], SSR = SSR[k])
      } else{
        currentVal <- list(x = xi[[k]], sigma2 = sigma2[k],
                           lambda = lambda[[k]], SSR = SSR[k])
      }

      if (r == 2){
        prevX <- reference.mean[[k]]
      } else{
        prevX <- xi.prev[[k]]
      }

      if (any(class(model) %in% c("truncatedN_incr"))){
        prop.result <- proposalFun_incr_cpp(dist.mat, currentVal, prevX, 
                                       tau[r], metric, hyperparList, n.incr, upper_bound)
      } else if (any(class(model) %in% c("truncatedT_incr"))){
        prop.result <- proposalFun_T_incr_cpp(dist.mat, currentVal, prevX, 
                                              tau[r], metric, hyperparList, n.incr, upper_bound)
      } else if (any(class(model) %in% c("truncatedSkewedN_incr"))){
        prop.result <- proposalFun_SN_incr_cpp(dist.mat, currentVal, prevX, 
                                          tau[r], metric, hyperparList, n.incr, upper_bound)
      }
      
      xi <- prop.result$x
      sigma2 <- prop.result$sigma2
      lambda <- prop.result$lambda
      if (any(class(model) %in% c("truncatedT_incr"))){
        g <- prop.result$g
      }
      psi <- prop.result$psi

      # save the previous particles' coordinates x
      xi.prev <- currentVal$x

      if (any(class(model) %in% c("truncatedT_incr"))){
        result <- list(xi = xi, sigma2 = sigma2, lambda = lambda, xi.prev = xi.prev, g = g, psi = psi)
      } else{
        result <- list(xi = xi, sigma2 = sigma2, lambda = lambda, xi.prev = xi.prev, psi = psi)
      }

      return(result)

    }

    if (n.core > 1){
      result.MCMC <- foreach(i = 1:K, .combine = 'c') %dopar% {
        temp <- MCMCmove(i)
        list(temp)
      }
      xi <- lapply(1:K, function(k){result.MCMC[[k]]$xi})
      sigma2 <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$sigma2}))
      lambda <- lapply(1:K, function(k){result.MCMC[[k]]$lambda})
      xi.prev <- lapply(1:K, function(k){result.MCMC[[k]]$xi.prev})
      psi <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$psi}))
      if (any(class(model) %in% c("truncatedT_incr"))){
        g <- lapply(1:K, function(k){result.MCMC[[k]]$g})
      }
    } else{
      result.MCMC <- lapply(1:K, MCMCmove)
      xi <- lapply(1:K, function(k){result.MCMC[[k]]$xi})
      sigma2 <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$sigma2}))
      lambda <- lapply(1:K, function(k){result.MCMC[[k]]$lambda})
      xi.prev <- lapply(1:K, function(k){result.MCMC[[k]]$xi.prev})
      psi <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$psi}))
      if (any(class(model) %in% c("truncatedT_incr"))){
        g <- lapply(1:K, function(k){result.MCMC[[k]]$g})
      }
    }

    sigma2.list[[r]] <- sigma2

    ## update the log marginal likelihood estimate
    logZ <- logZ + log(sum(exp(logw - logmax))) + logmax

    ## compute the rESS
    rESS[r] <- rESSFun(logW)

    # check the value for tau_r
    if (tau[r] == 1){
      # can end the algorithm
      # return the current particle population
      if (any(class(model) %in% c("truncatedT"))){
        output.list <- list(xi.output = xi, sigma2.output = sigma2, psi.output = psi,
                            lambda.output = lambda, g = g,
                            weight.output = W,
                            SSR.output = SSR, logZ = logZ)
      } else{
        output.list <- list(xi.output = xi, sigma2.output = sigma2, psi.output = psi,
                            lambda.output = lambda,
                            weight.output = W,
                            SSR.output = SSR, logZ = logZ)
      }
    } else if (rESS[r] < tuningparList$eps){
      # particle degeneracy is too severe, resampling is needed
      index <- multinomialResampleFun(W)
      xi <- xi[index]
      sigma2 <- sigma2[index]
      lambda <- lambda[index]
      sigma2.list[[r]] <- sigma2

      # reset particle weights
      W <- rep(1/K, K)
      logW <- log(W)

    }

  }

  ## Set the name for the class
  class(output.list) <- append(class(output.list),"BMDSParticles_incr")

  return(output.list)

}


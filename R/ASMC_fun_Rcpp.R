#' This script contains the main function to perform annealed SMC algorithm.
#'
#' @param model  likelihood model, including " " ...
#' @param hyperparList  hyper parameters
#' @param dist.mat      distance matrix
#' @param tuningparList  SMC tuning parameters
#' @param n.core         the number of cores
#' @param cmds.result    results from the classical MDS
#' @param metric         distance metric used in the model
#' @param delta.mat  delta matrix
#' @return results of weighted particles, marginal likelihood estimates
#' @examples
#' print(" ")
#' 

ASMC_Rcpp <- function(model, dist.mat, tuningparList, n.core, cmds.result, metric){

  # register the number of cores to use for parallel execution
  registerDoMC(n.core)
  
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

  if (any(class(model) %in% c("truncatedT"))){
    g <- list()
  }
  
  initialK <- function(k){
    #res <- initialFun(model, cmds.result, dist.mat, metric)
    if (any(class(model) %in% c("truncatedN"))){
      res <- initialFun_cpp(cmds.result, dist.mat, metric, hyperparList)
    } else if (any(class(model) %in% c("truncatedT"))){
      res <- initialFun_T_cpp(cmds.result, dist.mat, metric, hyperparList)
    } else if (any(class(model) %in% c("truncatedSkewedN"))){
      res <- initialFun_SN_cpp(cmds.result, dist.mat, metric, hyperparList)
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
    if (any(class(model) %in% c("truncatedT"))){
      g <- lapply(1:K, function(k){theta[[k]]$g})
    }
  } else{
    if (any(class(model) %in% c("truncatedN"))){
      theta <- lapply(1:K, function(i){initialFun_cpp(cmds.result, dist.mat, metric, hyperparList)})
    } else if (any(class(model) %in% c("truncatedT"))){
      theta <- lapply(1:K, function(i){initialFun_T_cpp(cmds.result, dist.mat, metric, hyperparList)})
    } else if (any(class(model) %in% c("truncatedSkewedN"))){
      theta <- lapply(1:K, function(i){initialFun_SN_cpp(cmds.result, dist.mat, metric, hyperparList)})
    }
    #theta <- lapply(1:K, function(i){initialFun_cpp(cmds.result, dist.mat, metric, hyperparList)})
    xi <- lapply(1:K, function(k){theta[[k]]$x})
    sigma2 <- unlist(lapply(1:K, function(k){theta[[k]]$sigma2}))
    psi <- unlist(lapply(1:K, function(k){theta[[k]]$psi}))
    lambda <- lapply(1:K, function(k){theta[[k]]$lambda})
    if (any(class(model) %in% c("truncatedT"))){
      g <- lapply(1:K, function(k){theta[[k]]$g})
    }
  }

  sigma2.list[[r]] <- sigma2

  reference.mean <- rep(list(cmds.result), K)


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
      if (any(class(model) %in% c("truncatedT"))){
        previousVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                           lambda = lambda[[k]], SSR = SSR[k], g = g[[k]])
      } else{
        previousVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                           lambda = lambda[[k]], SSR = SSR[k])
      }
      
      if (any(class(model) %in% c("truncatedN"))){
        lik.result <- likelihoodFun_cpp(dist.mat, previousVal, metric, hyperparList)
      } else if (any(class(model) %in% c("truncatedT"))){
        lik.result <- likelihoodFun_T_cpp(dist.mat, previousVal, metric, hyperparList)
      } else if (any(class(model) %in% c("truncatedSkewedN"))){
        lik.result <- likelihoodFun_SN_cpp(dist.mat, previousVal, metric, hyperparList)
      }
      
      #lik.result <- likelihoodFun(model,  dist.mat, proposal.result = previousVal, metric)
      #lik.result <- likelihoodFun_cpp(dist.mat, previousVal, metric, hyperparList)
      xi.list <- lapply(seq_len(nrow(xi[[k]])), function(i) matrix(xi[[k]][i,], nrow = 1))
      lambda.list <- rep(list(lambda[[k]]), n)
      # prior.d.x <- mapply(mvtnorm::dmvnorm, xi.list, x.mean.list, lambda.list, log = TRUE)
      prior.d.x <- mapply(dmvnrm_arma_fast, xi.list, x.mean.list, lambda.list, log = TRUE)
      reference.mean.k <- lapply(seq_len(nrow(reference.mean[[k]])), function(i) reference.mean[[k]][i,])
      # ref.d.x <- mapply(mvtnorm::dmvnorm, xi.list, reference.mean.k, rep(list(reference.x.sd), n), log = TRUE)
      ref.d.x <- mapply(dmvnrm_arma_fast, xi.list, reference.mean.k, rep(list(reference.x.sd), n), log = TRUE)
      logPriorRef <- sum(ref.d.x)
      logL <- lik.result$logposterior - logPriorRef
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

      if (any(class(model) %in% c("truncatedT"))){
        currentVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                           lambda = lambda[[k]], SSR = SSR[k], g = g[[k]])
      } else{
        currentVal <- list(x = xi[[k]], sigma2 = sigma2[k], psi = psi[k],
                           lambda = lambda[[k]], SSR = SSR[k])
      }

      if (r == 2){
        prevX <- reference.mean[[k]]
      } else{
        prevX <- xi.prev[[k]]
      }
      
      # prop.result <- list()
      # attempt <- 0
      # while(length(prop.result) == 0 & attempt < 5){
      #   try(
      #     if (any(class(model) %in% c("truncatedN"))){
      #       prop.result <- proposalFun_cpp(dist.mat, currentVal, prevX, 
      #                                      tau[r], metric, hyperparList)
      #     } else if (any(class(model) %in% c("truncatedT"))){
      #       prop.result <- proposalFun_T_cpp(dist.mat, currentVal, prevX, 
      #                                      tau[r], metric, hyperparList)
      #     } else if (any(class(model) %in% c("truncatedSkewedN"))){
      #       prop.result <- proposalFun_SN_cpp(dist.mat, currentVal, prevX, 
      #                                      tau[r], metric, hyperparList)
      #     }
      #     
      #     # prop.result <- proposalFun_cpp(dist.mat, currentVal, prevX, 
      #     #                                tau[r], metric, hyperparList)
      #   )
      #   attempt <- attempt + 1
      # }
      # print(attempt)
      
      if (any(class(model) %in% c("truncatedN"))){
        prop.result <- proposalFun_cpp(dist.mat, currentVal, prevX, 
                                       tau[r], metric, hyperparList)
      } else if (any(class(model) %in% c("truncatedT"))){
        prop.result <- proposalFun_T_cpp(dist.mat, currentVal, prevX, 
                                         tau[r], metric, hyperparList)
      } else if (any(class(model) %in% c("truncatedSkewedN"))){
        prop.result <- proposalFun_SN_cpp(dist.mat, currentVal, prevX, 
                                          tau[r], metric, hyperparList)
      }

      # prop.result <- proposalFun_cpp(dist.mat, currentVal, prevX,
      #                                annealingPar = tau[r], metric, hyperparList)
      # 
      xi <- prop.result$x
      sigma2 <- prop.result$sigma2
      lambda <- prop.result$lambda
      if (any(class(model) %in% c("truncatedT"))){
        g <- prop.result$g
      }
      psi <- prop.result$psi

      # save the previous particles' coordinates x
      xi.prev <- currentVal$x

      if (any(class(model) %in% c("truncatedT"))){
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
      if (any(class(model) %in% c("truncatedT"))){
        g <- lapply(1:K, function(k){result.MCMC[[k]]$g})
      }
    } else{
      result.MCMC <- lapply(1:K, MCMCmove)
      xi <- lapply(1:K, function(k){result.MCMC[[k]]$xi})
      sigma2 <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$sigma2}))
      lambda <- lapply(1:K, function(k){result.MCMC[[k]]$lambda})
      xi.prev <- lapply(1:K, function(k){result.MCMC[[k]]$xi.prev})
      psi <- unlist(lapply(1:K, function(k){result.MCMC[[k]]$psi}))
      if (any(class(model) %in% c("truncatedT"))){
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
                            lambda.output = lambda, g.output = g,
                            weight.output = W,
                            SSR.output = SSR, logZ = logZ, iteration = r)
      } else{
        output.list <- list(xi.output = xi, sigma2.output = sigma2, psi.output = psi,
                            lambda.output = lambda,
                            weight.output = W,
                            SSR.output = SSR, logZ = logZ, iteration = r)
      }
    } else if (rESS[r] < tuningparList$eps){
      # particle degeneracy is too severe, resampling is needed
      #print("resamped!")
      index <- multinomialResampleFun(W)
      xi <- xi[index]
      sigma2 <- sigma2[index]
      lambda <- lambda[index]
      psi <- psi[index]
      sigma2.list[[r]] <- sigma2

      # reset particle weights
      W <- rep(1/K, K)
      logW <- log(W)
    }

  }

  ## Set the name for the class
  class(output.list) <- append(class(output.list),"BMDSParticles")
  return(output.list)
}





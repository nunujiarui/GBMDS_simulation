
truncatedN <- function(hyperparList, p = 2, reference.x.sd){
  bmdsmodel <- list()
  bmdsmodel$hyperparList <- hyperparList
  bmdsmodel$p <- p
  bmdsmodel$reference.x.sd <- reference.x.sd
  ## Set the name for the class
  class(bmdsmodel) <- append(class(bmdsmodel), c("truncatedN", "BMDSModel"))
  return(bmdsmodel)
}


prior <- function(model) UseMethod("prior", model)

prior.truncatedN <- function(model){
  hyperparList <- model$hyperparList
  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  # get multiplier

  ### Initialize particles
  # sigma2
  sigma2.initial <- MCMCpack::rinvgamma(1, shape = a, scale = b)
  # lambda
  lambda.initial <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p)
    diag(lambda.initial)[j] <- MCMCpack::rinvgamma(1, shape = alpha, scale = beta[j])


  # x_i
  x.initial <- matrix(data = 0, nrow = n, ncol = p)
  for (i in 1:n)
    x.initial[i, ] <- mvtnorm::rmvnorm(1, mean = rep(0, p), sigma = lambda.initial)

  output <- list(sigma2 = sigma2.initial,
                 lambda = lambda.initial,
                 x = x.initial)
  return(output)
}


initialFun <- function(model, ...) UseMethod("initialFun", model)

initialFun.truncatedN <- function(model,cmds.result, dist.mat, metric){

  n.obj <- nrow(dist.mat)
  m <- n.obj*(n.obj - 1)/2

  ### Initialize particles
  ## x_i (from reference distribution)
  x.cmds <- cmds.result
  x.initial <- t(apply(x.cmds, 1, reference$x))

  # calculate delta matrix and d matrix
  d.mat <- as.matrix(dist.mat)
  if (dist.metric[dist.metric.index] == "cosine"){
    delta.mat <- 1-philentropy::distance(x.initial, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  } else{
    delta.mat <- philentropy::distance(x.initial, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  }

  # calculate SSR initial
  SSR.initial <- SSRFun(d.mat, delta.mat)

  hyperparList <- model$hyperparList
  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  ### Initialize particles
  ## sigma2 (from prior distribution)
  sigma2.initial <- MCMCpack::rinvgamma(1, shape = a, scale = b)
  ## lambda (from prior distribution)
  lambda.initial <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p){
    diag(lambda.initial)[j] <- MCMCpack::rinvgamma(1, shape = alpha, scale = beta[j])
  }

  output <- list(sigma2 = sigma2.initial,
                 lambda = lambda.initial,
                 x = x.initial)

  return(output)
}



proposalFun <- function(model,...) UseMethod("proposalFun", model)

proposalFun.truncatedN <-  function(model, currentVal, n, dist.mat,
                                    prevX, metric, annealingPar){

  hyperparList <- model$hyperparList
  n.obj <- n
  m <- n.obj *(n.obj - 1)/2

  # Extract values from parList
  a <- hyperparList$a
  b <- hyperparList$b
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta
  constant.multiple <- hyperparList$constant_multiple

  # Get current values
  x.cur <- currentVal$x
  sigma2.cur <- currentVal$sigma2
  lambda.cur <- currentVal$lambda
  SSR.cur <- currentVal$SSR

  ### Propose new values for parameters

  ##### lambda_j #####
  # The full conditional posterior distrbuiton of lambda_j is the inverse Gamma distribution
  # lambda_j ~ IG(alpha + n/2 , beta_j + s_j/2)
  # where s_j/n is the sample variance of the jth coordinates of x_i's
  lambda.proposal <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p){

    # calculate shape parameter
    lambda.shape <- alpha + n/2
    # calculate scale parameter
    # first calculate s_j/n
    sj.n <- var(x.cur[, j])*(n.obj-1)/n.obj
    lambda.scale <- beta[j] + sj.n*n.obj/2

    diag(lambda.proposal)[j] <- MCMCpack::rinvgamma(n = 1, shape = lambda.shape, scale = lambda.scale)

  }

  ##### x_i #####
  # A normal proposal density is used in the random walk Metropolis algorithm for
  # generation of x_i, i = 1, ... , n
  # Choose the variance of the normal proposal density to be a constant
  # multiple of sigma2/(n-1)
  x.proposal <- matrix(data = 0, nrow = n.obj, ncol = p)
  x.var <- constant.multiple * sigma2.cur/(n.obj - 1)
  for (i in 1:n.obj){
    x.proposal[i, ] <- mvtnorm::rmvnorm(n = 1, mean = x.cur[i, ], sigma = diag(rep(x.var, p)))
  }

  ##### sigma2 #####
  # A normal proposal density is used in the random walk Metropolis algorithm for
  # generation of sigma2
  # Choose the variance of the normal proposal density to be proportional to the
  # variance of IG(m/2+a, SSR/2+b)
  m <- n.obj*(n.obj-1)/2
  sigma2.sd <- sqrt(constant.multiple * (SSR.cur/2 + b)^2/((m/2+a-1)^2*(m/2+a-2)))
  sigma2.proposal <- exp(rnorm(n = 1, mean = log(sigma2.cur), sd = sigma2.sd))

  if (sigma2.proposal <=  1e-10 | sigma2.proposal == Inf){
    output <- list(lambda = lambda.cur,
                   x = x.cur,
                   sigma2 = sigma2.cur)
  } else {

    proposal <- list(lambda = lambda.proposal, x = x.proposal, sigma2 = sigma2.proposal)

    ### Determine accept or reject
    # Compute acceptance probability
    result.new <- likelihoodFun(model, dist.mat, proposal.result = proposal, metric)
    result.cur <- likelihoodFun(model, dist.mat, proposal.result = currentVal, metric)
    dproposal.new <- dproposalFun(model, n = n.obj, dist.mat,
                                  para.result.l = proposal, para.result.r = currentVal, metric)
    dproposal.cur <- dproposalFun(model, n = n.obj, dist.mat,
                                  para.result.l = currentVal, para.result.r = proposal, metric)

    probab <- exp(#result.new - result.cur +
                    annealingPar*(result.new$logposterior - result.cur$logposterior) +
                    (1-annealingPar)*logReferenceRatio(proposal$x, x.cur, prevX))

    # accept or reject step
    if (runif(1) < probab){
      # accept
      lambda.proposal <- proposal$lambda
      x.proposal <- proposal$x
      sigma2.proposal <- proposal$sigma2

    }else{
      # reject
      lambda.proposal <- lambda.cur
      x.proposal <- x.cur
      sigma2.proposal <- sigma2.cur
    }

    output <- list(lambda = lambda.proposal,
                   x = x.proposal,
                   sigma2 = sigma2.proposal)

  }

  return(output)

}



dproposalFun <- function(model,...) UseMethod("dproposalFun", model)

dproposalFun.truncatedN <- function(model, n, dist.mat,
                                              para.result.l, para.result.r, metric){

  hyperparList <- model$hyperparList
  n.obj <- n
  m <- n.obj *(n.obj - 1)/2

  # Extract values from parList
  a <- hyperparList$a
  b <- hyperparList$b
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta
  constant.multiple <- hyperparList$constant_multiple

  # Get parameter values
  x.l <- para.result.l$x
  sigma2.l <- para.result.l$sigma2
  lambda.l <- para.result.l$lambda
  x.r <- para.result.r$x
  sigma2.r <- para.result.r$sigma2
  lambda.r <- para.result.r$lambda

  # calculate delta matrix and d matrix
  d.mat.l <- dist.mat
  d.mat.r <- dist.mat
  if (dist.metric[dist.metric.index] == "cosine"){
    delta.mat.l <- 1-philentropy::distance(x.l, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat.l) =  rownames(delta.mat.l) <- rownames(d.mat.l)
    delta.mat.r <- 1-philentropy::distance(x.r, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat.r) =  rownames(delta.mat.r) <- rownames(d.mat.r)
  } else{
    delta.mat.l <- philentropy::distance(x.l, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat.l) =  rownames(delta.mat.l) <- rownames(d.mat.l)
    delta.mat.r <- philentropy::distance(x.r, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat.r) =  rownames(delta.mat.r) <- rownames(d.mat.r)
  }

  # calculate SSR
  SSR.l <- SSRFun(d.mat.l, delta.mat.l)
  SSR.r <- SSRFun(d.mat.r, delta.mat.r)

  ##### lambda_j #####
  # The full conditional posterior distrbuiton of lambda_j is the inverse Gamma distribution
  # lambda_j ~ IG(alpha + n/2 , beta_j + s_j/2)
  # where s_j/n is the sample variance of the jth coordinates of x_i's
  lambda.density.mat <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p){

    # calculate shape parameter
    shape.lambda <- alpha + n/2
    # calculate scale parameter
    # first calculate s_j/n
    sj.n <- var(x.r[, j])
    scale.lambda <- beta[j] + sj.n*n.obj/2

    diag(lambda.density.mat)[j] <- log(MCMCpack::dinvgamma(x = diag(lambda.l)[j],
                                                           shape = shape.lambda, scale = scale.lambda))

  }
  lambda.density <- sum(diag(lambda.density.mat))

  ##### x_i #####
  # A normal proposal density is used in the random walk Metropolis algorithm for
  # generation of x_i, i = 1, ... , n
  # Choose the variance of the normal proposal density to be a constant
  # multiple of sigma2/(n-1)
  x.density.mat <- numeric()
  x.var <- constant.multiple * sigma2.r/(n.obj - 1)
  for (i in 1:n.obj){
    x.density.mat[i] <- mvtnorm::dmvnorm(x = x.l[i, ], mean = x.r[i, ], sigma = diag(rep(x.var, p)),
                                         log = TRUE)
  }
  x.density <- sum(x.density.mat)

  ##### sigma2 #####
  # A normal proposal density is used in the random walk Metropolis algorithm for
  # generation of sigma2
  # Choose the variance of the normal proposal density to be proportional to the
  # variance of IG(m/2+a, SSR/2+b)
  sigma2.sd <- sqrt(constant.multiple * (SSR.r/2 + b)^2/((m/2+a-1)^2*(m/2+a-2)))
  sigma2.density <- dnorm(x = log(sigma2.l), mean = sigma2.r, sd = sigma2.sd,
                          log = TRUE)

  output <- lambda.density + x.density + sigma2.density

  return(output)

}



likelihoodFun <- function(model, ...) UseMethod("likelihoodFun", model)

likelihoodFun.truncatedN <- function(model,  dist.mat, proposal.result, metric){

  hyperparList <- model$hyperparList

  n.obj <- n
  m = n.obj *(n.obj - 1)/2

  # get xi
  # xi is the x.mat
  x.mat <- proposal.result$x

  # get sigma^2
  sigma2 <- proposal.result$sigma2

  # get lambda
  lambda <- proposal.result$lambda

  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  # calculate delta matrix and d matrix
  d.mat <- dist.mat
  if (dist.metric[dist.metric.index] == "cosine"){
    delta.mat <- 1-philentropy::distance(x.mat, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  } else{
    delta.mat <- philentropy::distance(x.mat, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  }


  # calculate SSR
  SSR <- SSRFun(d.mat, delta.mat)

  # calculate the term of sum over log of standard normal cdf
  delta.mat[!lower.tri(delta.mat)] <- NA
  div.mat <- delta.mat/sqrt(sigma2)
  log.normal.cdf <- log(pnorm(div.mat))
  sum.log.normal.cdf <- sum(log.normal.cdf[is.finite(log.normal.cdf)], na.rm = TRUE)

  ## calculate the log likelihood
  loglikelihood <- -(m/2)*log(sigma2) - SSR/(2*sigma2) - sum.log.normal.cdf

  ## calculate the log priors
  # prior for X
  x.logprior <- mvtnorm::dmvnorm(x.mat, mean = rep(0, p), sigma = lambda, log = TRUE)
  # prior for sigma2
  sigma2.logprior <- log(MCMCpack::dinvgamma(sigma2, shape = a, scale = b))
  # prior for lambda
  lambda.logprior <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p){
    diag(lambda.logprior)[j] <- log(MCMCpack::dinvgamma(diag(lambda)[j], shape = alpha, scale = beta[j]))
  }

  logprior <- sum(x.logprior) + sigma2.logprior + sum(diag(lambda.logprior))


  # calculate log posterior
  # As the sum of the log-likelihood and the log-prior.
  # This is equivalent to multiplying the probabilities and then taking the log
  logposterior <- loglikelihood + logprior

  output <- data.frame(logposterior = logposterior,
                       loglikelihood = loglikelihood,
                       logprior = logprior,
                       SSR = SSR)
  return(output)

}

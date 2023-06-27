

truncatedT_incr <- function(hyperparList, p=2, reference.x.sd){
  bmdsmodel <- list()
  bmdsmodel$hyperparList <- hyperparList
  bmdsmodel$p <- p
  bmdsmodel$reference.x.sd <- reference.x.sd
  ## Set the name for the class
  class(bmdsmodel) <- append(class(bmdsmodel), c("truncatedT_incr", "BMDSModel"))
  return(bmdsmodel)
}



prior <- function(model) UseMethod("prior", model)

prior.truncatedT_incr <- function(model){

  hyperparList <- model$hyperparList

  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  # get df
  df.nu <- hyperparList$df

  ### Initialize particles
  # sigma2
  sigma2.initial <- MCMCpack::rinvgamma(1, shape = a, scale = b)
  # lambda
  lambda.initial <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p){
    diag(lambda.initial)[j] <- MCMCpack::rinvgamma(1, shape = alpha, scale = beta[j])
  }

  # g_ij
  g.initial <- matrix(data = 0, nrow = n, ncol = n)
  g.initial[upper.tri(g.initial)] <- rgamma(n = n.obj*(n.obj - 1)/2, shape = df.nu, rate = df.nu)

  # x_i
  x.initial <- matrix(data = 0, nrow = n, ncol = p)
  for (i in 1:n){
    x.initial[i, ] <- mvtnorm::rmvnorm(1, mean = rep(0, p), sigma = lambda.initial)
  }

  output <- list(sigma2 = sigma2.initial,
                 lambda = lambda.initial,
                 x = x.initial,
                 g = g.initial)

  return(output)

}


initialFun <- function(model, ...) UseMethod("initialFun", model)

initialFun.truncatedT_incr <- function(model, prev.result, dist.mat, metric){

  hyperparList <- model$hyperparList

  n.obj <- nrow(dist.mat)
  m <- n.obj*(n.obj - 1)/2

  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  # get df
  df.nu <- hyperparList$df

  ### Initialize particles
  ## x_i (from reference distribution)
  x.prev <- prev.result$xi
  x.prev.ref <- t(apply(x.prev, 1, reference$x))
  n.incr <- nrow(dist.mat) - nrow(x.prev)
  x.incr <- reference$x.incr(n.incr = n.incr)
  x.initial <- rbind(x.prev.ref, x.incr)
  rownames(x.initial) <- rownames(dist.mat)

  ## sigma2 (from prior distribution)
  sigma2.initial <- MCMCpack::rinvgamma(1, shape = a, scale = b)
  ## lambda (from prior distribution)
  lambda.initial <- matrix(data = 0, nrow = p, ncol = p)
  for (j in 1:p)
    diag(lambda.initial)[j] <- MCMCpack::rinvgamma(1, shape = alpha, scale = beta[j])

  ## g (from prior distribution)
  g.initial <- matrix(data = 0, nrow = n.obj, ncol = n.obj)
  g.initial[upper.tri(g.initial)] <- rgamma(n = n.obj*(n.obj - 1)/2, shape = df.nu, rate = df.nu)

  output <- list(sigma2 = sigma2.initial,
                 lambda = lambda.initial,
                 x = x.initial,
                 g = g.initial)

  return(output)

}



proposalFun <- function(model,...) UseMethod("proposalFun", model)

proposalFun.truncatedT_incr <- function(model, currentVal, n, dist.mat,
                                           prevX, metric, annealingPar){

  hyperparList <- model$hyperparList
  n.obj <- n
  m <- n.obj *(n.obj - 1)/2

  # Extract values from parList
  a <- hyperparList$a
  b <- hyperparList$b
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta
  df.nu <- hyperparList$df
  constant.multiple <- hyperparList$constant.multiple

  # Get current values
  x.cur <- currentVal$x
  sigma2.cur <- currentVal$sigma2
  lambda.cur <- currentVal$lambda
  SSR.cur <- currentVal$SSR
  g.cur <- currentVal$g

  # calculate delta matrix and d matrix
  d.mat <- as.matrix(dist.mat)
  if (dist.metric[dist.metric.index] == "cosine"){
    delta.mat <- 1-philentropy::distance(x.cur, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  } else{
    delta.mat <- philentropy::distance(x.cur, method = dist.metric[dist.metric.index], mute.message = TRUE)
    colnames(delta.mat) =  rownames(delta.mat) <- rownames(d.mat)
  }

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

  ##### g_ij #####
  # g_ij condition on other variables follows a Gamma distribution
  # g_ij ~ Gamma((1+df)/2 , (delta_ij-d_ij)^2/(2*sigma2) + df/4)
  shape.par <- (1+df.nu)/2
  g.proposal <- matrix(data = 0, nrow = n.obj, ncol = n.obj)
  for (i in 1:n.obj){
    for (j in 1:n.obj){
      d.ij <- d.mat[i, j]
      delta.ij <- delta.mat[i, j]
      g.proposal[i, j] <- rgamma(n = 1, shape = shape.par,
                                 rate = (d.ij - delta.ij)^2/(2*sigma2.cur) + df.nu/4)
    }
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
                   sigma2 = sigma2.cur,
                   g = g.cur)
  } else {

    proposal <- list(lambda = lambda.proposal, x = x.proposal,
                     sigma2 = sigma2.proposal, g = g.proposal)

    ### Determine accept or reject
    # Compute acceptance probability
    result.new <- likelihoodFun(model, dist.mat, proposal.result = proposal, metric, n.incr)
    result.cur <- likelihoodFun(model, dist.mat, proposal.result = currentVal, metric, n.incr)
    dproposal.new <- dproposalFun(model, n = n.obj, dist.mat,
                                        para.result.l = proposal, para.result.r = currentVal, metric)
    dproposal.cur <- dproposalFun(model, n = n.obj, dist.mat,
                                        para.result.l = currentVal, para.result.r = proposal, metric)

    probab <- exp(annealingPar*(result.new$logposterior - result.cur$logposterior) +
                    (1-annealingPar)*logReferenceRatio(proposal$x, x.cur, prevX))

    # accept or reject step
    if (runif(1) < probab){
      # accept
      lambda.proposal <- proposal$lambda
      x.proposal <- proposal$x
      sigma2.proposal <- proposal$sigma2
      g.proposal <- proposal$g

    }else{
      # reject
      lambda.proposal <- lambda.cur
      x.proposal <- x.cur
      sigma2.proposal <- sigma2.cur
      g.proposal <- g.cur
    }

    output <- list(lambda = lambda.proposal,
                   x = x.proposal,
                   sigma2 = sigma2.proposal,
                   g = g.proposal)

  }

  return(output)

}



dproposalFun <- function(model,...) UseMethod("dproposalFun", model)

dproposalFun.truncatedT_incr <- function(model, n, dist.mat,
                                            para.result.l, para.result.r, metric){

  hyperparList <- model$hyperparList
  n.obj <- n
  m <- n.obj *(n.obj - 1)/2

  # Extract values from parList
  a <- hyperparList$a
  b <- hyperparList$b
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta
  df.nu <- hyperparList$df
  constant.multiple <- hyperparList$constant.multiple

  # Get parameter values
  x.l <- para.result.l$x
  sigma2.l <- para.result.l$sigma2
  lambda.l <- para.result.l$lambda
  x.r <- para.result.r$x
  sigma2.r <- para.result.r$sigma2
  lambda.r <- para.result.r$lambda
  g.l <- para.result.l$g
  g.r <- para.result.r$g

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
  #m <- n.obj*(n.obj-1)/2
  sigma2.sd <- sqrt(constant.multiple * (SSR.r/2 + b)^2/((m/2+a-1)^2*(m/2+a-2)))
  sigma2.density <- dnorm(x = log(sigma2.l), mean = sigma2.r, sd = sigma2.sd,
                          log = TRUE)

  ##### g #####
  # g_ij condition on other variables follows a Gamma distribution
  # g_ij ~ Gamma((1+df)/2 , (delta_ij-d_ij)^2/(2*sigma2) + df/4)
  shape.par <- (1+df.nu)/2
  g.density.mat <- matrix(data = NA, nrow = n.obj, ncol = n.obj)
  for (i in 1:n.obj){
    for (j in 1:n.obj){
      g.density.mat[i, j] <- dgamma(g.l[i, j], shape = shape.par,
                                    rate = (d.mat.l[i, j] - delta.mat.l[i, j])^2/(2*sigma2.l) + df.nu/4,
                                    log = TRUE)
    }
  }
  g.density <- sum(g.density.mat[upper.tri(g.density.mat)])


  output <- lambda.density + x.density + sigma2.density + g.density

  return(output)

}



likelihoodFun <- function(model, ...) UseMethod("likelihoodFun", model)

likelihoodFun.truncatedT_incr <- function(model, dist.mat, proposal.result, metric, n.incr){

  hyperparList <- model$hyperparList

  n.obj <- n
  m <- n.obj *(n.obj - 1)/2
  #n.incr.obj <- n.incr

  # get xi
  # xi is the x.mat
  x.mat <- proposal.result$x

  # get sigma^2
  sigma2 <- proposal.result$sigma2

  # get lambda
  lambda <- proposal.result$lambda

  # get g
  g <- proposal.result$g

  # get a, b
  a <- hyperparList$a
  b <- hyperparList$b

  # get alpha, beta
  alpha <- hyperparList$alpha
  beta <- hyperparList$beta

  # get df
  df.nu <- hyperparList$df

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
  # calculate SSR * g
  d.upper.mat <- d.mat[upper.tri(d.mat)]
  delta.upper.mat <- delta.mat[upper.tri(delta.mat)]
  g.upper.mat <- g[upper.tri(g)]
  SSR.g <- sum(g.upper.mat * (d.upper.mat-delta.upper.mat)^2)

  # calculate the term of sum over log of standard normal cdf
  div.mat <- delta.upper.mat*sqrt(g.upper.mat)/sqrt(sigma2)
  log.normal.cdf <- log(pnorm(div.mat))
  sum.log.normal.cdf <- sum(log.normal.cdf, na.rm = TRUE)

  ## calculate the log likelihood
  loglikelihood <- -(m/2)*log(sigma2) - SSR.g/(2*sigma2) - sum.log.normal.cdf + 0.5*sum(log(g.upper.mat))

  ## calculate the log priors
  # prior for X_incr
  x.prev.mat <- x.mat[1:(n.obj-n.incr),]
  x.logprior <- mvtnorm::dmvnorm(x.prev.mat, mean = rep(0, p), sigma = lambda, log = TRUE)
  # # prior for sigma2
  # sigma2.logprior <- log(MCMCpack::dinvgamma(sigma2, shape = a, scale = b))
  # # prior for lambda
  # lambda.logprior <- matrix(data = 0, nrow = p, ncol = p)
  # for (j in 1:p){
  #   diag(lambda.logprior)[j] <- log(MCMCpack::dinvgamma(diag(lambda)[j], shape = alpha, scale = beta[j]))
  # }
  # # prior for g
  # g.logprior <- matrix(data = 0, nrow = n.obj, ncol = n.obj)
  # for (i in 1:n.obj){
  #   for (j in 1:n.obj){
  #     g.logprior[i, j] <- dgamma(g[i, j], shape = df.nu, rate = df.nu, log = TRUE)
  #   }
  # }
  #
  # logprior <- sum(x.logprior) + sigma2.logprior + sum(diag(lambda.logprior)) + sum(g.logprior[upper.tri(g.logprior)])
  logprior <- x.logprior

  # calculate log posterior
  # As the sum of the log-likelihood and the log-prior.
  # This is equivalent to multiplying the probabilities and then taking the log
  logposterior <- loglikelihood + logprior

  output <- list(logposterior = logposterior,
                 loglikelihood = loglikelihood,
                 logprior = logprior,
                 SSR = SSR)

  return(output)

}

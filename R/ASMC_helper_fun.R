#' Calculate the stress value
#'
#' @param d.mat  distance matrix
#' @param delta.mat  delta matrix
#' @return the stress value
#' @examples
#' print("TODO: add an example")
stressFun <- function(d.mat, delta.mat){

  d.mat[!lower.tri(d.mat)] <- 0
  delta.mat[!lower.tri(delta.mat)] <- 0

  stress <- sqrt(sum((d.mat-delta.mat)^2)/(sum(d.mat^2)))

  return(STRESS = stress)

}

#' Calculate the SSR value
#'
#' @param d.mat  distance matrix
#' @param delta.mat  delta matrix
#' @return the SSR value
#' @examples
#' print("TODO: add an example")
SSRFun <- function(d.mat, delta.mat){

  d.mat[!lower.tri(d.mat)] <- 0
  delta.mat[!lower.tri(delta.mat)] <- 0

  ssr <- sum((d.mat-delta.mat)^2)

  return(SSR = ssr)

}

#' Calculate the ln(SSR) value
#'
#' @param d.mat  distance matrix
#' @param delta.mat  delta matrix
#' @return the ln(SSR) value
#' @examples
#' print("TODO: merge this function with the above function, and add an example")
SSR.ln.fun <- function(d.mat, delta.mat){

  upper.d.mat <- d.mat[upper.tri(d.mat)]
  upper.delta.mat <- delta.mat[upper.tri(delta.mat)]

  upper.d.mat[which(upper.d.mat == 0)] <- 1e-8

  ssr.ln <- sum((log(upper.d.mat)-upper.delta.mat)^2)

  return(SSR.ln = ssr.ln)

}




NextAnnealingParameter <- list()

#' Calculate relative conditional effective sample size
#'
#' @param W  particle weights
#' @param logL log likelihood
#' @param num  evaluated point
#' @param phi tuning parameter that controls the length of the annealing sequence
#' @return relative conditional effective sample size
#' @examples
#' print("TODO: add an example")
rCESSFun <- function(W, logL, num, phi){
  logw <- num*logL
  logmax <- max(logw)
  w <- exp(logw - logmax)
  rt <- (sum(w*W))^2/(sum(W*w^2)) - phi
  return(rt)
}


#' Calculate the next annealing parameter
#'
#' @param low  lower bound of the search interval
#' @param high upper bound of the search interval
#' @param W  particle weights
#' @param logL log likelihood
#' @param phi tuning parameter that controls the length of the annealing sequence
#' @return the incremental change of the two successive annealing parameters
#' @examples
#' print("TODO: add an example")
bisectionFun <- function(low, high, W, logL, phi) {
  n = 1000
  tol = 1e-7
  g.low <- rCESSFun(W, logL, low, phi)
  g.high <- rCESSFun(W, logL, high, phi)
  if (is.na(g.low*g.high)){
    stop("Bisection flawed")
  }
  if (g.low*g.high > 0){
    return(99e-3)
  }

  for (i in 1:n) {
    mid <- (low + high)/2
    g.mid <- rCESSFun(W, logL, mid, phi)
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if ((g.mid == 0) || ((high - low) / 2) < tol) {
      return(mid)
    }

    # If another iteration is required,
    # check the signs of the function at the points
    # reassign low or high accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(g.mid) == sign(g.low),
           low <- mid,
           high <- mid)
  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}



reference <- list()

#' Generate x from the reference distribution
#'
#' @param x  current values of x
#' @return proposed values of x from the reference distribution
#' @examples
#' print("TODO: add an example")
reference$x <- function(x){

  x.new <- mvtnorm::rmvnorm(1, mean = x, sigma = reference.x.sd)
  return(x.new)

}

#' Generate x from the reference distribution in the case of incremental inference
#'
#' @param n.incr the number of incremental dimensions of x
#' @return proposed values of the incremental part of x from the reference distribution
#' @examples
#' print("TODO: add an example")
reference$x.incr <- function(n.incr){

  x.incr <- mvtnorm::rmvnorm(n.incr, mean = rep(0, p), sigma = reference.x.sd)
  return(x.incr)

}

#' Calculate the density of x
#'
#' @param x current values of x
#' @param prev.x previous values of x
#' @return density of x
#' @examples
#' print("TODO: add an example")
reference$d.x <- function(x, prev.x){

  x.list <- lapply(seq_len(nrow(x)), function(i) x[i,])
  prev.x.list <- lapply(seq_len(nrow(prev.x)), function(i) prev.x[i,])

  reference.x.sd.list <- rep(list(reference.x.sd), n)

  d.x <- mapply(mvtnorm::dmvnorm, x.list, prev.x.list, reference.x.sd.list, log = TRUE)

  return(d.x)

}

class(reference) <- append(class(reference), c("referenceDistribution"))


# logReferenceRatio <- function(new.x, current.x, prev.x){
#
#   ratio <- sum(reference$d.x(new.x, current.x)) - sum(reference$d.x(current.x, prev.x))
#
#   return(ratio)
#
# }

#' Calculate the reference ratio in log scale
#'
#' @param new.x proposed values of x
#' @param current.x current values of x
#' @param prev.x previous values of x
#' @return log reference ratio
#' @examples
#' print("TODO: add an example")
logReferenceRatio <- function(new.x, current.x, prev.x){

  ratio <- sum(reference$d.x(new.x, prev.x)) - sum(reference$d.x(current.x, prev.x))

  return(ratio)

}




#' Calculate relative effective sample size
#'
#' @param logW  particle weights in log scale
#' @return relative effective sample size
#' @examples
#' print("TODO: add an example")
rESSFun <- function(logW){

  K <- length(logW)
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(K)

  return(exp(logRESS))

}


#' Obtain resampling index using multinomial resampling
#'
#' @param W  particle weights
#' @return index for particle after resampling
#' @examples
#' print("TODO: add an example")
multinomialResampleFun <- function(W){

  # resample particles using the normalized weights
  N <- length(W)
  index <- sample(1:N, size = N, replace = T, prob = W)

  return(index)

}

# log-sum-exponential evaluation of log(sum(w))
logsum <- function(logw){
  logmax = max(logw)
  log(sum(exp(logw-logmax)))+logmax
}


jaccardDist <- function(input){
  distmatrix <- stringdist::stringdistmatrix(tolower(input),
                                 useNames = FALSE, method = "jaccard")
  dis <- as.matrix(distmatrix)
  
  return(dis)
  
}




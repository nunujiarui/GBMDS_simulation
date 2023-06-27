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
library(ggrepel)

#' Plot CMDS result
#'
#' @param cmds.result a matrix with up to p columns whose rows give the coordinates of the points chosen to represent the dissimilarities from CMDS.
#' @param data.label labels of the points
#' @param plot.title title of the plot
#' @param rotate indicates whether rotation should be applied
#' @param rotate.angle degree counter clockwise rotation
#' @param flip indicates whether flipping should be applied
#' @param flip.h indicates whether flip is horizontal
#' @param flip.v indicates whether flip is vertical
#' @return a 2D scatterplot of points with labels
#' @examples
#' print("TODO: add an example")
plot.CMDS <- function(cmds.result, data.label,
                      plot.title = "",
                      rotate = FALSE, rotate.angle = NULL,
                      flip = FALSE, flip.h = FALSE, flip.v = FALSE){

  result <- cmds.result$points

  if (!(rotate | flip)){
    ## no rotation/flipping applied
    ## Make the plot
    cmds.xy.data <- data.frame(result, data.label)
    colnames(cmds.xy.data) <- c("X", "Y", "label")
    plot <- ggplot() +
      geom_point(data = cmds.xy.data, aes(x = X, y = Y, colour = data.label)) +
      geom_label_repel(data = cmds.xy.data, aes(x = X, y = Y, label = label)) +
      theme_bw() +
      labs(x = "X coordinate", y = "Y coordinate",
           title = plot.title) +
      theme(legend.position = "None")
  }

  if (rotate | flip){
    if (rotate){
      ## rotate the points
      t.matrix <- matrix(data = c(cos(rotate.angle), -sin(rotate.angle), sin(rotate.angle), cos(rotate.angle)),
                         nrow = 2, byrow = TRUE)
      result <- t(t.matrix %*% t(result))
    }
    if (flip){
      ## flip the points
      if (flip.h){
        result[, 1] <- result[, 1]*(-1)
      }
      if (flip.v){
        result[, 2] <- result[, 2]*(-1)
      }
    }
    cmds.xy.data <- data.frame(result, data.label)
    colnames(cmds.xy.data) <- c("X", "Y", "label")
    plot <- ggplot() +
      geom_point(data = cmds.xy.data, aes(x = X, y = Y, colour = data.label)) +
      geom_label_repel(data = cmds.xy.data, aes(x = X, y = Y, label = label)) +
      theme_bw() +
      labs(x = "Transformed X coordinate", y = "Transformed Y coordinate",
           title = plot.title) +
      theme(legend.position = "None")
  }

  print(plot)

}





#' Plot ASMC result
#'
#' @param asmc.result an object with class name "BMDSParticles" from ASMC() function
#' @param data.label labels of the points
#' @param ci.level confidence level
#' @param plot.title title of the plot
#' @param procruste.transform indicates whether Procruste transform should be applied
#' @param rotate indicates whether rotation should be applied
#' @param rotate.angle degree counter clockwise rotation
#' @param flip indicates whether flipping should be applied
#' @param flip.h indicates whether flip is horizontal
#' @param flip.v indicates whether flip is vertical
#' @return a 2D scatterplot of points with labels
#' @examples
#' print("TODO: add an example")
plot.BMDSParticles <- function(asmc.result, data.label, ci.level = 0.95,
                               plot.title = "",
                               procruste.transform = TRUE,
                               rotate = FALSE, rotate.angle = NULL,
                               flip = FALSE, flip.h = FALSE, flip.v = FALSE){
  # posterior inference
  index.asmc <- which.min(asmc.result$SSR.output)
  asmc.res <- asmc.result$xi.output[[index.asmc]]

  # resample so that all particles have the same weight
  # construct the confidence region from these particles
  particle.index <- multinomialResampleFun(asmc.result$weight.output)
  asmc.result$xi.output.resample <- list()
  for (i in 1:length(particle.index)){
    asmc.result$xi.output.resample[[i]] <- asmc.result$xi.output[[particle.index[i]]]
  }

  if (procruste.transform){
    ## perform the procruste transformation
    x.asmc.ptransf.all <- lapply(asmc.result$xi.output.resample, MCMCpack::procrustes, Xstar = asmc.res)
    # extract the transformed matrices
    x.asmc.ptransf <- lapply(x.asmc.ptransf.all , `[[`, 1)
  } else{
    ## do NOT perform the Procruste transformation
    # extract the transformed matrices
    x.asmc.ptransf <- asmc.result$xi.output.resample
  }

  if (rotate | flip){
    if (rotate){
      ## rotate the points
      t.matrix <- matrix(data = c(cos(rotate.angle), -sin(rotate.angle), sin(rotate.angle), cos(rotate.angle)),
                         nrow = 2, byrow = TRUE)
      x.asmc.ptransf <- lapply(x.asmc.ptransf, function(x) t(t.matrix %*% t(x)))
    }
    if (flip){
      ## flip the points
      if (flip.h){
        x.asmc.ptransf <- lapply(x.asmc.ptransf, function(x) t(t(x) * c(-1,1)))
      }
      if (flip.v){
        x.asmc.ptransf <- lapply(x.asmc.ptransf, function(x) t(t(x) * c(1,-1)))
      }
    }
  }

  # reshape the data
  bind.ith.rows.asmc.ptransf <- function(i) do.call(rbind, lapply(x.asmc.ptransf, "[", i, TRUE))
  x.posterior.asmc.ptransf = lapply(1:n, bind.ith.rows.asmc.ptransf)
  x.credible.region.asmc.ptransf <- list()
  x.credible.region.asmc.df.ptransf <- list()
  for (i in 1:length(x.posterior.asmc.ptransf)){
    x.asmc.all.ptransf <- as.matrix(x.posterior.asmc.ptransf[[i]])
    ellipse.coor.ptransf <- car::dataEllipse(x.asmc.all.ptransf[, 1], x.asmc.all.ptransf[, 2],
                                             levels = ci.level, draw = FALSE)
    x.credible.region.asmc.ptransf[[i]] <- ellipse.coor.ptransf
    x.credible.region.asmc.df.ptransf[[i]] <- data.frame(ellipse.coor.ptransf)
    colnames(x.credible.region.asmc.df.ptransf[[i]]) <- c("X", "Y")
  }

  # get the posterior median as our estimation
  X.median <- matrix(data = NA, nrow = n, ncol = p)
  for (i in 1:n){
    X.median[i, ] <- apply(x.posterior.asmc.ptransf[[i]] , 2, median)
  }

  # plot the estimated coordinates and its credible region after procrustes transformation
  x.posterior.mode.df.asmc.ptransf <- data.frame(X.median)
  colnames(x.posterior.mode.df.asmc.ptransf) <- c("X", "Y")
  x.posterior.mode.df.asmc.ptransf$id <- data.label
  names(x.credible.region.asmc.df.ptransf) <- data.label
  plot.ellipse.asmc.ptransf <- dplyr::bind_rows(x.credible.region.asmc.df.ptransf, .id = "id")

  id <- 1:length(data.label)
  x.label <- ifelse(procruste.transform | rotate | flip, "Transformed X coordinate",
                    "X coordinate")
  y.label <- ifelse(procruste.transform | rotate | flip, "Transformed Y coordinate",
                    "Y coordinate")
  ggplot() +
    geom_point(data = x.posterior.mode.df.asmc.ptransf, aes(x = X, y = Y, color = id)) +
    geom_label_repel(data = x.posterior.mode.df.asmc.ptransf, aes(x = X, y = Y, label = id)) +
    geom_path(data = plot.ellipse.asmc.ptransf, aes(x = X, y = Y, group = id, color = id)) +
    labs(x = x.label, y = y.label,
         title = plot.title) +
    theme_bw() +
    theme(legend.position = "None")

}











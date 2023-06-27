# This script contains codes for real data application
# This example is for section 5.3 NIH text data

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
library(tidytext)
library(textmineR)
library(text2vec)
library(philentropy)

set.seed(1402)

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

dist.metric <- c("manhattan", "euclidean", "chebyshev", "cosine")
dist.metric.index <- 4

## read in data
# load nih_sample data set from textmineR
data(nih_sample)

# create a document term matrix
dtm <- CreateDtm(doc_vec = nih_sample$ABSTRACT_TEXT, # character vector of documents
                 doc_names = nih_sample$APPLICATION_ID, # document names
                 ngram_window = c(1, 2), # minimum and maximum n-gram length
                 stopword_vec = c(stopwords::stopwords("en"), # stopwords from tm
                                  stopwords::stopwords(source = "smart")), # this is the default value
                 lower = TRUE, # lowercase - this is the default value
                 remove_punctuation = TRUE, # punctuation - this is the default
                 remove_numbers = TRUE, # numbers - this is the default
                 verbose = FALSE, # Turn off status bar for this demo
                 cpus = 2) # default is all available cpus on the system

# construct the matrix of term counts to get the IDF vector
tf_mat <- TermDocFreq(dtm)

# TF-IDF and cosine similarity
tfidf <- t(dtm[ , tf_mat$term ]) * tf_mat$idf
tfidf <- t(tfidf)
csim <- tfidf / sqrt(rowSums(tfidf * tfidf))
csim <- csim %*% t(csim)
cdist <- as.dist(1 - csim)
dis <- as.matrix(cdist)

# # small data
# data <- readRDS(file = "/data/NIPS_data_small.rds")
# dis <- 1-philentropy::distance(t(data), method = dist.metric[dist.metric.index])
# dis.df <- data.frame(dis = dis[upper.tri(dis)])

ggplot(dis.df, aes(x = dis)) +
  geom_histogram(bins = 40) +
  labs(x = "Cosine dissimilarity") +
  theme_bw()


## general settings
p <- 2
n <- nrow(dis)

cmds.result <- cmdscale(d = dis, k = p,
                        eig = TRUE, add = FALSE, x.ret = FALSE)
class(cmds.result) <- append(class(cmds.result), "CMDS")


## set hyperparameters
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

tuningparList <- list(K = 100, phi = 0.8, eps = 0.5)
n.core <- detectCores()-1

## run GBMDS-ASMC
model <- truncatedN(hyperparList, p, reference.x.sd)
#model <- truncatedT(hyperparList, p, reference.x.sd)
#model <- truncatedSkewedN(hyperparList, p, reference.x.sd)


start.time <- Sys.time()
asmc.result <- ASMC(model = model, dist.mat = dis,
                    tuningparList, n.core, cmds.result = cmds.result$points,
                    metric = dist.metric[dist.metric.index])
end.time <- Sys.time()
end.time - start.time


# log marginal likelihood estimate
asmc.result$logZ


# posterior inference
index.asmc <- which.min(asmc.result$SSR.output)
asmc.res <- asmc.result$xi.output[[index.asmc]]

## compute STRESS values
cmds.stress <- stressFun(d.mat = dis,
                         delta.mat = philentropy::distance(cmds.result$points, method = "euclidean"))
asmc.stress <- stressFun(d.mat = dis,
                         delta.mat = philentropy::distance(asmc.res, method = "euclidean"))

cat("STRESS value for CMDS is", round(cmds.stress, 4))
cat("STRESS value for BMDS with ASMC is", round(asmc.stress, 4))



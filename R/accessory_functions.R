## load some libraries
if(!require(tidyverse)) install.packages('tidyverse'); library(tidyverse)
if(!require(magrittr)) install.packages('magrittr'); library(magrittr)
if(!require(reshape2)) install.packages('reshape2'); library(reshape2)
if(!require(RColorBrewer)) install.packages('RColorBrewer'); library(RColorBrewer)
if(!require(viridis)) install.packages('viridis'); library(viridis)

##---
## functions we'll need for later
##---
## standardize columns func
norm.mat <- function(mat){
  ms <- apply(mat, 2, function(col) (col - mean(col)) / var(col))
  ms
}

## compute euclid dist and standrdize to similarities
get.euc <- function(mat){
  eu <- dist(mat, method = 'euclidean') %>% as.matrix()
  ieu <- 1 / eu
  diag(ieu) <- 0
  ieu
}

## function to threshold the normalized distance mat
threshold <- function(mat, top_k = 10){
  
  thr <- mat
  tnr <- nrow(thr)
  
  ## set similarities outside of the top k to 0
  for(i in 1:tnr){
    thr[i, !order(thr[i, ], decreasing = T) %in% 1:top_k] <- 0
  }
  
  for(i in 1:tnr){
    for(j in 1:tnr){
      if(thr[i, j] == 0 & thr[j, i] != 0){
        thr[i, j] <- thr[j, i]
      }
    }
  }
  
  thr
  
}

## function to calculate norm laplacian
get.laplac <- function(mat){
  
  L <- -mat
  S <- rowSums(mat)
  nL <- L / S
  diag(nL) <- 1
  nL
  
}


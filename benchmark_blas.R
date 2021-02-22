#!/usr/bin/Rscript

sessionInfo()
set.seed(42); x <- matrix(rnorm(100000), nrow=200)
x[1:10]
prcomp(x, scale=TRUE)$x[1:5, 1:5]
prcomp(x, scale=TRUE)$x[1:5, 1:5]
prcomp(x, scale=TRUE)$x[1:5, 1:5]



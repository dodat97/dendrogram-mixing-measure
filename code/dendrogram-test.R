library(mclust)

n = 10000
X = c(rnorm(n%/%2), rnorm(n%/%2) + 3)

mc <- Mclust(X, modelNames = 'E', G=4)
mc$parameters$pro
mc$parameters$mean
# mc$parameters$variance

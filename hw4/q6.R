############################################################ 
# Filename: q6.R
# Purpose : 113-2 Statistical Computation and Simulation, HW4, Q6
# Author  : Potumas Liu
# Date    : 2025/05/14
############################################################ 

library(MASS)
library(parallel)

n <- 1000000
nIter <- 100
nrv <- 5
th <- 21.6

# ensure reproducibility when using parallel computing
RNGkind("L'Ecuyer-CMRG")
seed <- 2025
set.seed(seed)

c1 <- makeCluster(12)
clusterSetRNGStream(c1, seed)
clusterExport(c1, ls())


# monte carlo
mcRes <- parSapply(c1, 1:nIter, FUN = function (i) {
        simX <- sapply(1:nrv, FUN = function (j) { return(j*rexp(n)) })
        simX <- rowSums(simX)
        return(mean(simX >= th))
    })
mcRes.mu <- mean(mcRes); mcRes.var <- var(mcRes)
m <- "Monte Carlo"
cat(m, ", estimate: ", mcRes.mu, ", variance: ", mcRes.var, "\n", sep="")


# important sampling
# assume p(x) is a half normal: f(x) = sqrt(2/pi) * exp{-x^2/2}
isRes <- parSapply(c1, 1:nIter, FUN = function (i) {
        #simX <- sapply(1:nrv, function (j) { return(abs(rnorm(n))) })
        simX <- sapply(1:nrv, function (j) { return(rexp(n, rate=0.8)) })
        u <- rep(0, n); for (k in 1:nrv) u <- u + k*simX[, k]
        exppdf  <- apply(dexp(simX), 1, prod)
        exp2pdf <- apply(dexp(simX, rate=0.8), 1, prod)
        #normpdf <- apply(2*dnorm(simX), 1, prod)
        return( mean((u >= th) * exppdf / exp2pdf) )
        #return( mean((u >= th) * exppdf / normpdf) )
    })
isRes.mu <- mean(isRes); isRes.var <- var(isRes)
m <- "Important Sampling"
cat(m, ", estimate: ", isRes.mu, ", variance: ", isRes.var, "\n", sep="")



# Control Variable
cvRes <- parSapply(c1, 1:nIter, FUN = function (i) {
        simU <- sapply(1:nrv, function (j) { return(runif(n)) })
        simExp  <- -log(simU)
        simX <- rep(0, n); for (k in 1:nrv) simX <- simX + k*simExp[, k]
        simX <- simX >= th

        yth <- 1
        simY <- apply(simExp < yth, 1, prod)
        a <- 1 / ((1-exp(-1))^5 - (1-exp(-1))^10)

        return( mean(simX - (0.1)*(simY-(1-exp(-1))^5)) )
    })
cvRes.mu <- mean(cvRes); cvRes.var <- var(cvRes)
m <- "Control Variable"
cat(m, ", estimate: ", cvRes.mu, ", variance: ", cvRes.var, "\n", sep="")




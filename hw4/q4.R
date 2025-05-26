############################################################ 
# Filename: q4.R
# Purpose : 113-2 Statistical Computation and Simulation, HW4, Q4
# Author  : Potumas Liu
# Date    : 2025/05/05
############################################################ 

library(parallel)

n <- 1000000
nIter <- 1000
xs <- c(-2, -2.5, -3, -3.5, -4, -5, -6) 

RNGkind("L'Ecuyer-CMRG")
seed <- 2025
set.seed(seed)

c1 <- makeCluster(12)
clusterSetRNGStream(c1, seed)
clusterExport(c1, ls())


# list theoretical F(x)
for (x in xs) cat("pnorm(", x, ") = ", pnorm(x), "\n", sep="")


# original monte carlo
mcRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            return(mean(rnorm(n) < x)) 
        })
        m <- "Monte Carlo"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        return(c(var(simMean), mean(simMean)))
    })
mcRes


# important sampling
# assume p(x) is a double exponential, lambda=1
isRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            simExp <- rexp(n) 
            simX <- simExp - x
            simX <- (-x<simX) * dnorm(simX) / dexp(simExp)
            return(mean(simX)) 
        })
        m <- "Important Sampling"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        return(c(var(simMean), mean(simMean)))
    })
isRes


# antithetic variable
avRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            simX1 <- rnorm(n)
            simX2 <- -simX1
            return( mean((x > simX1) + (x > simX2)) / 2 )
        })
        m <- "Antithetic Variable"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        return(c(var(simMean), mean(simMean)))
    })
avRes


# control variate
cvRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            simX <- rnorm(n)
            simX1 <- simX < x
            simX2 <- simX < -1.96
            E.X1 <- ifelse(x <= -3, 0.01, 0.03)
            a <- E.X1 / 0.025
            return( mean( simX1 - a*(simX2-0.025) ) )
        })
        m <- "Control Variable"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        return(c(var(simMean), mean(simMean)))
    })
cvRes

stopCluster(c1)

# export results
resMtx <- cbind(t(mcRes), t(isRes), t(avRes), t(cvRes))
colnames(resMtx) <- c("mc.var", "mc.mean", "ip.var", "ip.mean",
                      "av.var", "av.mean", "cv.var", "cv.mean")
rownames(resMtx) <- xs
fileName <- paste0("_n-", n, "_runs-", nIter, "_seed-", seed)
write.csv(resMtx, paste0("q4-results", fileName, ".csv"))

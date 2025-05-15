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
#mcRes <- parSapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            #m <- "Monte Carlo"
            #cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
            return(mean(rnorm(n) < x)) 
        })
        m <- "Monte Carlo"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        return(c(var(simMean), mean(simMean)))
    })
mcRes


# important sampling
# assume p(x) is a double exponential, lambda=1
#isRes <- sapply(xs, FUN = function(x) {
isRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            #m <- "Important Sampling"
            #cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
            #simX <- rexp(n) * ((runif(n)>0.5)*2-1)
            #simX <- (x>simX) * dnorm(simX) / (dexp(abs(simX))/2)
            simExp <- rexp(n) 
            simX <- simExp - x
            simX <- (-x<simX) * dnorm(simX) / dexp(simExp)
            return(mean(simX)) 
        })
        m <- "Important Sampling"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        #return(var(simMean))
        return(c(var(simMean), mean(simMean)))
    })
isRes


# antithetic variable
#avRes <- sapply(xs, FUN = function(x) {
avRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            #m <- "Antithetic Variable"
            #cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
            simX1 <- rnorm(n)
            simX2 <- -simX1
            return( mean((x > simX1) + (x > simX2)) / 2 )
        })
        m <- "Antithetic Variable"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        #return(var(simMean))
        return(c(var(simMean), mean(simMean)))
    })
avRes


# control variate
#cvRes <- sapply(xs, FUN = function(x) {
cvRes <- sapply(xs, FUN = function(x) {
        simMean <- parSapply(c1, 1:nIter, FUN = function(i){
            #m <- "Control Variable"
            #cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
            # generate 2 iid Unif variable s.t. U^2+V^2<1
            #U <- runif(n*1.4, min=-1); V <- runif(n*1.4, min=-1)
            #S <- U^2 + V^2
            #idx <- which(S < 1)
            #if (length(idx) < n) stop()
            #idx <- idx[1:n]
            #U <- U[idx]; S <- S[idx]
            #simX1 <- ( U * sqrt((-2*log(S))/S) ) < x
            #simX2 <- U < -0.94
            #a <- -0.2 / 0.0291
            simX <- rnorm(n)
            simX1 <- simX < x
            simX2 <- simX < -1.96
            E.X1 <- ifelse(x <= -3, 0.01, 0.03)
            a <- E.X1 / 0.025
            return( mean( simX1 - a*(simX2-0.025) ) )
        })
        m <- "Control Variable"
        cat(m, ", x = ", x, ", pnorm(", x, ") = ", mean(simMean), ", s.e. = ", var(simMean), "\n", sep="")
        #return(var(simMean))
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

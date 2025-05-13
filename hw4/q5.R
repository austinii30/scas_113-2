############################################################ 
# Filename: q5.R
# Purpose : 113-2 Statistical Computation and Simulation, HW4, Q5
# Author  : Potumas Liu
# Date    : 2025/05/12
############################################################ 

library(MASS)
library(parallel)

n <- 1000000
nIter <- 50
rhos <- c(0.7, 0.3, 0, -0.5, -0.9)
ks <- c(0, 1, 2, 3, 4)
mu1 <- mu2 <- 0
Sigma <- matrix(c(1, 0, 0, 1), nrow=2)

set.seed(2025)

c1 <- makeCluster(12)
clusterExport(c1, ls())


# original monte carlo
#mcRes <- sapply(xs, FUN = function(x) {
Res <- parLapply(c1, rhos, fun = function(rho) {
        library(MASS)
        kRes <- sapply(ks, FUN = function (k) {
            simMean <- sapply(1:nIter, FUN = function(i) {
                m <- "Monte Carlo"
                cat(m, ", k=", k, ", rho=", rho, 
                    ", iter ", i, "    \r", sep="")
                Sigma[0, 1] <- Sigma[1, 0] <- rho
                simX <- mvrnorm(n, mu=c(mu1, mu2), Sigma=Sigma)
                return( mean(simX < k) )  
            })
            cat("\nP(X+Y<", k,") = ", mean(simMean), "\n", sep="")
            return(c(var(simMean), mean(simMean)))
        })
        cat("\nVar=", kRes, "\n", sep="")
        return (kRes)
    })
mcRes <- list()
mcRes$mean <- mcRes$var <- 
    matrix(NA, nrow=length(rhos), ncol=length(ks), dimnames=("rho", "k"))
colnames(mcRes$mean) <- colnames(mcRes$var) <- ks
rownames(mcRes$mean) <- rownames(mcRes$var) <- rhos
for (i in 1:length(Res)) {
   mcRes$var[i, ]  <- Res[[i]][1, ]
   mcRes$mean[i, ] <- Res[[i]][2, ]
}
mcRes


# important sampling
# assume p(x) is a folded exponential, lambda=1
#ipRes <- sapply(xs, FUN = function(x) {
#ipRes <- parSapply(c1, xs, FUN = function(x) {
#        simMean <- sapply(1:nIter, FUN = function(i){
#            m <- "Important Sampling"
#            cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
#            simX <- rexp(n) * ((runif(n)>0.5)*2-1)
#            simX <- (x>simX) * dnorm(simX) / (dexp(abs(simX))/2)
#            return(mean(simX)) 
#        })
#        cat("\npnorm(", x, ") = ", mean(simMean), "\n", sep="")
#        return(var(simMean))
#    })
#ipRes




stopCluster(c1)

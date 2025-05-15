############################################################ 
# Filename: q5.R
# Purpose : 113-2 Statistical Computation and Simulation, HW4, Q5
# Author  : Potumas Liu
# Date    : 2025/05/12
############################################################ 

library(MASS)
library(parallel)
library(VGAM, quietly=TRUE)  # bivariate normal evaluation

n <- 1000000
nIter <- 50
rhos <- c(0.7, 0.3, 0, -0.5, -0.9)
ks <- c(0, 1, 2, 3, 4)
mu1 <- mu2 <- 0
Sigma <- matrix(c(1, 0, 0, 1), nrow=2)
rhonames <- c("0.7_var", "0.7_mean", "0.3_var", "0.3_mean",
              "0_var", "0_mean", "-0.5_var", "-0.5_mean",
              "-0.9_var", "-0.9_mean")

RNGkind("L'Ecuyer-CMRG")
seed <- 2025
set.seed(seed)

c1 <- makeCluster(12)
clusterSetRNGStream(c1, seed)
clusterExport(c1, ls())


# original monte carlo
#Res <- sapply(xs, FUN = function(x) {
Res <- parLapply(c1, rhos, fun = function(rho) {
        library(VGAM)
        kRes <- sapply(ks, FUN = function (k) {
            simMean <- sapply(1:nIter, FUN = function(i) {
                #m <- "Monte Carlo"
                #cat(m, ", k=", k, ", rho=", rho, 
                #    ", iter ", i, "    \r", sep="")
                #Sigma[0, 1] <- Sigma[1, 0] <- rho
                #simX <- mvrnorm(n, mu=c(mu1, mu2), Sigma=Sigma)
                #return( mean(simX < k) )  
                simX <- rbinorm(n, cov12=rho)
                xy <- rowSums(simX)

                return( mean(xy < k) )
            })
            cat("\nP(X+Y<", k,") = ", mean(simMean), "\n", sep="")
            return(c(var(simMean), mean(simMean)))
        })
        cat("\nVar=", kRes, "\n", sep="")
        return (kRes)
    })
cat("Monte Carlo\n")
mcRes <- matrix(NA, nrow=length(rhos), ncol=length(ks)*2)
colnames(mcRes) <- colnames(mcRes) <- rhonames
rownames(mcRes) <- rownames(mcRes) <- ks
for (i in 1:length(ks)) {
   mcRes[, i*2-1] <- Res[[i]][1, ]
   mcRes[, i*2]   <- Res[[i]][2, ]
}


# important sampling
#Res <- sapply(xs, FUN = function(x) {
Res <- parLapply(c1, rhos, fun = function(rho) {
        library(VGAM)
        kRes <- sapply(ks, FUN = function (k) {
            simMean <- sapply(1:nIter, FUN = function(i) {
                #m <- "Important Sampling"
                #cat(m, ", k=", k, ", rho=", rho, 
                #    ", iter ", i, "    \r", sep="")

                simX <- rexp(n) * ((runif(n)>0.5)*2-1)
                simY <- rexp(n) * ((runif(n)>0.5)*2-1)
                xy <- simX + simY

                simDat <- (xy<k) * dbinorm(simX, simY, cov12=rho) / 
                        (dexp(abs(simX))*dexp(abs(simY))/4)

                return( mean(simDat) )
            })
            cat("\nP(X+Y<", k,") = ", mean(simMean), "\n", sep="")
            return(c(var(simMean), mean(simMean)))
        })
        cat("\nVar=", kRes, "\n", sep="")
        return (kRes)
    })
cat("Important Sampling\n")
isRes <- matrix(NA, nrow=length(rhos), ncol=length(ks)*2)
colnames(isRes) <- colnames(isRes) <- rhonames
rownames(isRes) <- rownames(isRes) <- ks
for (i in 1:length(ks)) {
   isRes[, i*2-1] <- Res[[i]][1, ]
   isRes[, i*2]   <- Res[[i]][2, ]
}



# Antithetic Variable
#Res <- sapply(xs, FUN = function(x) {
Res <- parLapply(c1, rhos, fun = function(rho) {
        library(VGAM)
        kRes <- sapply(ks, FUN = function (k) {
            simMean <- sapply(1:nIter, FUN = function(i) {
                #m <- "Antithetic Variable"
                #cat(m, ", k=", k, ", rho=", rho, 
                #    ", iter ", i, "    \r", sep="")

                
                simX <- rbinorm(n, cov12=rho)
                xy <- rowSums(simX)

                return( mean((k > xy) + (k > (-xy))) / 2 )
            })
            cat("\nP(X+Y<", k,") = ", mean(simMean), "\n", sep="")
            return(c(var(simMean), mean(simMean)))
        })
        cat("\nVar=", kRes, "\n", sep="")
        return (kRes)
    })
cat("Antithetic Variable\n")
avRes <- matrix(NA, nrow=length(rhos), ncol=length(ks)*2)
colnames(avRes) <- colnames(avRes) <- rhonames
rownames(avRes) <- rownames(avRes) <- ks
for (i in 1:length(ks)) {
   avRes[, i*2-1] <- Res[[i]][1, ]
   avRes[, i*2]   <- Res[[i]][2, ]
}


# Control Variable
#Res <- lapply(rhos, FUN = function(rho) {
Res <- parLapply(c1, rhos, fun = function(rho) {
        library(VGAM)#; library(MASS)
        kRes <- sapply(ks, FUN = function (k) {
            simMean <- sapply(1:nIter, FUN = function(i) {
                #m <- "Antithetic Variable"
                #cat(m, ", k=", k, ", rho=", rho, 
                #    ", iter ", i, "    \r", sep="")
                
                #Sigma[1, 2] <- Sigma[2, 1] <- rho
                #simX <- mvrnorm(n, mu=c(0, 0), Sigma=Sigma)
                simX <- rbinorm(n, cov12=rho)
                xy <- rowSums(simX)

                xy1 <- xy < k
                xy2 <- xy < 0
                #print(mean(xy2 > 0.5))
                return( mean(xy1 - (xy2-0.5)) )
            })
            cat("\nP(X+Y<", k,") = ", mean(simMean), "\n", sep="")
            return(c(var(simMean), mean(simMean)))
        })
        cat("\nVar=", kRes, "\n", sep="")
        return (kRes)
    })
cvRes <- matrix(NA, nrow=length(rhos), ncol=length(ks)*2)
colnames(cvRes) <- colnames(cvRes) <- rhonames
rownames(cvRes) <- rownames(cvRes) <- ks
for (i in 1:length(ks)) {
   cvRes[, i*2-1] <- Res[[i]][1, ]
   cvRes[, i*2]   <- Res[[i]][2, ]
}
cat("Control Variable\n")

allRes <- list(mcRes, isRes, avRes, cvRes)

stopCluster(c1)

# export results
methods <- c("Monte-Carlo", "Important-Sampling", 
             "Antithetic-Variable", "Control-Variable")
for (idx in 1:length(methods)) {
    mRes <- allRes[[idx]]; m <- methods[idx]
    fileName <- paste0("_", m, "_n-", n, "_runs-", nIter, "_seed-", seed)
    write.csv(mRes, paste0("q5-results", fileName, ".csv"))
}

############################################################ 
# Filename: q4.R
# Purpose : 113-2 Statistical Computation and Simulation, HW4, Q4
# Author  : Potumas Liu
# Date    : 2025/05/05
############################################################ 

n <- 10000000
nIter <- 10
xs <- c(-2, -2.5, -3, -3.5, -4, -5, -6) 

set.seed(2025)

# list theoretical F(x)
for (x in xs) cat("pnorm(", x, ") = ", pnorm(x), "\n", sep="")


# original monte carlo
#mcRes <- sapply(xs, FUN = function(x) {
#        simMean <- sapply(1:nIter, FUN = function(i){
#            m <- "Monte Carlo"
#            cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
#           return(mean(rnorm(n) < x)) 
#        })
#        cat("\npnorm(", x, ") = ", mean(simMean), "\n", sep="")
#        return(var(simMean))
#    })
#mcRes


# important sampling
# assume p(x) is a folded exponential, lambda=1
#ipRes <- sapply(xs, FUN = function(x) {
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


# antithetic variable
#avRes <- sapply(xs, FUN = function(x) {
#        simMean <- sapply(1:nIter, FUN = function(i){
#            m <- "Antithetic Variable"
#            cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
#            simX1 <- rnorm(n)
#            simX2 <- -simX1
#            return( mean((x > simX1) + (x > simX2)) / 2 )
#        })
#        cat("\npnorm(", x, ") = ", mean(simMean), "\n", sep="")
#        return(var(simMean))
#    })
#avRes


# control variate
cvRes <- sapply(xs, FUN = function(x) {
        simMean <- sapply(1:nIter, FUN = function(i){
            m <- "Control Variate"
            cat(m, ", x = ", x, ", iter ", i, "    \r", sep="")
            simX1 <- rnorm(n)
            #simX2 <- rexp(n)
            simX2 <- rexp(n) * ((runif(n)>0.5)*2-1)
            return( mean( (simX1<x) - simX2 ) )
        })
        cat("\npnorm(", x, ") = ", mean(simMean), "\n", sep="")
        return(var(simMean))
    })
cvRes



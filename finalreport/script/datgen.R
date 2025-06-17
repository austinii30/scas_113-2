###########################################################
# Filename: datgen.R
# Purpose : Generate spatial data for the final report 
# Author  : Austin Liu
# Date    : 2025/06/14
###########################################################

source(file.path(this.path::this.dir(), "utils.R"))

set.seed(2025)

#N <- 700      # size for each sample
N <- 500
nsub <- 1000 
gridsize <- 10  # size of study area
errsig <- 0.5  # sigma for the random error
#errsig <- 0.1  # sigma for the random error
#errsig <- 1  # sigma for the random error
#errsig <- 1.5  # sigma for the random error
#errsig <- 2  # sigma for the random error
#errsig <- 2.5  # sigma for the random error
#errsig <- 3  # sigma for the random error
varnames <- c("response", "x", "y", "x1", "x2", "x3", "error")
lmparas <- c(0.3, 0.7, -0.8, 0.5)  # betas for LM
gwrpar0 <- function(x, y) { return( (x-y)/10 ) }
gwrpar1 <- function(x, y) { return( (sin(x/0.5)+cos(y/0.8))/2 ) }
gwrpar2 <- function(x, y) { return( exp(-((x-2.5)^2+(y-3)^2)/2) + exp(-((x-8)^2+(y-6)^2)/4) - exp(-((x-7)^2+(y-2)^2)/0.8) - exp(-((x-4)^2+(y-8)^2)/1.5)) }
gwrpar3 <- function(x, y) { return( exp(-(y-2*(x-1))^2/3) + exp(-(3*(y-3)-x)^2/10) - 1) }


# --------------------------------------------
# Linear Global Regression
# --------------------------------------------
datlm <- datgwr <- lapply(1:nsub, function(i) {return(NULL)} )
for (i in 1:nsub) {
    x <- runif(N)*gridsize; y <- runif(N)*gridsize  # location
    # covariates
    x0 <- rep(1, N)  
    x1 <- runif(N)
    x2 <- runif(N)
    x3 <- runif(N)
    covariates <- cbind(x0, x1, x2, x3)
    # error
    error <- rnorm(N, mean=0, sd=errsig)
    # response 
    response.lm <- covariates %*% lmparas + error
     
    # store the generated data
    subdat.lm <- cbind(response.lm, x, y, x1, x2, x3, error)
    colnames(subdat.lm) <- varnames
    datlm[[i]] <- subdat.lm

    # gwr parameters
    gp0 <- sapply(1:N, function(idx) { return(gwrpar0(x[idx], y[idx])) } )
    gp1 <- sapply(1:N, function(idx) { return(gwrpar1(x[idx], y[idx])) } )
    gp2 <- sapply(1:N, function(idx) { return(gwrpar2(x[idx], y[idx])) } )
    gp3 <- sapply(1:N, function(idx) { return(gwrpar3(x[idx], y[idx])) } )
    gwrparas <- cbind(gp0, gp1, gp2, gp3)
    # response 
    response.gwr <- rowSums(covariates * gwrparas) + error
     
    # store the generated data
    subdat.gwr <- cbind(response.gwr, x, y, x1, x2, x3, error, 
                    gp0, gp1, gp2, gp3)
    colnames(subdat.gwr) <- c(varnames, "gp0", "gp1", "gp2", "gp3")
    datgwr[[i]] <- subdat.gwr
}
save(datlm, file=datpath(paste0("datlm-sigma-", errsig, ".RData")))
save(datgwr, file=datpath(paste0("datgwr-sigma-", errsig, ".RData")))

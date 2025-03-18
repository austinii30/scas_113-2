
source("code.R")


# 1.(a)
set.seed(2025)

nIter <- 1000
n <- 100
theta <- c(0, 0.05, 0.1, 0.15, 0.2)

simDat <- list()
# for each theta, draw 1000 random samples, each with n=100
for (th in theta) {
    dat <- c()
    for (i in 1:nIter)
        dat <- c(dat, arima.sim(n=n, list(ar=c(th, th))))
    # a column is a sample (1000 columns in total)
    simDat[[length(simDat)+1]] <- matrix(dat, ncol=nIter)
}


# calculate correlation
simDatCor <- list()
for (dat in simDat) {
    # returns a matirx with 2 rows, each is a correlation
    simCor <- sapply(1:ncol(dat), FUN = function(sIdx) {
            sample <- dat[, sIdx]
            cor1 <- cor(sample[-1], sample[-n])
            cor2 <- cor(sample[-c(1, 2)], sample[-c(n-1,n)])
            return(c(cor1, cor2))
        }) 
    # transpose so row: a sample, col: cor1 & cor2
    simDatCor[[length(simDatCor)+1]] <- t(simCor)
}














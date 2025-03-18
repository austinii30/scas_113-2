
# 1.(a)
set.seed(2025)

nIter <- 1000
n <- 100
theta <- c(0, 0.05, 0.1, 0.15, 0.2)

simDat <- list()
for (th in theta) {
    dat <- c()
    for (i in 1:nIter)
        dat <- c(dat, arima.sim(n=n, list(ar=c(th, th))))
    simDat[[length(simDat)+1]] <- matrix(dat, nrow=nIter, byrow=T)
}

# calculate correlation
simCor <- sapply(simDat, FUN = function(x) {
        
    })





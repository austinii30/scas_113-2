############################################################ 
# Filename: code.R
# Purpose : 113-2 Statistical Computation and Simulation
# Author  : Potumas Liu
# Date    : 2025/03/06
############################################################ 



# 1.(a)
# ================================================================
rm_1a <- function(n, seed=123, rmSeed=NULL) {
    # endogenous initial value
    # ex: seed=123, x=1.23; seed=34, x=3.4
    k = 0
    while (seed %/% (10^k) != 0)
       k = k+1 
    x <- seed / (10^(k-1))
    
    # initial value by runif()
    if (!is.null(rmSeed)) {
        set.seed(rmSeed)
        x <- runif(1)
    }

    res <- c()
    for (i in 1:n) {
        x <- exp(x)
        x <- x - floor(x)
        res <- c(res, x)
    }

    return(res)
}




# 1.(b)
# ================================================================
rm_1b <- function(n, seed=123, rmSeed=NULL) {
    modX <- 30269; modY <- 30307; modZ <- 30323
    mltX <- 171  ; mltY <- 172  ; mltZ <- 170

    # endogenous initial value
    # ex: seed=123, x=1.23; seed=34, x=3.4
    x <- seed/mltX+modX
    y <- seed/mltY+modY  
    z <- seed/mltZ+modZ

    # initial value by runif()
    if (!is.null(rmSeed)) {
        set.seed(rmSeed)
        x <- runif(1); y <- runif(1); z <- runif(1)
    }

    res <- c()
    for (i in 1:n) {
        x <- (mltX*x) %% modX
        y <- (mltY*y) %% modY
        z <- (mltZ*z) %% modZ

        u <- (x/modX + y/modY + z/modZ) %% 1

        res <- c(res, u)
    }

    return(res)
}

# ================================================================



# 2.(a)
# ================================================================ 
fibN <- function(n, m, sameInit=F, seed=123, rmSeed=NULL) {
    if (n < m+1) stop("'n' must at least be 'm+1'!")

    # m+1 initial values must be given first
    # endogenous initial value
    if (sameInit) 
        res <- rep(rm_1b(1, seed=seed), m+1)
    else
        res <- rm_1b(m+1, seed=seed)
    
    # initial value by runif()
    if (!is.null(rmSeed)) {
        set.seed(rmSeed)
        if (sameInit)
            res <- rep(runif(1), m+1)
        else
            res <- runif(m+1)
    }

    i <- m + 2 
    while (i <= (n+m+1)) {
        x <- (res[i-1] + res[i-1-m]) %% 1
        res <- c(res, x)
        i <- i + 1
    }

    return(res[-(1:(m+1))])
}



# 2.(b)
# ================================================================ 
rm_2b <- function(n, k, method, rmSeed){
    if (!method %in% c("sample", "ceiling"))
        stop("'method' should only be 'sample' or 'ceiling'")

    set.seed(rmSeed)

    if (method == "sample")
        res <- sample(1:k, n, replace=T)     

    if (method == "ceiling")
        res <- ceiling(runif(n, min=0, max=k))

    return(res)
}









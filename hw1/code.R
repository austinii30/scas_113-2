
# 1.(a)
# ================================================================
rm <- function(n, seed=123) {
    # initial value ex: seed=123, x=1.23; seed=34, x=3.4
    k = 0
    while (seed %/% (10^k) != 0)
       k = k+1 
    x <- seed / (10^(k-1))
    
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
rm_1b <- function(n, seed=123) {
    modX <- 30269; modY <- 30307; modZ <- 30323
    mltX <- 171  ; mltY <- 172  ; mltZ <- 170

    x <- seed/mltX+modX
    y <- seed/mltY+modY  
    z <- seed/mltZ+modZ

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

rm(10)
rm_1b(10)




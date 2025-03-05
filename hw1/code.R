
# 1.(a)
# ================================================================
rm <- function(n, seed=123) {
    x <-  
    
    res <- c()
    for (i in 1:n) {
        x <- exp(x)
        x <- x - floor(x)
        res <- c(res, x)
    }

    return res
}




# 1.(b)
# ================================================================
rm_1b <- function(n, x0, y0, z0) {
    x <- x0; y <- y0; z <- z0

    modX <- 30269; modY <- 30307; modZ <- 30323
    mltX <- 171  ; mltY <- 172  ; mltZ <- 170

    res <- c()
    for (i in 1:n) {
        x <- (mltX*x) %% modX
        y <- (mltY*y) %% modY
        z <- (mltZ*z) %% modZ

        u <- (x/modX + y/modY + z/modZ) %% 1

        res <- c(res, u)
    }

    return res
}


# ================================================================




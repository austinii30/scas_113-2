

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







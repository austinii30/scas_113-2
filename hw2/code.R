









# 3.(a)
rN.1 <- function (n, seed=123) {
    set.seed(seed)

    Nvec <- c()
    for (i. in 1:n) {
        x <- runif(1)
        N <- 1 
        while (x <= 1) {
            x <- x + runif(1)
            N <- N + 1
        }
        Nvec <- c(Nvec, N)
    }

    return(Nvec)
}

rN.2 <- function (n, seed=123) {
    set.seed(seed)

    Nvec <- sapply(1:n, FUN = function(i.) {
            x <- runif(1)
            N <- 1 
            while (x <= 1) {
                x <- x + runif(1)
                N <- N + 1
            }
            return(N)
        })

    return(Nvec)
}







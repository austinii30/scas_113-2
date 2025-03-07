############################################################ 
# Filename: main.R
# Purpose : 113-2 Statistical Computation and Simulation
# Author  : Potumas Liu
# Date    : 2025/03/06
############################################################ 

library(randtoolbox)  # for Gap test
source("./code.R")


# 1 (a)
# ================================================================
x.1a <- rm_1a(10000, rmSeed=2025)
print(ks.test(x.1a, "punif"))


# 1 (b)
# ================================================================
x.1b <- rm_1b(10000, rmSeed=2025)
print(ks.test(x.1b, "punif"))


stop()

# 2 (a)
# ================================================================
x.2a <- list()
for (i in 1:20)
    x.2a[[length(x.2a)+1]] <- fibN(10000, m=i, rmSeed=2025)

testResult <- lapply(x.2a, FUN = function(x) {
        # uniform
        uni.p <- ks.test(x, "punif")$p
        # independent tests
#       gap.p <- (x)
        #per.p <- ()
        #uad.p <- ()
        
#        res <- c(uni.p, gap.p, per.p, uad.p)
#        return(res)
    })



# 2 (b)
# ================================================================
n = 100000
x.s <- x.c <- c()

K <- 5:50
K <- c(50, 100, 200, 300, 400, 500)
K <- 300:500
for (k in K) {
    x.s <- c(x.s, rm_2b(n, k, "sample", rmSeed=2025))
    x.c <- c(x.c, rm_2b(n, k, "ceiling", rmSeed=2025))
}

# each column is a random number sample of size n
x.s <- matrix(x.s, nrow = n)
x.c <- matrix(x.c, nrow = n)
x.all <- cbind(x.s, x.c)

pValues <- c()
for (i in 1:ncol(x.all)) {
    p <- chisq.test(table(x.all[, i]))$p.value
    pValues <- c(pValues, p)
}
pValues <- matrix(pValues, ncol=2, 
                  dimnames = list(k=K, 
                                  method=c("sample", "ceiling")))
print(pValues)



# 3 
# ================================================================
stop()
N <- 1000
n10 <- 10; n50 <- 50; n100 <- 100

res <- c()
#for () {
#    res <- c(res, rnorm((n10+n50+n100)*N))

#}

t10.10  <- rt(10, df = 10)
t50.10  <- rt(50, df = 10)
t100.10 <- rt(100, df = 10)
t10.20  <- rt(10, df = 20)
t50.20  <- rt(50, df = 20)
t100.20 <- rt(100, df = 20)

















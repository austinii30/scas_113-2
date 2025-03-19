
source("code.R")
set.seed(2025)
nIter <- 10

# 3.(a)
n <- c(1000, 2000, 5000, 10000, 100000)
# each element is a 'n'
NSamples <- lapply(n, FUN = function(x) {
    # each row is a sample
    return(
        t(sapply(1:nIter, FUN = function(i.) return(rN.2(x))))
    )
})

Nmean <- lapply(NSamples, FUN = function(mtx) return(rowMeans(mtx)))
hist(Nmean)

Nsd <- lapply(Nmean, FUN = sd)






# 3.(b)



stop()

# test execution time (rN.2 is much faster)
beg <- Sys.time()
rN.1(100000)
end <- Sys.time()
cat("rN.1, time: ", beg-end, "\n")

beg <- Sys.time()
rN.2(100000)
end <- Sys.time()
cat("rN.2, time: ", beg-end, "\n")


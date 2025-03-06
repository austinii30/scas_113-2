#
source("./code.R")

# 1 (a)
# ================================================================
x.1a <- rm_1a(10000, rmSeed=2025)
ks.test(x.1a, "punif")

# 1 (b)
# ================================================================
x.1b <- rm_1b(10000, rmSeed=2025)
ks.test(x.1b, "punif")



# 2 (a)
# ================================================================
x.2a <- list()
for (i in 1:20)
    x.2a[[length(x.2a)+1]] <- fibN(10000, m=i, rmSeed=2025)

testResult <- lapply(x.2a, FUN = function(x) {
        # uniform
        uni.p <- ks.test(x, "punif")$p
        # independent tests
        gap.p <- ()
        per.p <- ()
        uad.p <- ()
        
        res <- c(uni.p, gap.p, per.p, uad.p)
        return(res)
    })


# 2 (b)
# ================================================================




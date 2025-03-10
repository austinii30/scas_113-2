############################################################ 
# Filename: code.R
# Purpose : 113-2 Statistical Computation and Simulation
# Author  : Potumas Liu
# Date    : 2025/03/06
############################################################ 


library(nortest)



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
fibN <- function(n, m, rmSeed=NULL) {
    if (n < m+1) stop("'n' must at least be 'm+1'!")

    set.seed(rmSeed)
    res <- runif(m+1)

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



# Gap Test
# ================================================================ 
gapTest <- function(a, b, data) {
    alpha <- a; beta  <- b

    binary_data <- ifelse(data >= alpha & data <= beta, 1, 0)
    gaps <- c()  

    count <- 0  
    for (i in 1:length(binary_data)) {
        if (binary_data[i] == 1) {
            gaps <- c(gaps, count)  
            count <- 0               
        } else {
            count <- count + 1       
        }
    }

    #如果最後結束時還有未儲存的 gap，則加入 (避免末尾未處理)
    #if (count > 0) {
    #    gaps <- c(gaps, count)
    #}

    gap_table  <- table(gaps)
    gap_freq   <- as.numeric(gap_table)
    gap_values <- as.numeric(names(gap_table))
    p <- beta - alpha  # success probability
    expected <- length(gaps) * (1 - p)^(gap_values) * p  # theoretical expectations for each group

    threshold <- 5  
    cb_gap_freq <- gap_freq  # realized frequency
    cb_expected <- expected  # theoretical frequency

    # ensure all groups have NP>5
    # if more than one group has NP<5, keep looping
    # for the first group with NP<5, combine it with the next group, until its NP>5
    # ex: group 5 has NP<5 --> combine group 5 with group 6, 7, ... until NP>5
    # (if only the last group remains NP<5, it will be dealt with later)
    while(sum(cb_expected < threshold) > 1) {
        first_small_idx <- which(expected < threshold)[1]
        
        # keep binding until NP>5
        smallExp <- cb_expected[first_small_idx]  # the expectation with NP<5
        i <- 0  # how many subsequent groups have been added
        while(smallExp < threshold) {
            i <- i + 1
            if (length(cb_expected) < (first_small_idx+i)) {
                i <- i - 1  # make sure i don't exceed the last group
                break
            }

            smallExp <- smallExp + cb_expected[first_small_idx+i]
        }

        cb_expected[first_small_idx] <- smallExp
        delete_idx <- (1:i) + first_small_idx 
        cb_expected <- cb_expected[-delete_idx]

        cb_gap_freq[first_small_idx] <- sum(c(cb_gap_freq[first_small_idx], cb_gap_freq[delete_idx]))
        cb_gap_freq <- cb_gap_freq[-delete_idx]
    }

    final <- length(cb_expected)
    if (cb_expected[final] < threshold) {
        cb_expected[final-1] <- sum(cb_expected[(final-1):final])
        cb_expected <- cb_expected[-final]

        cb_gap_freq[final-1] <- sum(cb_gap_freq[(final-1):final])
        cb_gap_freq <- cb_gap_freq[-final]
    }

#        # 合併到最後一格
#        cb_gap_freq[length(gap_freq)] <- sum(gap_freq[small_idx])
#        cb_expected[length(expected)] <- sum(expected[small_idx])
#        # 刪除過小的部分 (除了最後一格)
#        cb_gap_freq <- cb_gap_freq[-small_idx[-length(small_idx)]]
#        cb_expected <- cb_expected[-small_idx[-length(small_idx)]]

    chi_test <- chisq.test(cb_gap_freq, p = cb_expected / sum(cb_expected))

    return(chi_test)
}



uadTest <- function(data) {
    N <- length(data)

    ups_and_downs <- diff(data) > 0
    runs <- sum(diff(ups_and_downs) != 0)
    R <- runs+1
    E_R <- (2*N - 1) / 3
    Var_R <- (16*N - 29) / 90
    Z <- (R - E_R) / sqrt(Var_R)

    cd <- pnorm(Z)  # pnorm(): cdf
    if (cd > 0.5) return((1-cd)*2)
    else          return(cd*2)
}



# 3
# ================================================================ 
func3 <- function(n, nIter=1000, seed=2025, alpha=0.05) {
    set.seed(seed)
    res <- matrix(0, nrow=4, ncol=3)

    for (i in 1:nIter) { 
        cat("n = ", n, ", iter = ", i, "\r", sep="")
        # generate random samples
        nSample   <- rnorm(n)
        t10Sample <- rt(n, df = 10)
        t20Sample <- rt(n, df = 20)
        samples <- list(nSample, t10Sample, t20Sample)

        # conduct normality tests
        # store the results in matrix
        testRes <- sapply(samples, FUN = function(x) {
                ksP <- ks.test(x, "pnorm")$p
                adP <- ad.test(x)$p
                swP <- shapiro.test(x)$p
                llP <- lillie.test(x)$p
                return(c(ksP, adP, swP, llP))
            })

        testRes <- testRes < alpha  # reject or not
        res <- res + testRes
    }
    cat("\n")

    rownames(res) <- c("KS", "AD", "SW", "LL")
    colnames(res) <- c("N(0,1)", "t(10)", "t(20)")
    print(res)

    return(res)
}






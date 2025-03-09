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

    gap_table <- table(gaps)
    gap_freq <- as.numeric(gap_table)
    gap_values <- as.numeric(names(gap_table))
    p <- beta - alpha
    expected <- length(gaps) * (1 - p)^(gap_values) * p

    threshold <- 5  
    combined_gap_freq <- gap_freq
    combined_expected <- expected

    #cat("\n\n")
    #print(combined_expected)
    # 合併
    while(sum(combined_expected < threshold) > 1) {
    #if (any(expected < threshold)) {
        # 找出過小的部分
        first_small_idx <- which(expected < threshold)[1]
        
        smallExp <- combined_gap_freq[first_small_idx]
        i <- 0
        while(smallExp < threshold) {
            i <- i + 1
            if (length(combined_gap_freq) < (first_small_idx+i)) {
                i <- i - 1
                break
            }

            smallExp <- smallExp + combined_gap_freq[first_small_idx+i]
        }

        combined_gap_freq[first_small_idx] <- smallExp
        deleted_idx <- (1:i) + first_small_idx 
        combined_gap_freq <- combined_gap_freq[-deleted_idx]

        combined_expected[first_small_idx] <- sum(c(combined_expected[first_small_idx], combined_expected[deleted_idx]))
        combined_expected <- combined_expected[-deleted_idx]
    }

    final <- length(combined_gap_freq)
    if (combined_expected[final] < threshold) {
        combined_gap_freq[final-1] <- sum(combined_gap_freq[(final-1):final])
        combined_gap_freq <- combined_gap_freq[-final]

        combined_expected[final-1] <- sum(combined_expected[(final-1):final])
        combined_expected <- combined_expected[-final]
    }
    #print(combined_expected)

#        # 合併到最後一格
#        combined_gap_freq[length(gap_freq)] <- sum(gap_freq[small_idx])
#        combined_expected[length(expected)] <- sum(expected[small_idx])
#        # 刪除過小的部分 (除了最後一格)
#        combined_gap_freq <- combined_gap_freq[-small_idx[-length(small_idx)]]
#        combined_expected <- combined_expected[-small_idx[-length(small_idx)]]

    chi_test <- chisq.test(combined_gap_freq, p = combined_expected / sum(combined_expected))

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


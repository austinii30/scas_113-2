############################################################ 
# Filename: code.R
# Purpose : 113-2 Statistical Computation and Simulation, HW1
#           Some handy functions
# Author  : Potumas Liu
# Date    : 2025/03/06
############################################################ 


library(nortest)


# 1.(a)
# ================================================================
rm_1a <- function(n, seed=2025) {
    set.seed(seed)

    x <- runif(1)
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
rm_1b <- function(n, seed=2025) {
    modX <- 30269; modY <- 30307; modZ <- 30323
    mltX <- 171  ; mltY <- 172  ; mltZ <- 170

    # initial value 
    set.seed(seed)
    x <- runif(1); y <- runif(1); z <- runif(1)

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



# 2.(a)
# ================================================================ 
fibN <- function(n, m, seed=2025) {
    if (n < m+1) stop("'n' must at least be 'm+1'!")

    set.seed(seed)
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
rm_2b <- function(n, k, method, seed=2025) {
    if (!method %in% c("sample", "ceiling"))
        stop("'method' should only be 'sample' or 'ceiling'")

    set.seed(seed)

    if (method == "sample")
        res <- sample(1:k, n, replace=T)     

    if (method == "ceiling")
        res <- ceiling(runif(n, min=0, max=k))

    return(res)
}



# Gap Test (returns 'chisq.test()' object)
# ================================================================ 
gapTest <- function(a, b, data) {
    alpha <- a; beta <- b
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
    p <- beta - alpha
    expected <- length(gaps) * (1 - p)^(gap_values) * p
    
    threshold <- 5  
    combined_gap_freq <- gap_freq
    combined_expected <- expected
    
    if (any(expected < threshold)) {
        # 找出過小的部分索引
        small_idx <- which(expected < threshold)
        
        # 從尾端開始逐步合併直到所有的 expected >= threshold
        for (i in rev(small_idx)) {
            if (expected[i] < threshold) {
                # 合併當前格到前一格
                if (i > 1) {  # 確保不超出範圍
                    gap_freq[i - 1] <- gap_freq[i - 1] + gap_freq[i]
                    expected[i - 1] <- expected[i - 1] + expected[i]
                }
            }
        }
    
        # 刪除已合併的部分
        combined_gap_freq <- gap_freq[expected >= threshold]
        combined_expected <- expected[expected >= threshold]
    } else {
        combined_gap_freq <- gap_freq
        combined_expected <- expected
    }
    
    chi_test <- chisq.test(combined_gap_freq, p = combined_expected / sum(combined_expected))

    return(chi_test)
}



# Up and Down Test (returns the p-value)
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

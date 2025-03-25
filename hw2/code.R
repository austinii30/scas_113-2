








# 3.(a)
# Don't use this funciton!
# 'for' is much more slower than 'apply'
rN.1 <- function (n) {
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

rN.2 <- function (n) {
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





# Gap Test
gapTest <- function(data, a=0.25, b=0.7) {
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

    #print("-----------------")
    #print(combined_gap_freq)
    #print(combined_expected)
    
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
    
    #print(combined_gap_freq)
    if (length(combined_gap_freq) <= 1)
        return(NA)
        #stop("Must have at least 2 groups to perform a chi-square test (May arise from small sample size).")

    chi_test <- chisq.test(combined_gap_freq, p = combined_expected / sum(combined_expected))

    return(chi_test$p.value)
}



# Up and Down Test
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



# Permutation Test
pmtTest <- function (dat, k=3) {
    # remove the last remaining m numbers (m < k)
    dat <- dat[1:(length(dat) - (length(dat) %% k))]
    # each row is a group
    data_matrix <- matrix(dat, ncol=k, byrow=TRUE)
    rank_matrix <- t(apply(data_matrix, 1, function(x) order(x)))
    # change the order (ex: 1,2,3 -> 123) to calculate frequency
    rank_strings <- apply(rank_matrix, 1, paste, collapse = "")
    rank_freq <- table(rank_strings)

    # add 0 to make sure length of rank_freq is same as factorial(k)
    rank_freq <- c(rank_freq, rep(0, factorial(k)-length(rank_freq)))
    chi_test <- chisq.test(rank_freq, p=rep(1/factorial(k), factorial(k)))
    return (chi_test$p.value)
}



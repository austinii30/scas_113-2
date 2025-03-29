###########################################################
# Filename: code.R
# Course  : 113-2 Statistical Computation and Simulation, HW2,
#           Some usable functions.
# Author  : Potumas Liu 
# Date    : 2025/03/24
###########################################################

source("../../myPKG/R/EDA.R")  # for HDB()



HDB <- function(dat, varName="variable", ct="Set1", limits=NULL) {
    "
    Histogram, density, and box plots for numeric and integer variables.
    [Args]
        dat (num)     : the data vector
        varName (char): the name of the variable
        ct (char)     : the color theme
        limits (num)  : limit of the x-axis (c(min, max))
    [Return] 
        A ggplot object ('print(_obj_)' to plot it)
    "
    #devtools::install_github("psyteachr/introdataviz")
    #install.packages("Hmisc")
    cat("Drawing univariate raincloud plots...\n")

    if (! class(dat) %in% c("integer", "numeric"))
        stop("'dat' can only be 'integer' or 'numeric'!")

    df <- data.frame(dat)
    colnames(df) <- c(varName)
    
    colortheme <- ct     # 'Set1', 'Dark2'
    legendPos <- "none"
    titleText <- 16
    axisText  <- 21
    rain_height <- 0.1
    
    # add a variable 'dummy' which all rows have equal value
    df$dummy <- factor("same")

    # draw each plot separately
    for (varName in colnames(df)) {
        if (varName == "dummy") next

        dynamicCode <- "
var <- df$`__var__`
varName <- \"__var__\"

fig <- 
ggplot(df, aes(x = \"\", y = __var__, fill = dummy)) +

# clouds
# 'alpha' controls the transparancy of the graph (1: opague, 0: transparent)
introdataviz::geom_flat_violin(trim=TRUE, alpha = 0.4,
                               position = position_nudge(x = rain_height+.05)) +
# rain
geom_point(aes(colour = dummy), size = 2, alpha = .5, show.legend = FALSE, 
           position = position_jitter(width = rain_height, height = 0)) +
# boxplots
geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
             outlier.shape = NA,
             position = position_nudge(x = -rain_height*2)) +
# mean and SE point in the cloud
#stat_summary(fun.data = mean_cl_normal, mapping = aes(color = dummy), show.legend = FALSE,
stat_summary(fun = \"mean\", mapping = aes(color = dummy), show.legend = FALSE,
             position = position_nudge(x = rain_height * 3)) +
# adjust layout
#ggtitle(varName) + 
    scale_x_discrete(name = \"\", expand = c(rain_height*3, 0, 0, 0.7)) +
    scale_y_continuous(name = \"\", limits = limits) + 
    coord_flip() +
    # custom colours and theme
    scale_fill_brewer(palette = colortheme, name = \"Location\") +
    scale_colour_brewer(palette = colortheme) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = legendPos,
          legend.background = element_rect(fill = \"white\", color = \"white\"),
          text = element_text(size = 14), 
          plot.title = element_text(size = titleText,  hjust = 0.5),
          axis.text = element_text(size = axisText))
"
        dynamicCode <- gsub("__var__", varName, dynamicCode)
        return(eval(parse(text=dynamicCode)))
    }
}



# Don't use this funciton!
# 'for' is much more slower than 'apply'
# Generate N for N=min sum(U)>1 for U ~ unif(0,1)
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



rnorm_RejExp<- function(n, lambda=1) {
    " 
    Sample normal samples by rejection method from an exponential distribution.
    [Args]
        n (num)     : amount of samples
        lambda (num): parameter for the exponential distribution

    [Return]
        numeric: the numerical vector
    "
    normSample <- c()

    C <- lambda/sqrt(2*pi) * exp(1 / (2 * lambda^2))
    
    i <- 0
    k <- 0
    while (i < n) {
        u1 <- runif(1); u2 <- runif(1); u3 <- runif(1)
        k <- k + 1

        x <- (-1) * lambda* log(u1)
        if (u2 <= (dnorm(x) / (dexp(x, rate=(1/lambda)) * C)))
            i <- i + 1
        else
            next
        if (u3 < 0.5) x <- x * (-1)

        normSample <- c(normSample, x)
    }

    res <- list(sample=normSample, nIter=k, acceptRate=(n/k), C=C)
    return(res)
}

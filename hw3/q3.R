######################
# Question 2 for hw3
######################

set.seed(113354)

alpha <- 0.05
lb <- 2; ub <- 10
nIter <- 10000

cv <- matrix(nrow=0, ncol=length(lb:ub))

for (n1 in lb:ub) {
    cvRow <- c()
    for (n2 in lb:ub) {
        allU <- c()
        cat("N1: ", n1, ", N2: ", n2, "\n", sep="")
        for (i in 1:nIter) {
            # draw samples 
            x1 <- rexp(n1)
            #x1 <- runif(n1)
            x2 <- rexp(n2)
            #x2 <- runif(n2)
        
            # rank the observations
            rk <- rank(c(x1, x2))
            # calculate U
            u1 <- sum(rk[1:n1])    - (n1*(n1+1)/2)
            u2 <- sum(rk[n1+1:n2]) - (n2*(n2+1)/2)
            u <- min(u1, u2)
            allU <- c(allU, u)
        }
        # get the critical for alpha
        allU <- sort(allU)
        pctg <- order(allU) / length(allU)
        thisCV <- allU[which(pctg >= (alpha/2))[1]]
        
        cvRow <- c(cvRow, thisCV)
    }
    cv <- rbind(cv, cvRow)
}

colnames(cv) <- rownames(cv) <- lb:ub
cv

write.csv(cv, paste0("./q3_table_a", alpha, ".csv"))

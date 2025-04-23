###########################################################
# Filename: main_4.R
# Course  : 113-2 Statistical Computation and Simulation, HW2,
#           Question 4.
# Author  : Potumas Liu 
# Date    : 2025/03/24
###########################################################

source("code.R")

set.seed(2025)
n <- 10000
lambda <- 1:100

# draw normal samples (can be commented later)
res <- lapply(lambda, function(x) rnorm_RejExp(n=n, lambda=x))
save(res, file="4-NormalSampleResults.RData")


load("4-NormalSampleResults.RData")


# examine the pattern
# acceptance rate for different lambda
acpt <- sapply(res, FUN = function(x) return(x$acceptRate))

pdf("4-acceptRate.pdf", width=7, height=5)

par(mar = c(4, 4, 1, 1))
par(mgp = c(3, 1, 0))
plot(x=lambda, y=acpt, 
     ylab="Acceptance Rate", xlab="Lambda",
     las=1, pch=21, bg="skyblue", type="o", bty="l",
     cex=1.2 , cex.lab=1.3, cex.axis=1.3)

dev.off()


# C for each lambda
Cs <- sapply(res, FUN = function(x) return(x$C))
pdf("4-C.pdf", width=7, height=5)

par(mar = c(4, 4, 1, 1))
par(mgp = c(3, 1, 0))
plot(x=lambda, y=Cs, 
     ylab="C", xlab="Lambda",
     las=1, pch=21, bg="orange", type="o", bty="l",
     cex=1.2 , cex.lab=1.3, cex.axis=1.3)

dev.off()


# theoretical and sampled distribution of the samples
# plot the smoothed density
pickLambda <- c(1, 100)
pRes <- matrix(ncol=2, nrow=0)
colnames(pRes) <- c("gofTest.p", "indTest.p")

for (l in pickLambda) {
    # plot density of each sample
    pdf(paste0("4-normal_lambda-", l, ".pdf"), width=7, height=5)
    
    par(mar = c(4, 4, 1, 1))
    dens <- density(res[[l]]$sample)
    plot(dens, las=1, bty="l", 
         xlab = "X", ylab = "Density", main="",
         cex=1.8, cex.lab=1.3, cex.axis=1.3,
         col = "blue", lwd = 1.5, lty=2)
    curve(dnorm(x), from=-5, to=5, add=T, col="red")
    legend("topright", legend=c("Theoretical", "Simulation"), col=c("red", "blue"), lty=c(1, 2), lwd=1.5, bty="n", cex=1.3)
    text(3, 0.25, paste0("Acceptance Rate:\n", round(res[[l]]$acceptRate, 4)), cex = 1.3)
    
    dev.off()
    
    # goodness-of-fit test
    gof.p <- ks.test(res[[l]]$sample, "pnorm")$p.value
    # independent test (Permutation test)
    ind.p <- pmtTest(res[[l]]$sample)
    
    pRes <- rbind(pRes, c(gof.p, ind.p))
}

write.csv(pRes, "4-testResults.csv")

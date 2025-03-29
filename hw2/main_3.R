###########################################################
# Filename: main_3.R
# Course  : 113-2 Statistical Computation and Simulation, HW2,
#           Question 3(a).
# Author  : Potumas Liu 
# Date    : 2025/03/24
###########################################################

source("code.R")

set.seed(2025)
k <- 100  # how many samples for each N

n <- c(1000, 2000, 5000, 10000, 100000)
# each element is a 'n'
NSamples <- lapply(n, FUN = function(x) {
    # each row is a sample
    return(
        t(sapply(1:k, FUN = function(i.) return(rN.2(x))))
    )
})

# sample mean
Nmean <- lapply(NSamples, FUN = function(mtx) return(rowMeans(mtx)))

# sample standard error
Nsd <- sapply(Nmean, FUN = sd)


# calculate and plot the mean of Nbar
Nmeanavg <- sapply(Nmean, FUN = mean)
pdf(paste0("3.a_mean(Nbar)_n", n, "_k", k, ".pdf"), width=7, height=5)

par(mar = c(4, 6.3, 1, 1))
par(mgp = c(4, 1, 0))
plot(y=Nmeanavg, x=1:5,
     ylab="", xlab="", main="",
     pch=21, col="black", bg="orange", type="o",
     cex=1.8, cex.lab=1.3, cex.axis=1.3,
     yaxt="n", xaxt="n", las=1, bty="l")
# custom x-axis labels and axis
axis(1, at=1:5, labels=c("1,000", "2,000", "5,000", "10,000", "100,000"), cex.axis=1.3)  # Custom x-axis labels
mtext("n", side=1, line=3, cex=1.3) 
# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.3, cex.lab=1.3)
mtext(expression(paste("Mean of ", bar("N"))), side=2, line=4.8, cex=1.3) 

dev.off()


# plot sd(N) ~ n
pdf(paste0("3.a_sd(N)_n", n, "_k", k, ".pdf"), width=7, height=5)

par(mar = c(4, 5.2, 1, 1))
par(mgp = c(4, 1, 0))
plot(y=Nsd, x=1:5, 
     ylab="Standard Error of N", xlab="", main="",
     pch=21, col="black", bg="skyblue", 
     cex=1.8, cex.lab=1.3, cex.axis=1.3, type="o",
     xaxt="n", las=1, bty="l")
axis(1, at=1:5, labels=c("1,000", "2,000", "5,000", "10,000", "100,000"), cex.axis=1.3)  # Custom x-axis labels
mtext("n", side=1, line=3, cex=1.3) 

dev.off()


# draw HDBs for each plot
Nmean <- matrix(unlist(Nmean), nrow=k)
allHDB <- list()
for (i in 1:ncol(Nmean)) {
    dat <- Nmean[, i]
    varName <- "variable"
    allHDB[[length(allHDB)+1]] <- HDB(dat, varName, limits=c(2.62, 2.78))
}
pdf(paste0("3(a)_meanHDB_k", k, ".pdf"), width = 16, height = 15)
grid.arrange(grobs = allHDB, ncol = 2, nrow = 3)
dev.off()


# export the estimated results of E(Nbar) and se(N)
Nexp <- data.frame(Nmeanavg, Nsd)
TheoNse <- sqrt(3*exp(1)-exp(2)) / sqrt(n)
Nexp <- cbind(Nexp, TheoNse)
rownames(Nexp) <- c("n=1,000", "n=2,000", "n=5,000", "n=10,000", "n=100,000")
write.csv(Nexp, "3(a)-estimates.csv")


stop()

# test execution time for generating N (rN.2 is much faster)
beg <- Sys.time()
rN.1(100000)
end <- Sys.time()
cat("rN.1, time: ", beg-end, "\n")

beg <- Sys.time()
rN.2(100000)
end <- Sys.time()
cat("rN.2, time: ", beg-end, "\n")

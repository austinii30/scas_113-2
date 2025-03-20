
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


Nsd <- sapply(Nmean, FUN = sd)


# calculate the mean of Nbar
Nmeanavg <- sapply(Nmean, FUN = mean)
pdf("3.a_mean(Nbar).pdf", width=7, height=5)

par(mar = c(4, 6.3, 1, 1))
par(mgp = c(4, 1, 0))
plot(y=Nmeanavg, x=1:5,
     ylab="", xlab="", main="",
     pch=21, col="black", bg="orange", 
     cex=1.3, cex.lab=1.3, cex.axis=1.3,
     yaxt="n", xaxt="n", las=1, bty="l")
# custom x-axis labels and axis
axis(1, at=1:5, labels=c("1,000", "2,000", "5,000", "10,000", "100,000"), cex.axis=1.3)  # Custom x-axis labels
mtext("n", side=1, line=3, cex=1.3) 
# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.3, cex.lab=1.3)
mtext(expression(paste("Mean of ", bar("N"))), side=2, line=4.8, cex=1.3) 

dev.off()


# plot the distribution of Nbar


# plot sd(N) ~ n
pdf("3.a_sd(N).pdf", width=7, height=5)

par(mar = c(4, 5.2, 1, 1))
par(mgp = c(4, 1, 0))
plot(y=Nsd, x=1:5, 
     ylab="Standard Error of N", xlab="", main="",
     pch=21, col="black", bg="skyblue", 
     cex=1.3, cex.lab=1.3, cex.axis=1.3,
     xaxt="n", las=1, bty="l")
axis(1, at=1:5, labels=c("1,000", "2,000", "5,000", "10,000", "100,000"), cex.axis=1.3)  # Custom x-axis labels
mtext("n", side=1, line=3, cex=1.3) 

dev.off()







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


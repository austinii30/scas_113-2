########################
# Question 4 for hw3
########################

set.seed(2025)

nIntervals <- 1000
nBootstraps <- 500
bsSize <- 300
upper <- 0.975
lower <- 0.025

nSample <- rnorm(25)
tSample <- rt(n=25, df=5)
mean(nSample)
sd(nSample)
mean(tSample)
sd(tSample)

res <- matrix(NA, nrow=nIntervals, ncol=8)
containMean <- c()

for (expType in 1:4) {
    cat("Executing type ", expType, ".\n", sep="")
    dat <- nSample
    if (expType > 2) dat <- tSample

    parametric <- TRUE
    if (expType %% 2 == 0) parametric <- FALSE

    # get 1000 intervals
    meanWithin <- 0
    for (i in 1:nIntervals) {
        cat("Interval ", i, "\r", sep="")
        # do 500 bootstraps
        if (parametric)
            bsSample <- 
                rnorm(bsSize*nBootstraps, mean=mean(dat), sd=sd(dat))
        else
            bsSample <- 
                sample(x=dat, size=(bsSize*nBootstraps), replace=T)

        bsSample <- matrix(bsSample, ncol=bsSize, nrow=nBootstraps)
        bsAvg <- sort(rowMeans(bsSample))
        pctg <- order(bsAvg) / length(bsAvg)
        ub <- bsAvg[which(pctg >= upper)[1]]
        lb <- bsAvg[which(pctg >= lower)[1]]

        res[i, c(expType*2-1, expType*2)] <- c(ub, lb)
        if (ub >= 0 & 0 >= lb) meanWithin <- meanWithin + 1
    }
    cat("\n")
    containMean <- c(containMean, meanWithin)
}

colnames(res) <- 
    c("Norm_Para_UB", "Norm_Para_LB", "Norm_NPar_UB", "Norm_NPar_LB", 
      "t_Para_UB", "t_Para_LB", "t_NPar_UB", "t_NPar_LB")

write.csv(res, "./q4_simRes.csv")

names(containMean) <- c("Norm_Para", "Norm_NPar", "t_Para", "t_NPar")
containMean


# data
dat <- nSample

pdf("q4-normal-T&S.pdf", width=7, height=5)

# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.text <- 1.6
par(mar = c(4, 5.5, 1, 1))

# density for the data
dens <- density(dat)
plot(dens, las=1, bty="l", 
     xlab="X", ylab="", main="",
     cex=2, cex.lab=2, cex.axis=1.8,
     col="blue", lwd=1.8, lty=2,
     yaxt="n", xlim=c(-4, 4), ylim=c(0, 0.5),
     panel.first = {
       # Add grid lines; by default, use axis tick locations.
       grid(nx=NULL, ny=NULL, col="gray90", lty=1, lwd=1)
    })

# custom y-axis labels and axis
#axis(2)  # no customization
axis(2, las=1, cex.axis=cex.axis)
mtext("Density", side=2, line=3.8, cex=cex.lab) 

# theoretical density
curve(dnorm(x), from=-5, to=5, add=TRUE, col="red")

# legend
legend("topright", legend=c("Theoretical", "Sample"), 
       col=c("red", "blue"), lty=c(1, 2), lwd=1.5, 
       bty="n", cex=cex.text)

dev.off()


# data
dat <- tSample

pdf("q4-t5-T&S.pdf", width=7, height=5)

# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.text <- 1.6
par(mar = c(4, 5.5, 1, 1))

# density for the data
dens <- density(dat)
plot(dens, las=1, bty="l", 
     xlab="X", ylab="", main="",
     cex=2, cex.lab=2, cex.axis=1.8,
     col="blue", lwd=1.8, lty=2,
     yaxt="n", ylim=c(0, 0.5), xlim=c(-4, 4),
     panel.first = {
       # Add grid lines; by default, use axis tick locations.
       grid(nx=NULL, ny=NULL, col="gray90", lty=1, lwd=1)
    })

# custom y-axis labels and axis
#axis(2)  # no customization
axis(2, las=1, cex.axis=cex.axis)
mtext("Density", side=2, line=3.8, cex=cex.lab) 

# theoretical density
curve(dt(x, df=5), from=-5, to=5, add=TRUE, col="red")

# legend
legend("topright", legend=c("Theoretical", "Sample"), 
       col=c("red", "blue"), lty=c(1, 2), lwd=1.5, 
       bty="n", cex=cex.text)

dev.off()

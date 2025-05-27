###########################################################
# Filename: q4.R
# Course  : 113-2 Statistical Computation and Simulation, HW5, 
#           Question 4 
# Author  : Potumas Liu 
# Date    : 2025/05/27
###########################################################


library(MCMCpack)

library(showtext)  # display chinese in .pdf files
showtext_auto()


dat <- read.csv("dat-q4.csv")
#str(dat)
riders <- dat$riders_registered
temp_feel <- dat$temp_feel

# no missing values
print(nrow(dat[!complete.cases(dat), ]))

# 1. fit the linear regression by OLS
model <- lm(riders ~ temp_feel, data=dat)

summary(model)


xlab <- "體感溫度"
ylab <- "腳踏車比賽報名人數"


# 2. estimation by MCMC
# b0, B0, sigma.mu, sigma.var
posterior <- MCMCregress(riders ~ temp_feel, data=dat, 
                         burnin=20000, mcmc=30000)# verbose=5000)
pdf("q4-mcmcRes.pdf")
plot(posterior)
dev.off()

raftery.diag(posterior)
summary(posterior)

#str(summary(posterior))
a0 <- summary(posterior)$statistics["(Intercept)", "Mean"]
b0 <- summary(posterior)$statistics["temp_feel", "Mean"]


# 3. plot the scatter plot and the regression
# --- Line Plot --------------------------------------------
# Purpose: Line plot for 1/2 column, A4 size.
# ----------------------------------------------------------
# data
yDat <- riders
xDat <- temp_feel

pdf("q4-lm-plot.pdf", width=7, height=5, fonts="CNS1")

# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.dot  <- 0.7
par(mar=c(4.5, 6.5, 1, 1), mgp=c(3, 1.2, 0))

plot(
    y=yDat, x=xDat,
    ylab="", xlab=xlab, main="", 
    las=1, bty="l",
    pch=21, col="black", bg="skyblue", type="p",
    cex=cex.dot, cex.lab=cex.lab, cex.axis=cex.axis,
    yaxt="n", #xaxt="n",  # to custom axis
    panel.first = {
        # Add grid lines; by default, use axis tick locations.
        grid(nx=NULL, ny=NULL, col="gray80", lty=1, lwd=1)
    })

abline(model, lwd=2, col="red")

# custom x-axis labels and axis
#axis(1)  # no customization
#axis(1, at=1:length(xDat), labels=1:10, cex.axis=cex.axis)  
#mtext("n", side=1, line=3.3, cex=cex.lab)  # 'xlab'
# line: distance between the xlab and the x-axis

# custom y-axis labels and axis
axis(2, las=1, cex.axis=cex.axis)
mtext(ylab, side=2, line=4.9, cex=cex.lab) 

dev.off()
# ---------------------------------------------------------


# --- Line Plot --------------------------------------------
# Purpose: Line plot for 1/2 column, A4 size.
# ----------------------------------------------------------
# data
yDat <- riders
xDat <- temp_feel

pdf("q4-MCMCreg-plot.pdf", width=7, height=5, fonts="CNS1")

# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.dot  <- 0.7
par(mar=c(4.5, 6.5, 1, 1), mgp=c(3, 1.2, 0))

plot(
    y=yDat, x=xDat,
    ylab="", xlab=xlab, main="", 
    las=1, bty="l",
    pch=21, col="black", bg="skyblue", type="p",
    cex=cex.dot, cex.lab=cex.lab, cex.axis=cex.axis,
    yaxt="n", #xaxt="n",  # to custom axis
    panel.first = {
        # Add grid lines; by default, use axis tick locations.
        grid(nx=NULL, ny=NULL, col="gray80", lty=1, lwd=1)
    })

abline(a=a0, b=b0, lwd=2, col="red")

# custom x-axis labels and axis
#axis(1)  # no customization
#axis(1, at=1:length(xDat), labels=1:10, cex.axis=cex.axis)  
#mtext("n", side=1, line=3.3, cex=cex.lab)  # 'xlab'
# line: distance between the xlab and the x-axis

# custom y-axis labels and axis
axis(2, las=1, cex.axis=cex.axis)
mtext(ylab, side=2, line=4.9, cex=cex.lab) 

dev.off()
# ---------------------------------------------------------

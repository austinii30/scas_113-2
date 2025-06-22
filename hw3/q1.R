########################
# Question 1 for hw3
########################

library(showtext)  # display chinese in .pdf files
showtext_auto()

source("../../myPKG/R/dfManipulation.R")

# data
dataPath <- "./Q1_data.csv"
data <- read.csv(dataPath)
data <- recognizeRowNames(data, "year")
train <- data[1:15,  ]
test  <- data[16:20, ]

# PCA
data.cor <- cor(train)
data.eigen <- eigen(data.cor)

data.eigen
stop()

ev <- data.eigen$values
ev <- ev/ncol(train)
varExplained <- 1 - ev

# plot explained variation for each component
pdf("./q1_eigenVarExplained.pdf", width=7, height=5, fonts="CNS1")

yDat <- 1-varExplained
xDat <- 1:length(varExplained)
# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.dot  <- 2
par(mar=c(4.5, 6.5, 1, 1), mgp=c(3, 1.2, 0))

plot(
    y=yDat, x=xDat,
    ylab="", xlab="", main="", 
    las=1, bty="l",
    pch=21, col="black", bg="orange", type="o",
    cex=cex.dot, cex.lab=cex.lab, cex.axis=cex.axis,
    yaxt="n", xaxt="n",  # to custom axis
    panel.first = {
        grid(nx=NULL, ny=NULL, col="gray80", lty=1, lwd=1)
    })
# custom x-axis labels and axis
axis(1, at=1:length(xDat), labels=xDat, cex.axis=cex.axis)  
mtext("主成份", side=1, line=3.3, cex=cex.lab)  # 'xlab'
# custom y-axis labels and axis
axis(2, las=1, cex.axis=cex.axis)
mtext("解釋總變異比例", side=2, line=4.5, cex=cex.lab) 
dev.off()

# plot eigenvalues vs. 1
pdf("./q1_eigenValues.pdf", width=7, height=5, fonts="CNS1")

yDat <- data.eigen$values
xDat <- 1:length(xDat)
# key parameters
cex.lab  <- 2
cex.axis <- 1.8
cex.dot  <- 2
par(mar=c(4.5, 6.5, 1, 1), mgp=c(3, 1.2, 0))

plot(
    y=yDat, x=xDat,
    ylab="", xlab="", main="", 
    las=1, bty="l",
    pch=21, col="black", bg="orange", type="o",
    cex=cex.dot, cex.lab=cex.lab, cex.axis=cex.axis,
    yaxt="n", xaxt="n",  # to custom axis
    panel.first = {
        grid(nx=NULL, ny=NULL, col="gray80", lty=1, lwd=1)
    })
# custom x-axis labels and axis
axis(1, at=1:length(xDat), labels=xDat, cex.axis=cex.axis)  
mtext("主成份", side=1, line=3.3, cex=cex.lab)  # 'xlab'
# custom y-axis labels and axis
axis(2, las=1, cex.axis=cex.axis)
mtext("特徵值", side=2, line=4.5, cex=cex.lab) 
abline(a = 1, b = 0, col="red", lwd=3)
dev.off()

# error for 1~5 components
errors <- matrix(NA, nrow=17, ncol=2)

for (i in 1:16) {
    pc <- data.eigen$vectors[, 1:i]
    newCords <- data.matrix(test) %*% pc
    
    # reconstruct the testing data by PCs of training data
    reconstruct <- matrix(nrow=0, ncol=17)
    for (r in 1:nrow(newCords)) {
        vec <- newCords[r, ]
        res <- pc %*% matrix(vec, ncol=1)
        reconstruct <- rbind(reconstruct, as.vector(res))
    }
    
    se.pca <- sum((reconstruct - test)^2) / 60
    errors[i, 1] <- se.pca
    cat("Mean Squared Error of PCA (", i," PCs): ", se.pca, "\n", sep="")
}


# SVD
svdRes <- svd(train) 

for (i in 1:ncol(svdRes$v)) {
    v <- svdRes$v[, 1:i]
    newCords.svd <- data.matrix(test) %*% v
    
    reconstruct <- matrix(nrow=0, ncol=17)
    for (r in 1:nrow(newCords.svd)) {
        vec <- newCords.svd[r, ]
        res <- v %*% matrix(vec, ncol=1)
        reconstruct <- rbind(reconstruct, as.vector(res))
    }
    
    se.svd <- sum((reconstruct - test)^2) / 60
    errors[i, 2] <- se.svd
    cat("Mean Squared Error of SVD (", i, " DIMs): ", se.svd, "\n", sep="")
}

rownames(errors) <- 1:17
colnames(errors) <- c("MSE.pca", "MSE.svd")
errors
log10(errors)

# plot mse of PCA & SVD for 1~5 components
dat1 <- log10(errors[, 1])
dat2 <- log10(errors[, 2])
k <- length(dat1)

pdf("./q1_mse.pdf", width = 10, height = 6, fonts="CNS1")

cex.lab  <- 1.7
cex.axis <- 1.5
cex.dot  <- 1.5
par(mar = c(4, 4.5, 0.5, 1.5), mgp = c(2.5, 1, 0))  

plot(x=(1:k), y=dat1, type="n", bty="l", pch=16, 
    xlab = "主成份/維度個數", las=1,
    ylab = expression(paste("平均平方誤差（取", "log"[10], "）")),
    ylim = c(2, 8), 
    cex.axis=cex.axis, cex.lab=cex.lab,
    panel.first = {
       # Add grid lines; by default, use axis tick locations.
       grid(nx=NULL, ny=NULL, col="gray90", lty=1, lwd=1)
    })

lines(x=(1:k), y=dat1, type="o", pch=21, bg="orange", cex=cex.dot)
lines(x=(1:k), y=dat2, type="o", pch=24, bg="skyblue", cex=cex.dot)
legend("bottomleft", legend=c("PCA", "SVD"),
       lty=c(1, 1), pch=c(21, 24), pt.bg=c("orange", "skyblue"),
       bg="white", cex=cex.lab, pt.cex=cex.dot)

dev.off()

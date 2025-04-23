########################
# Question 1 for hw3
########################

library(EFAtools)
#library(Rfast) # for colVars

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

ev <- data.eigen$values
ev <- ev/ncol(train)
varExplained <- round(cumsum(ev), 2)
# plot cumulative explained variation
plot(varExplained, type="o")
# plot eigenvalues vs. 1
plot(data.eigen$values, pch = 19, type = "o", col="blue",
     ylim = c(-0.5, 5),
     xlab = "Principal Component", ylab = "Eigenvalue")
abline(a = 1, b = 0, col="red")

# 3 components
pc3 <- data.eigen$vectors[, 1:3]
pc <- pc3
newCords <- data.matrix(test) %*% pc

# reconstruct the testing data by PCs of training data
reconstruct <- matrix(nrow=0, ncol=17)
for (r in 1:nrow(newCords)) {
    vec <- newCords[r, ]
    res <- pc %*% matrix(vec, ncol=1)
    reconstruct <- rbind(reconstruct, as.vector(res))
}

se.pca <- sum((reconstruct - test)^2) / 60
cat("Mean Squared Error of PCA: ", se.pca, "\n", sep="")


# SVD
#mu <- colMeans(train)
#ct <- scale(train, center = TRUE, scale = FALSE)

#for (i in 1:length(mu)) {
#    mmm <- mu[i]
#    test[, i] <- test[, i] - mmm
#}

svdRes <- svd(train) 
v3 <- svdRes$v[, 1:3]

newCords.svd <- data.matrix(test) %*% v3

reconstruct <- matrix(nrow=0, ncol=17)
for (r in 1:nrow(newCords.svd)) {
    vec <- newCords.svd[r, ]
    res <- v3 %*% matrix(vec, ncol=1)
    reconstruct <- rbind(reconstruct, as.vector(res))
}

se.svd <- sum((reconstruct - test)^2) / 60
cat("Mean Squared Error of SVD: ", se.svd, "\n", sep="")

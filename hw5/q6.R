###########################################################
# Filename: q6.R
# Course  : 113-2 Statistical Computation and Simulation, HW5, 
#           Question 6 
# Author  : Potumas Liu 
# Date    : 2025/05/27
###########################################################


library(MCMCpack)

data("birthwt", package="MASS")

dat <- birthwt
dat <- dat[, -10]

# low, race, smoke, ht, ui
factorVar <- c(4, 5, 7, 8)
for (idx in factorVar) { 
    dat[, idx] <- as.factor(dat[, idx])
}

str(dat)

for (idx in 1:ncol(dat)) {
    cat(strrep("-", 30), "\n")
    print(colnames(dat)[idx])
    print(table(dat[, idx]))
    cat(strrep("-", 30), "\n\n")
}

# no missing values
print(nrow(dat[!complete.cases(dat), ]))


glmfit <- glm(low ~ ., binomial, dat)

summary(glmfit)


posterior <- MCMClogit(low ~ ., data=dat, 
                       burnin=20000, mcmc=30000)
                       #verbose=5000)

pdf("q6-mcmcRes.pdf")
plot(posterior)
dev.off()

raftery.diag(posterior)
summary(posterior)

mcmcest <- summary(posterior)$statistics
print(mcmcest)

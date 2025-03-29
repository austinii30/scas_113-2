###########################################################
# Filename: main_1.R
# Course  : 113-2 Statistical Computation and Simulation, HW2, 
#           Question 1
# Author  : Potumas Liu 
# Date    : 2025/03/24
###########################################################

source("../../myPKG/R/figure.R")
source("../../myPKG/R/EDA.R")
source("../../myPKG/R/dfManipulation.R")
source("code.R")


# 1.(a)
set.seed(2025)

k <- 1000
#n <- 1000
#n <- 100
#n <- 200
n <- 500
theta <- c(0, 0.05, 0.1, 0.15, 0.2)

# for each theta, draw k random samples, each with n=100
simDat <- list()
for (th in theta) {
    dat <- c()
    for (i in 1:k)
        dat <- c(dat, arima.sim(n=n, list(ar=c(th, th))))
    # a column is a sample (k columns in total)
    simDat[[length(simDat)+1]] <- matrix(dat, ncol=k)
}


# correlation
simDatCor <- list()
for (dat in simDat) {
    # returns a matirx with 2 rows, each is a correlation
    simCor <- sapply(1:ncol(dat), FUN = function(sIdx) {
            sample <- dat[, sIdx]
            cor1 <- cor(sample[-1], sample[-n])
            cor2 <- cor(sample[-c(1, 2)], sample[-c(n-1,n)])
            return(c(cor1, cor2))
        }) 
    # transpose so row: a sample, col: cor1 & cor2
    simDatCor[[length(simDatCor)+1]] <- t(simCor)
}
names(simDatCor) <- c("0", "0.05", "0.1", "0.15", "0.2")


# summary statistics for autocorrelations 
d <- simDatCor[[1]]
for (i in 2:length(simDatCor))
    d <- cbind(d, simDatCor[[i]])
colnames(d) <- c("0_1", "0_2", "0.05_1", "0.05_2", 
                 "0.1_1", "0.1_2", "0.15_1", "0.15_2", 
                 "0.2_1", "0.2_2")
d <- data.frame(d)
ss <- summaryStatistics(d)
for (i in 1:ncol(ss)) ss[, i] <- round(ss[, i], 3)
write.csv(ss, paste0("1(a)_summaryStatistics_n", n, "_k", k, ".csv"))


# draw HBDs 
allHBD <- list()
for (i in 1:ncol(d)) {
    dat <- d[, i]
    varName <- colnames(d)[i]

    if (i %% 2 == 0) p <- HDB(dat, varName, "Dark2", limits=c(-0.5, 0.5))
    else             p <- HDB(dat, varName, limits=c(-0.5, 0.5))
    allHBD[[length(allHBD)+1]] <- p
}
pdf(paste0("1(a)_corHBD_n", n, "_k", k, ".pdf"), width = 18, height = 25)
grid.arrange(grobs = allHBD, ncol = 2, nrow = 5)
dev.off()

warnings()



# 1.(b)
a <- 0.05
pRes <- matrix(nrow=3, ncol=0)
pNAs <- matrix(nrow=3, ncol=0)
for (i in 1:length(simDat)) {
    dat <- simDat[[i]]
    p.gap <- apply(dat, 2, FUN = gapTest)
    p.uad <- apply(dat, 2, FUN = uadTest)
    p.pmt <- apply(dat, 2, FUN = pmtTest)

    pRes <- cbind(pRes, c(sum(p.gap<a, na.rm=TRUE), 
                          sum(p.uad<a, na.rm=TRUE), 
                          sum(p.pmt<a, na.rm=TRUE)))
    pNAs <- cbind(pNAs, c(sum(is.na(p.gap)), sum(is.na(p.uad)), sum(is.na(p.pmt))))
}
colnames(pRes) <- colnames(pNAs) <- names(simDatCor)
rownames(pRes) <- rownames(pNAs) <- c("Gap", "UaD", "Pmt")

write.csv(pRes, paste0("1(b)_rejH0Times_n", n, "_k", k, "_a", a, ".csv"))
write.csv(pNAs, paste0("1(b)_cannotTest_n", n, "_k", k, "_a", a, ".csv"))

warnings()

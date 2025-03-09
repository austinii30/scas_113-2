############################################################ 
# Filename: main.R
# Purpose : 113-2 Statistical Computation and Simulation
# Author  : Potumas Liu
# Date    : 2025/03/06
############################################################ 

library(showtext)  # display chinese in .pdf files
showtext_auto()

library(randtoolbox)  # for Gap test
source("./code.R")

library(parallel)


# 1 (a)
# ================================================================
x.1a <- rm_1a(10000, rmSeed=2025)
print(ks.test(x.1a, "punif"))


# 1 (b)
# ================================================================
x.1b <- rm_1b(10000, rmSeed=2025)
print(ks.test(x.1b, "punif"))


pdf("1a-hist.pdf", width = 8, height = 6, fonts="CNS1")
#pdf("1b-hist.pdf", width = 8, height = 6, fonts="CNS1")

par(mar = c(4.5, 4.5, 0.5, 1.5))  # Lower the bottom margin (the 1st value)
par(font.main = 1)  # 1: normal, 2: bold
par(mgp = c(3.2, 1, 0))

data <- x.1a
#data <- x.1b

hist(data, prob=T,
    las=1, col="orange",
    main="", xlab=expression(x[i]), ylab="相對次數",
    cex.main=3, cex.axis=1.5, cex.lab=1.6)

dev.off()



# 2 (a)
# ================================================================
alpha <- 0.05
maxM <- 50
n <- 10000
x.2a <- list()
for (i in 1:maxM)
    x.2a[[length(x.2a)+1]] <- fibN(n, m=i, rmSeed=2025)

testResult <- sapply(x.2a, FUN = function(x) {
        # uniform test
        uni.p <- ks.test(x, "punif")$p

        # independent tests
        gap.p <- gapTest(0.35, 0.8, x)$p.value
        uad.p <- uadTest(x)
        
        res <- c(uni.p, gap.p, uad.p)

        return(res)
    })

testResult <- t(testResult)
colnames(testResult) <- c("uni.p", "gap.p", "uad.p")

save(testResult, 
     file=paste("2a", "_maxM", maxM, "_n", n, 
                "_testResults.RData"))

uniTest <- testResult[, "uni.p"]     < alpha

gapTest.05 <- testResult[, "uad.p"]  < alpha
uadTest.05 <- testResult[, "uni.p"]  < alpha
gapTest.025 <- testResult[, "uad.p"] < (alpha/2)
uadTest.025 <- testResult[, "uni.p"] < (alpha/2)


stop()



# 2 (b)
# ================================================================
#n <- 10000
n <- 100000
#nIter <- 100
nIter <- 1000
seed <- 2025
#maxK <- 10   
maxK <- 1000
alpha <- 0.05

exampleK <- c(2, 5, 10, 100, 500, 1000)

begTime <- Sys.time() 
comment(paste0("begin parallel, time: ", begTime))

func2b <- function(k) {
    library(showtext)  # display chinese in .pdf files
    showtext_auto()

    set.seed(seed)
    x.s <- x.c <- c()
    for (nIter in 1:nIter) {
        x.s <- c(x.s, sample(1:k, n, replace=T))
        x.c <- c(x.c, ceiling(runif(n, min=0, max=k)))
    }
    x.s <- matrix(x.s, nrow = n)
    x.c <- matrix(x.c, nrow = n)
    x.all <- cbind(x.s, x.c)

    pValues <- c()
    for (i in 1:ncol(x.all)) {
        p <- chisq.test(table(x.all[, i]))$p.value
        pValues <- c(pValues, p)
    }
    pValues <- matrix(
        pValues, ncol=2, 
        dimnames = list(iter=1:nIter, method=c("sample", "ceiling"))
    )

    rejects <- pValues < alpha
    rej.s <- sum(rejects[, "sample"])
    rej.c <- sum(rejects[, "ceiling"])

    #return(c(rej.s, rej.c))

    if (k %in% exampleK) {
        acpIdx.s <- which(rejects[, "sample"]==F)[1]
        acpIdx.c <- which(rejects[, "ceiling"]==F)[1]
        rejIdx.s <- which(rejects[, "sample"]==T)[1]
        rejIdx.c <- which(rejects[, "ceiling"]==T)[1]
        idxVec <- c(acpIdx.s, acpIdx.c, rejIdx.s, rejIdx.c)
        if (sum(is.na(idxVec))!=0) stop(idxVec)

        dat.acp.s <- x.s[, acpIdx.s] 
        dat.acp.c <- x.c[, rejIdx.s]
        dat.rej.s <- x.s[, acpIdx.c]
        dat.rej.c <- x.c[, rejIdx.c]
        datList <- list(dat.acp.s, dat.acp.c, dat.rej.s, dat.rej.c)
        datIdx <- c("dat.acp.s", "dat.acp.c", "dat.rej.s", "dat.rej.c")

        for (i in 1:length(datIdx)) {
            dat <- datList[[i]]; datName <- datIdx[i]
            idx <- idxVec[i]

            if (grepl(".s", datName)) method <- "sample"
            else method <- "ceiling"

            if (grepl(".rej", datName)) decision <- "rej"
            else decision <- "acp"

            fileName <- paste0("2b_k=", k, "_", method, "_", decision, "_idx", idx,".pdf")
            pdf(fileName, width = 8, height = 6, fonts="CNS1")

            par(mar = c(4.5, 4.5, 0.5, 1.5))
            par(mgp = c(3.2, 1, 0))
            
            save(datList, file = paste0(fileName, ".RData"))

            if (k <= 10) {
                barplot(table(dat), prob=T,
                    las=1, col="orange",
                    main="", xlab=expression(x[i]), ylab="相對次數",
                    cex.main=3, cex.axis=1.5, cex.lab=1.6)
            } else {
                hist(dat, prob=T,
                    las=1, col="orange",
                    main="", xlab=expression(x[i]), ylab="相對次數",
                    cex.main=3, cex.axis=1.5, cex.lab=1.6)
            }

            dev.off()
        }
    }

    return(c(rej.s, rej.c))
}

# make a cluster, pass all objects in the global into the cluster
coreNum <- ifelse(detectCores()-2 > 1, detectCores()-2, 1)
c1 <- makeCluster(coreNum)  # specify number of nodes
clusterExport(c1, ls()) 

# 'parLapply(*cluster, *inputs, *function)'
allKResults <- parLapply(c1, 2:maxK, func2b)

stopCluster(c1)  # stop the cluster, free up space (important!!!)
endTime <- Sys.time()

# export results
res.s <- res.c <- c()
for (i in 1:length(allKResults)) {
    res.s <- c(res.s, allKResults[[i]][1])
    res.c <- c(res.c, allKResults[[i]][2])
}
names(res.s) <- 2:maxK
names(res.c) <- 2:maxK
clusterRes <- list(res.s = res.s, res.c = res.c, 
                   timeSpent = round(endTime-begTime, 2))

RDataName <- "2(b)"
save(clusterRes, file = paste0(
    "./", RDataName, "_seed", seed, "_n", n, 
    "_nIter", nIter, "_maxK", maxK, "_alpha", alpha, 
    ".RData"))

# plot the results
pdf("2bPlot.pdf", width = 10, height = 6, fonts="CNS1")

par(mar = c(4.5, 4.5, 0.5, 1.5))  # Lower the bottom margin (the 1st value)
par(mgp = c(2.5, 1, 0))  

plot(x=2:maxK, y=res.s, type="n", bty="l", pch=16, 
    xlab = "k",
    ylab = expression(paste("拒絕", "H"[0], "次數（共", 10, "次）")),
    ylim = c(0, max(res.s, res.c)), 
    cex.axis=1.5, cex.lab=1.5,
    panel.first = {
       # Draw a rectangle covering the entire plotting region with a light gray color.
       rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "white", border = NA)
       # Add grid lines; by default, use axis tick locations.
       grid(nx=NULL, ny=NULL, col="gray90", lty=1, lwd=1)
    }
)

lines(x=2:maxK, y=res.s, col="red", type="o", pch=16)
lines(x=2:maxK, y=res.c, col="blue", type="o", pch=16)
legend("topleft", legend=c("sample()", "ceiling(runif())"),
       col=c("red", "blue"), lty=c(1, 1), pch=c(16, 16), 
       bg="white", cex=1.5, pt.cex=1)

dev.off()

stop()

# ---------- old -------------------------------------------------
#for (k in K) {
#    x.s <- c(x.s, rm_2b(n, k, "sample", rmSeed=2025))
#    x.c <- c(x.c, rm_2b(n, k, "ceiling", rmSeed=2025))
#}

# each column is a random number sample of size n
#x.s <- matrix(x.s, nrow = n)
#x.c <- matrix(x.c, nrow = n)
#x.all <- cbind(x.s, x.c)

#pValues <- c()
#for (i in 1:ncol(x.all)) {
#    p <- chisq.test(table(x.all[, i]))$p.value
#    pValues <- c(pValues, p)
#}
#pValues <- matrix(pValues, ncol=2, 
#                  dimnames = list(k=K, 
#                                  method=c("sample", "ceiling")))
#print(pValues)
#
#stop()
# ---------- old -------------------------------------------------

# graphical tools
par(mfrow = c(2, 3))
for (k in exampleK) {
}

# 3 
# ================================================================
N <- 1000
n10 <- 10; n50 <- 50; n100 <- 100

res <- c()
#for () {
#    res <- c(res, rnorm((n10+n50+n100)*N))

#}

t10.10  <- rt(10, df = 10)
t50.10  <- rt(50, df = 10)
t100.10 <- rt(100, df = 10)
t10.20  <- rt(10, df = 20)
t50.20  <- rt(50, df = 20)
t100.20 <- rt(100, df = 20)

















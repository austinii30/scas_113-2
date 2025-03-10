n <- 100000
nIter <- 100
#nIter <- 1000
seed <- 2025
#maxK <- 10   
maxK <- 1000
alpha <- 0.05


library(showtext)  # display chinese in .pdf files
showtext_auto()

load("./2(b)_seed2025_n1e+05_nIter100_maxK1000_alpha0.05.RData")

res.s <- clusterRes[[1]]
res.c <- clusterRes[[2]]

# plot the results
pdf("2bPlot.pdf", width = 10, height = 6, fonts="CNS1")

par(mar = c(4.5, 4.5, 0.5, 1.5))  # Lower the bottom margin (the 1st value)
par(mgp = c(2.5, 1, 0))  

plot(x=2:maxK, y=res.s, type="n", bty="l", pch=16, 
    xlab = "k",
    ylab = expression(paste("拒絕", "H"[0], "次數（共檢定", 100, "次）")),
    ylim = c(0, 15), 
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


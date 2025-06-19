###########################################################
# Filename: analyze.R
# Purpose : Analyze the organized resutls
# Author  : Austin Liu
# Date    : 2025/06/18
###########################################################

source(file.path(this.path::this.dir(), "utils.R"))
import("func.R")
import("myGWRfunc.R")

library(showtext)  # display chinese in .pdf files
showtext_auto()

sigmas <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
tmodels <- c("LM", "GWR")
nsub <- 1000
alpha <- 0.05
aicth <- 3
secth <- 1

MoransI    <- LRT <- AICc <- matrix(NA, nrow=2, ncol=length(sigmas))
MonteCarlo <- SEC <- matrix(NA, nrow=8, ncol=length(sigmas))

colnames(MoransI) <- colnames(LRT) <- 
    colnames(AICc) <- colnames(MonteCarlo) <- sigmas
rownames(MoransI) <- rownames(LRT) <- rownames(AICc) <- tmodels
rownames(MonteCarlo) <- c("LM.0", "LM.1", "LM.2", "LM.3",
                          "GWR.0", "GWR.1", "GWR.2", "GWR.3")

sink("Analyze.txt", split=TRUE)
for (tmidx in 1:length(tmodels)) {
    pdf(outpath(paste0(tmidx, "-truemodel.pdf")), width=12, height=16)  # 13.3
    par(mar = c(3, 4, 0.5, 1), mgp = c(4, 1, 0), mfrow=c(7, 4))
for (sidx in 1:length(sigmas)) {
    s <- sigmas[sidx]

    #pdf(outpath(paste0(s, "-sec.pdf")), width=12, height=3.8)
    #, fonts="CNS1")
    #par(mar = c(3, 4, 0.5, 1), mgp = c(4, 1, 0), mfrow=c(2, 4))

        tm <- tmodels[tmidx]
        load(file=outpath(paste0("organized-", tm, "-sigma-", s, ".RData")))
        dat <- combres

        # bandwidth
        
        # AICc differences (global - local)
        AICc[tmidx, sidx] <- sum(dat$aiccdiff > aicth)

        # Moran's I p-value
        MoransI[tmidx, sidx] <- sum(dat$morani < alpha)
        #moranires <- sum(dat$morani < alpha)
        #MoransI[tmidx, sidx] <- moranires
        #cat("[Moran's I]  " moranires, "  ( ", pctg(moranires/nsub), " )\n", sep="")

        # Monte-Carlo p-value
        from <- (tmidx-1)*4 + 1; to <- tmidx*4
        MonteCarlo[from:to, sidx] <- colSums(dat$montecarlo < alpha)
        #mcres <- colSums(dat$montecarlo < alpha)

        # LRT p-value
        LRT[tmidx, sidx] <- sum(dat$lrt < alpha)
        #lrtres <- sum(dat$lrt < alpha)
        #cat("[LRT]  " lrtres, "  ( ", pctg(lrtres/nsub), " )\n", sep="")

        # compare SE (local / global)
        SEdat <- dat$secompare
        for (cidx in 1:ncol(SEdat)) {
            sedat <- SEdat[, cidx]
            
            cex.axis <- 2
            cex.lab  <- 2
            
            hist(sedat, prob=TRUE,
                las=1, col="skyblue", 
                main="", xlab="", ylab="", xaxt="n", 
                cex.axis=cex.axis, cex.lab=cex.lab)
            box(bty="l",  # specifically for histograms
                panel.first = {
                    # Add grid lines; by default, use axis tick locations.
                    grid(nx=NULL, ny=NULL, col="gray80", lty=1, lwd=1)
                })
            hist(sedat, prob=TRUE,
                las=1, col="skyblue",
                main="", xlab="", ylab="", 
                cex.axis=cex.axis, cex.lab=cex.lab, add=TRUE)
            # custom y-axis labels and axis
            axis(1, cex.axis=cex.axis)
            mtext("", side=1, line=3, cex=cex.lab) 
            abline(v=1, col="red", lwd=2.5)
        }
        SEdat <- colSums(SEdat > secth)
        SEC[from:to, sidx] <- SEdat
    }
    dev.off()
    #stop()
}
sink()

write.csv(MoransI, outpath("Morans.I.csv"))
write.csv(MonteCarlo, outpath("MonteCarlo.csv"))
write.csv(SEC, outpath("seCompare.csv"))
write.csv(LRT, outpath("LRT.csv"))
write.csv(AICc, outpath("AICc.csv"))

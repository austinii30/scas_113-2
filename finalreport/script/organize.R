###########################################################
# Filename: organize.R
# Purpose : Organize the results of 'test.R'
# Author  : Austin Liu
# Date    : 2025/06/18
###########################################################

source(file.path(this.path::this.dir(), "utils.R"))
import("func.R")
import("myGWRfunc.R")

sigmas <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
tmodels <- c("LM", "GWR")

for (s in sigmas) {
    for (tm in tmodels) {
        resname <- paste0("simRes-", tm, "-1000-sigma-", s, ".RData")
        load(datpath(resname)); dat <- simres
        
        # bandwidth
        bw <- sapply(dat, `[[`, 1)

        # AICc differences (global - local)
        #print(dat)
        aiccdiff <- sapply(dat, `[[`, 2)
        #print(aiccdiff)
        aiccdiff <- aiccdiff[1, ] - aiccdiff[2, ]

        # Moran's I p-value
        morani <- sapply(dat, `[[`, 3)[4, ]

        # Monte-Carlo p-value
        montecarlo <- sapply(dat, function (d) { return(as.vector(d[[4]])) })
        montecarlo <- t(montecarlo)

        # LRT p-value
        lrt <- sapply(dat, `[[`, 5)[2, ]

        # compare SE (local / global)
        secompare <- sapply(dat, function (d) {
            mtx <- d[[6]]
            ratio <- mtx[, 1] / mtx[, 2]
            return(as.vector(ratio))
        })
        secompare <- t(secompare)

        combres <- list(bw=bw, aiccdiff=aiccdiff, morani=morani, montecarlo=montecarlo,
                           lrt=lrt, secompare=secompare)

        save(combres, file=outpath(paste0("organized-", tm, "-sigma-", s, ".RData")))
    }
}

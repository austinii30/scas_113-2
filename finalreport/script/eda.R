###########################################################
# Filename: eda.R
# Purpose : Examine generated spatial data
# Author  : Austin Liu
# Date    : 2025/06/14
###########################################################

source(file.path(this.path::this.dir(), "utils.R"))
import("func.R")

load(datpath("datlm.RData"))
load(datpath("datgwr.RData"))

xlim <- c(0, 10)
ylim <- c(0, 10)


# --------------------------------------------
# sample location
# --------------------------------------------
pdf(outpath("dat-location.pdf"))
for (i in 1:length(datlm)) {
    plot(x=datlm[[i]][, "x"], y=datlm[[i]][, "y"]) 
}
dev.off()


# --------------------------------------------
# LM: sample response
# --------------------------------------------
pdf(outpath("lm-response.pdf"), width=8, height=6.5)
for (i in 1:length(datlm)) {
    x <- datlm[[i]][, "x"]
    y <- datlm[[i]][, "y"]
    z <- datlm[[i]][, "response"]
    plot2d(x, y, z, xlim, ylim)
}
dev.off()


# --------------------------------------------
# GWR: sample response
# --------------------------------------------
pdf(outpath("gwr-response.pdf"), width=8, height=6.5)
for (i in 1:length(datgwr)) {
    x <- datgwr[[i]][, "x"]
    y <- datgwr[[i]][, "y"]
    z <- datgwr[[i]][, "response"]
    plot2d(x, y, z, xlim, ylim)
}
dev.off()


# --------------------------------------------
# GWR: covariates
# --------------------------------------------
for (ct in c("gp0", "gp1", "gp2", "gp3")) {
    pdf(outpath(paste0("gwr-", ct, ".pdf")), width=8, height=6.5)
    for (i in 1:length(datgwr)) {
        x <- datgwr[[i]][, "x"]
        y <- datgwr[[i]][, "y"]
        z <- datgwr[[i]][, ct]
        plot2d(x, y, z, xlim, ylim)
    }
    dev.off()
}

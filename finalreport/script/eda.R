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
zlim <- c(-1, 1)

gwrpar0 <- function(x, y) { return( (x-y)/10 ) }
gwrpar1 <- function(x, y) { return( (sin(x/0.5)+cos(y/0.8))/2 ) }
gwrpar2 <- function(x, y) { return( exp(-((x-2.5)^2+(y-3)^2)/2) + exp(-((x-8)^2+(y-6)^2)/4) - exp(-((x-7)^2+(y-2)^2)/0.8) - exp(-((x-4)^2+(y-8)^2)/1.5)) }
gwrpar3 <- function(x, y) { return( exp(-(y-2*(x-1))^2/3) + exp(-(3*(y-3)-x)^2/10) - 1) }
gwrparfuncs <- c(gwrpar0, gwrpar1, gwrpar2, gwrpar3)

# --------------------------------------------
# exact plots for each covariate
# --------------------------------------------
idx <- 0
for (gf in gwrparfuncs) {
    pdf(outpath(paste0("gwr-covariate-", idx, ".pdf")), width=8, height=6.5)
    x <- seq(0, 10, length.out=200)
    y <- seq(0, 10, length.out=200)
    z <- outer(x, y, gf)
    print(max(z))
    print(min(z))
    plot2d(x, y, z, xlim, ylim, zlim, ke=FALSE) 
    dev.off()
    idx <- idx + 1
}
stop()

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

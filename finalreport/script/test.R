###########################################################
# Filename: test.R
# Purpose : Test justification approaches 
# Author  : Austin Liu
# Date    : 2025/06/15
###########################################################

library(GWmodel)
#devtools::install_version("MuMIn", version = "1.47.1", repos = "https://cran.r-project.org")
library(MuMIn)  # for AICc()
library(ape)  # for moran()

source(file.path(this.path::this.dir(), "utils.R"))
import("func.R")
import("myGWRfunc.R")

load(datpath("datlm.RData"))
datlm <- lapply(datlm, data.frame) 
datlm.gw <- lapply(datlm, function(d) { coordinates(d) <- c("x", "y"); return(d) }) 

load(datpath("datgwr.RData")) 
datgwr <- lapply(datgwr, data.frame)
datgwr.gw <- lapply(datgwr, function(d) { coordinates(d) <- c("x", "y"); return(d) })

dat    <- list(datlm, datgwr)
dat.gw <- list(datlm.gw, datgwr.gw)

model <- formula(response ~ x1 + x2 + x3)
kernel <- "gaussian"

set.seed(2025)

results <- lapply(1:(2*length(datlm)), function(x) return(NULL))
datnames <- c("LM", "GWR")

# data for LM or GWR
for (i in 1:2) {
    # for each subdata
    #for (j in 1:length(datlm)) {
    #simres <- lapply(1:length(datlm), function(j) { 
    simres <- lapply(1:10, function(j) { 
        sink(outpath("junk"))
        # data for my code and GWmodel
        dlm <- dat[[ i ]][[ j ]]
        dgw <- dat.gw[[ i ]][[ j ]]

        # -----------------------------  
        # Fit global regression
        # -----------------------------  
        modellm <- lm(model, data=dlm)

        # -----------------------------  
        # Optimal bandwidth by AICc
        # -----------------------------  
        bwAICc <- bw.gwr(model, data=dgw, approach="AIC", kernel=kernel)

        # -----------------------------  
        # Fit GWR
        # -----------------------------  
        mygwr <- GWR(model, dlm, c("x", "y"), bwAICc, seCompare=TRUE, AICs=TRUE, doLRT=TRUE)
        #print(mygwr)

        gwgwr <- gwr.basic(model, dgw, bw=bwAICc, kernel=kernel)
        print(gwgwr)
        #print(str(gwgwr))

        # -----------------------------  
        # AICc scores
        # -----------------------------  
        AICc.gwr <- mygwr$AICs[1, 1]
        #AICc.gwr <- gwgwr$GW.diagnostic$AICc 
        AICc.lm  <- AICc(modellm)
        #print(AICc.lm)
        AICs <- c(AICc.lm, AICc.gwr); names(AICs) <- c("lm", "gwr")

        # -----------------------------  
        # BIC scores (unused)
        # -----------------------------  
        #bic <- bbmle::BIC(modellm)
        #print(bic)
        
        # -----------------------------  
        # Moran's I
        # -----------------------------  
        morani <- unlist(Moran.I(mygwr$residuals, invDist(distanceMtx(dlm[, 2:3]))))
        
        # -----------------------------  
        # Monte Carlo test (p-values)
        # -----------------------------  
        mcAICc <- gwr.montecarlo(model, dgw, 99, kernel, bw=bwAICc)

        # -----------------------------  
        # Approximate LRT
        # -----------------------------  
        lrt <- mygwr$LRT

        # -----------------------------  
        # Compare sigma
        # -----------------------------  
        seCompare <- mygwr$seCompare 
        
        res <- list(bwAICc=bwAICc, AICcs=AICs, MoransI=morani, montecarlo=mcAICc, LRT=lrt, seCompare=seCompare)
        #print(res); stop()
        idx <- (i-1)*length(datlm) + j
        sink()
        cat(datnames[i], ", data ", idx, " done.    \r", sep="")
        #print(res)
        #stop()
        return(res)
        #results[[idx]] <- res
    })
    cat("\n")
        
    #from <- (i-1)*length(datlm)
    #to   <- i*length(datlm)
    #results[from:to] <- simres
    save(simres, file=datpath(paste0("simRes-", datnames[i], ".RData")))
}

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
datlm.gw <- lapply(datlm, function(d) { rownames(d) <- 1:nrow(d); coordinates(d) <- c("x", "y"); return(d) }) 

load(datpath("datgwr.RData")) 
datgwr <- lapply(datgwr, data.frame)
datgwr.gw <- lapply(datgwr, function(d) { rownames(d) <- 1:nrow(d); coordinates(d) <- c("x", "y"); return(d) })

dat    <- list(datlm, datgwr)
dat.gw <- list(datlm.gw, datgwr.gw)

model <- formula(response ~ x1 + x2 + x3)
kernel <- "gaussian"

set.seed(2025)

# data for LM or GWR
for (i in 1:2) {
    # for each subdata
    for (j in 1:length(datlm)) {
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
        AICc.gwr <- gwgwr$GW.diagnostic$AICc 
        AICc.lm  <- AICc(modellm)
        print(AICc.lm)
        #aic <- AIC(modellm)
        #n <- nobs(modellm); k <- length(coef(modellm))
        #aic <- aic + ((2 * k * (k+1)) / (n - k - 1))
        #print(aic)
        #stop()
        #AICc <- mygwr$AICs[1, 1]
        AICs <- c(AICc.lm, AICc.gwr); names(AICs) <- c("lm", "gwr")

        # -----------------------------  
        # BIC scores
        # -----------------------------  
        #bic <- bbmle::BIC(modellm)
        #print(bic)
        #stop()
        
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
        print(res)
    
        stop()
    }
}

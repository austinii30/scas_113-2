###########################################################
# Filename: test.R
# Purpose : Test justification approaches 
# Author  : Austin Liu
# Date    : 2025/06/15
###########################################################

library(GWmodel)
#devtools::install_version("MuMIn", version = "1.47.1", repos = "https://cran.r-project.org")
library(MuMIn)  # for AICc()
library(ape)    # for Moran.I()
library(parallel)

source(file.path(this.path::this.dir(), "utils.R"))
import("func.R")
import("myGWRfunc.R")

RNGkind("L'Ecuyer-CMRG")
seed <- 2025
set.seed(seed)
c1 <- makeCluster(7)
clusterSetRNGStream(c1, seed)

sigmas <- c("0.1", "0.5", "1", "1.5", "2", "2.5", "3")
for (s in sigmas) {
    cat("Sigma: ", s, "    \n", sep="")
    
    load(datpath(paste0("datlm-sigma-", s, ".RData")))
    datlm <- lapply(datlm, data.frame) 
    datlm.gw <- lapply(datlm, function(d) { coordinates(d) <- c("x", "y"); return(d) }) 
    
    load(datpath(paste0("datgwr-sigma-", s, ".RData")))
    datgwr <- lapply(datgwr, data.frame)
    datgwr.gw <- lapply(datgwr, function(d) { coordinates(d) <- c("x", "y"); return(d) })
    
    dat    <- list(datlm, datgwr)
    dat.gw <- list(datlm.gw, datgwr.gw)
    
    model <- formula(response ~ x1 + x2 + x3)
    kernel <- "gaussian"
    datnames <- c("LM", "GWR")

    exetimes <- list()
    
    # data for LM or GWR
    for (i in 1:2) {
        cat("Begin simulation for ", datnames[i], ". \n", sep="")
        begtime <- Sys.time()
    
        clusterExport(c1, ls(), envir=environment()) 
        clusterExport(c1, ls(envir=.GlobalEnv), envir=.GlobalEnv) 
    
        # for each subdata
        #for (j in 1:length(datlm)) {
        #simres <- lapply(1:length(datlm), function(j) { 
        simres <- parLapply(c1, 1:length(datlm), function(j) { 
            #sink(outpath("junk"))
            # data for my code and GWmodel
            dlm <- dat[[ i ]][[ j ]]
            dgw <- dat.gw[[ i ]][[ j ]]
            library(GWmodel)
            library(MuMIn)  # for AICc()
            library(ape)    # for moran()
    
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
    
            #gwgwr <- gwr.basic(model, dgw, bw=bwAICc, kernel=kernel)
            #print(gwgwr)
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
            # Moran's I
            # -----------------------------  
            morani <- unlist(Moran.I(modellm$residuals, invDist(distanceMtx(dlm[, 2:3]))))
            
            # -----------------------------  
            # Monte Carlo test (p-values)
            # -----------------------------  
            mcAICc <- gwr.montecarlo(model, dgw, 199, kernel, bw=bwAICc)
    
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
            #sink()
            #cat(datnames[i], ", data ", idx, " done.    \r", sep="")
            return(res)
        })
        #cat("\n")
    
        endtime <- Sys.time()
        exetimes[[i]] <- duration(endtime - begtime)
            
        save(simres, file=datpath(paste0("simRes-", datnames[i], "-", length(datlm), "-sigma-", s, ".RData")))
        cat("Done simulating ", datnames[i], ", time spent: ", exetimes[[i]],
            " seconds. \n\n", sep="")
    }
    save(exetimes, file=datpath(paste0("exetime-sigma-", s, ".RData")))
}

stopCluster(c1)

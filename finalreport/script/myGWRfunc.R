##################################################
# Filename: myFunc.R
# Purpose : Hand written functions for GWR (WP1)
# Author  : Liu, Chih-Tse
# Date    : 2025/04/26
##################################################

library(parallel)
library(matrixcalc)  # for 'isSingular()'


# calculate population variance
popVar <- function(x) {
    return( mean((x-mean(x))^2) ) 
}


# correted AIC
ms.AICc <- function (n, sigmaHat, S) {
    return( 2*n*log(sigmaHat) + n*log(2*pi) + n*(n+tr(S))/(n-2-tr(S)) )
}


# AIC
ms.AIC <- function (n, sigmaHat, S) {
    return( 2*n*log(sigmaHat) + n*log(2*pi) + n + tr(S) )
}


# trace function
tr <- function (X) {
    return(sum(diag(X)))
}


# check if matrix is singular
isSingular <- function (X) {
    if (sum(is.na(X))!=0) return(TRUE)  # if NA exist
    if (!qr(X)$rank == ncol(X)) return(TRUE)
    if (min(abs(eigen(X)$values)) < 1e-10) return(TRUE)

    return(FALSE)
}


# Gaussian weighting function
wk.Gaussian <- function (dij, bandwidth, step=FALSE) {
    if (step) {
        if (dij > bandwidth) return( 0 )
    }
    return( exp(-(1/2) * (dij/bandwidth)^2) )
}


# alternative Gaussian weighting function
wk.Gaussian.1 <- function (dij, bandwidth, step=FALSE) {
    if (step) {
        if (dij > bandwidth) return( 0 )
    }
    return( exp(-(dij/bandwidth)^2) )
}


# bisquare weighting function
wk.biSquare <- function (dij, bandwidth) {
    if (dij > bandwidth) return( 0 )
    return( (1 - (dij/bandwidth)^2)^2 )
}


kernelFunction <- function (kernel) {
    "
    Select existing kernel functions (return error if 'kernel' does not exist).
    [Args]
        kernel (str): kernel function (Gaussian, biSquare, etc.)
    [Return]
        (function): selected kernel function
    "
    validKernels <- 
        list(Gaussian=wk.Gaussian, Gaussian.alt=wk.Gaussian.1, 
             biSquare=wk.biSquare, adaptKNN=wk.biSquare)
    kf <- validKernels[[kernel]]
    if (is.null(kf)) stop("Invalid kernel!")

    return (kf)
}


weighting <- function (distanceVec, kernel, bandwidth, adapt=FALSE) {
    "
    Calculate weighting for each observation.
    [Args]
        distanceVec (numeric): vector of distances
        kernel (str): weighting kernel
        bandwidth (float): pre-defined bandwidth
        adapt (bool): adaptive kernel or not
    [Return]
        (numeric): weighting for each observation
    "
    kf <- kernelFunction(kernel)
    if (adapt) bandwidth <- sort(distanceVec)[bandwidth]
    Wi <- sapply(distanceVec, kf, bandwidth=bandwidth)

    return(Wi)
}


distanceMtx <- function (coord, xIdx=1, yIdx=2) {
    "
    Calculate (Euclidean) distance matrix for a matrix of coordinates.
    [Args]
        coord (matrix): data (column: covariates, row: observations)
        xIdx (int): index of the X coordinates
        yIdx (int): index of the Y coordinates
    [Return]
        (matrix): distance matrix
    "
    x <- coord[, xIdx] 
    y <- coord[, yIdx] 
    n <- nrow(coord)

    D <- matrix(0, nrow=n, ncol=n)
    for (i in 1:(n-1)) 
        for (j in (i+1):n)
            D[i, j] <- sqrt( (x[i]-x[j])^2 + (y[i]-y[j])^2 )

    return( D + t(D) )
}

goldenSectionSearch <- function (a0, b0, scoreFunc, mStr, eps=1, discrete=FALSE) {
    # use golden section search to find the local minimum
    #eps <- 1  # epsilon, the threshold for stopping the algorithm
    r <- (sqrt(5)-1)/2  # the golden ratio

    if (discrete) { a0 <- floor(a0); b0 <- ceiling(b0) }
    x1 <- r*a0 + (1-r)*b0
    x2 <- (1-r)*a0 + r*b0
    if (discrete) { x1 <- floor(x1); x2 <- ceiling(x2) }
    fx1 <- scoreFunc(x1)
    fx2 <- scoreFunc(x2)

    #if (discrete) { xvec <- rep(NA, length(a0:b0)) }
    #xvec[1, 2] <- c(a0, b0)
    lasta0 <- a0; lastb0 <- b0
    
    nIter <- 0
    while ((b0-a0) > eps) {
        nIter <- nIter + 1

        # compare f(x1) and f(x2), determine which interval to zoom in
        if (fx1 < fx2) {
            b0 <- x2
            # old x1 becomes new x2
            x2 <- x1; fx2 <- fx1
            x1 <- r*a0 + (1-r)*b0
            if (discrete) x1 <- floor(x1)
            fx1 <- scoreFunc(x1)
        } else {
            a0 <- x1
            # old x2 becomes new x1
            x1 <- x2; fx1 <- fx2
            x2 <- (1-r)*a0 + r*b0
            if (discrete) x2 <- ceiling(x2)
            fx2 <- scoreFunc(x2)
        }
        cat(mStr, ", Iteration: ", nIter, ", Interval: [ ", a0, " , ", b0, " ]\n", sep="")

        if (a0 == lasta0 & b0 == lastb0) {
            xs <- a0:b0
            fxs <- sapply(xs, scoreFunc)
            # assign fx1 as the minimum
            fx1 <- sort(fxs)[1]
            x1  <- xs[order(fxs)][1]
            fx2 <- sort(fxs)[length(fxs)]
            x2  <- xs[order(fxs)][length(fxs)]
            break
        } else {
            lasta0 <- a0; lastb0 <- b0
        }
    }
    
    # return the (local) optimal bandwidth
    # NOTE: why choose a0 or b0 in the original code?
    x <- ifelse(fx1 <= fx2, x1, x2)
    s <- ifelse(fx1 <= fx2, fx1, fx2)
    cat("Optimal: ", x, 
        ", ", mStr, ": ", s, "\n", sep="")

    res <- list(nIter=nIter, result=x, score=s)
    return (res)
}


GWR <- function (formula, data, coordIdx, bandwidth, 
                 kernel="Gaussian", adapt=FALSE, 
                 hatMatrix=FALSE, seCompare=FALSE, AICs=FALSE,
                 doLRT=FALSE, alphaCI=NULL, nCluster=NULL) {
    "
    Basic geographically weighted regression.
    [Args]  
        formula (formula): model formula
        data (data.frame): the data frame
        coordIdx (character): columns of coordinates (eg. c('x', 'y'))
        bandwidth (float): bandwidth for the weighting kernel
        kernel (str): kernel function
        adapt (bool): adaptive bandwidth or fixed bandwidth
        hatMatrix (bool): to give information about the hat matrix or not
        seCompare (bool): to compare standard errors or not 
        AICs (bool): to calculate AICc and AIC scores or not 
        doLRT (bool): to do likelihood ratio test or not
        alphaCI (float): alpha for confidence interval of the coefficients
        nCluster (int): how many clusters to assign for parallel computing
    [Return]
        (list): list of objects
    "
    # check input arguments
    getCI <- ifelse(is.null(alphaCI), FALSE, TRUE)
    
    # data preprocessing
    #data <- data[complete.cases(data), ]
    n <- nrow(data)
    data$Wi <- rep(NA, n)

    # fit the global model
    model.g <- lm(formula=formula, data=data)
    beta.g   <- summary(model.g)$coefficients[, "Estimate"]
    betaSE.g <- summary(model.g)$coefficients[, "Std. Error"]
    sigma.g  <- summary(model.g)$sigma
    X.g <- model.matrix(formula, data)
    #H.g <- X.g %*% solve(t(X.g) %*% X.g) %*% t(X.g)
    H.g <- diag(hatvalues(model.g))
    nX <- length(beta.g)
    yName <- colnames(model.g$model)[1]

    # output
    allBeta <- matrix(NA, nrow=n, ncol=nX) # include intercept
    betaSummary <- matrix(NA, nrow=nX, ncol=6)
    rownames(betaSummary) <- names(beta.g)
    colnames(betaSummary) <-
        c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.", "Global")
    
    # calculate distance matrix
    coord <- data[, coordIdx]
    D <- distanceMtx(coord)

    if (hatMatrix | AICs)
        S <- matrix(NA, nrow=n, ncol=n)

    # for each location, get 'W' and fit the WLS
    yihat <- rep(NA, n)
    if (getCI) allCI <- as.list(rep(NA, n))

    # function to fit gwr
    gwrfit <- function (i) {
        if (getCI) {
            CIi <- matrix(NA, nrow=nX, ncol=2)  # lower/upper bond
            rownames(CIi) <- names(beta.g)
            colnames(CIi) <- c("lower", "upper")
        }

        # calculate weights for all points
        data$Wi <- weighting(D[, i], kernel, bandwidth, adapt)

        # fit the WLS
        model.i <- lm(formula=formula, weights=Wi, data=data)
        beta.i  <- as.vector(summary(model.i)$coefficients[, "Estimate"])
        yih <- predict(model.i, data[i, ])

        res <- list(beta.i=beta.i, yihat=yih)

        if (hatMatrix | AICs | getCI) {
            X <- model.matrix(formula, data)  # design matrix
            W.i <- diag(data$Wi)
            Ci <- solve(t(X) %*% W.i %*% X) %*% t(X) %*% W.i
            Si <- X[i, ] %*% Ci
            res$Si <- Si
        }

        if (getCI) {
            varBetai <- diag(Ci %*% t(Ci) * sigma.g^2)
            CIi[, 1] <- beta.i - abs(qnorm(alphaCI/2))*sqrt(varBetai)
            CIi[, 2] <- beta.i + abs(qnorm(alphaCI/2))*sqrt(varBetai)
            res$CI <- CIi
        }

        #res <- list(beta.i=beta.i, CI=CIi, yihat=yih, Si=Si)
        return(res)
    }
    
    # fit gwr for each observation
    if (!is.null(nCluster)) {
        c1 <- makeCluster(nCluster)
        # default 'envir' for 'clusterExport()' is .GlobalEnv!!!
        clusterExport(c1, ls(envir=environment()), envir=environment())
        clusterExport(c1, ls(envir=.GlobalEnv), envir=.GlobalEnv)
        #print(clusterEvalQ(c1, ls()))  # list available objects in each node
        gwrres <- parLapply(c1, 1:n, gwrfit)
        stopCluster(c1)
    } else {
        gwrres <- lapply(1:n, gwrfit)
    }

    # extract each component 
    allBeta <- do.call(rbind, lapply(gwrres, `[[`, 1))
    yihat <- sapply(gwrres, `[[`, 2)
    if (hatMatrix | AICs | getCI) 
        { S <- do.call(rbind, lapply(gwrres, `[[`, 3)) } 
    if (getCI) 
        { allCI <- lapply(gwrres, `[[`, (3+(hatMatrix | AICs | getCI))) } 

    RSS <- sum((data[, yName] - yihat)^2)

    # construct beta summary
    betaSummary[, 1:5] <- t(apply(allBeta, 2, fivenum))
    betaSummary[, 6] <- as.vector(beta.g)

    # results to return
    res <- list(allBeta=allBeta, summary=betaSummary, residuals=(data[, yName]-yihat))

    if (hatMatrix | doLRT) {
        effP  <- tr(S)  # including the intercept
        effDF <- n - 2*tr(S) + tr(t(S) %*% S)
        res$hatMatrix <- c(RSS, effP, effDF)
        names(res$hatMatrix) <- 
            c("Residual SS", "Effective P", "Effective Residual df")
    }

    if (seCompare) {
        seCompare <- matrix(NA, nrow=length(beta.g), ncol=2)
        colnames(seCompare) <- c("Local SE", "Global SE")
        rownames(seCompare) <- names(beta.g)
        seCompare[, 1] <- as.numeric(apply(allBeta, 2, sd))
        seCompare[, 2] <- betaSE.g
        res$seCompare <- seCompare
    }

    if (AICs) {
        sigmaHat <- sqrt(RSS/n)
        AICc <- ms.AICc(n, sigmaHat, S)
        AIC  <- ms.AIC(n, sigmaHat, S)
        AICc.g <- ms.AICc(n, sigma.g, H.g)
        AIC.g  <- ms.AIC(n, sigma.g, H.g)
        AICs <- matrix(c(AICc, AIC, AICc.g, AIC.g), nrow=2, byrow=TRUE)
        colnames(AICs) <- c("AICc", "AIC")
        rownames(AICs) <- c("GWR", "Global")
        res$AICs <- AICs
    }

    if (getCI) { res$CIs <- allCI }

    if (doLRT) {
        df1 <- (n-nX); df2 <- effDF
        f0 <- ((sigma.g^2*df1) / df1) / (RSS / df2)
        pValue <- 1-pf(f0, df1=df1, df2=df2)
        LRTres <- c(f0, pValue, df1, df2)
        names(LRTres) <- c("f0", "p-value", "df1", "df2")
        res$LRT <- LRTres
    }

    return(res)
}


cvssBandwidth <- function (formula, data, coordIdx, 
                           kernel="Gaussian", adapt=FALSE, by=1) {
    "
    Bandwidth selection by cross-validation. Loop over all possible bandwidths and return the optimal one. If the hat matrix is singular, return NA.
    [Args]
        by (float): increament of the sqauence
    [Return]
        (float): optimal bandwidth
        (.csv): a matrix of kernel, bandwidth, and scores
    "
    # get the maximum distance between any 2 locations
    coord <- data[, coordIdx]
    n <- nrow(data)
    D <- distanceMtx(coord)
    minBW <- min(D[D > 0])
    maxBW <- max(D)
    if (adapt) { minBW <- 1; maxBW <- n }

    # calculate CV score for bandwidth 'bw'
    cvScore <- function (bw) {
        error2s <- sapply(1:n, FUN = function(i) {
            data.i <- data[-i, ] 

            # calculate weights for all points
            data.i$Wi <- weighting(D[, i][-i], kernel, bw, adapt)

            # fit the WLS
            model.i <- lm(formula=formula, weights=Wi, data=data.i)
            yName <- colnames(model.i$model)[1]

            # predict for location i, calculate squared error
            error2 <- (data[i, yName] - predict(model.i, data[i, ]))^2
            return(error2)
        })
        return( sum(error2s) )
    }

    gss <- goldenSectionSearch(minBW, maxBW, cvScore, "CV", discrete=adapt)
    bw <- gss$result
    return(bw)
}


gcvBandwidth <- function (formula, data, coordIdx, 
                          kernel="Gaussian", adapt=FALSE, by=1) {
    "
    Bandwidth selection by generalized cross-validation. Loop over all possible bandwidths and return the optimal one. If the hat matrix is singular, return NA.
    [Return]
        (float): optimal bandwidth
        (.csv): a matrix of kernel, bandwidth, and scores
    "
    # get the maximum distance between any 2 locations
    coord <- data[, coordIdx]
    D <- distanceMtx(coord)
    n <- nrow(data)
    minBW <- min(D[D > 0])
    maxBW <- max(D)
    if (adapt) { minBW <- 1; maxBW <- n }

    # gcv score for bandwidth 'bw'
    gcvScore <- function (bw) {
        # returns the hat matrix and an additional error column
        # NOTE: each column in 'S_' is a sapply() returned value
        S_ <- sapply(1:n, FUN = function (i) {
            # calculate weights for all points
            data$Wi <- weighting(D[, i], kernel, bandwidth=bw, adapt)
            data$Wi[data$Wi < 1e-10] <- 0

            # fit the WLS
            model.i <- lm(formula=formula, weights=Wi, data=data)
            yName <- colnames(model.i$model)[1]

            X <- model.matrix(formula, data)  # design matrix X
            W.i <- diag(data$Wi)

            # if the psudo inverse is singular, return NA
            psuInv <- t(X) %*% W.i %*% X
            #print(psuInv)
            if (isSingular(psuInv)) return(rep(NA, (n+1)))

            # get row i for the hat matrix
            Hi <- X[i, ] %*% solve(psuInv) %*% t(X) %*% W.i

            # predict for location i, calculate squared error
            error2 <- (data[i, yName] - predict(model.i, data[i, ]))^2
            return(c(Hi, error2))
        })
        S_ <- t(S_) 
        S <- S_[, 1:n]; error2s <- S_[, (n+1)]
        RSS <- sum(error2s); sigmaHat <- sqrt(RSS/n)
        v1 <- tr(S)

        return(sum(error2s) * (n / (n-v1)^2))
    }

    gss <- goldenSectionSearch(minBW, maxBW, gcvScore, "GCV", discrete=adapt)
    bw <- gss$result
    return(bw)
}


aicBandwidth <- function (formula, data, coordIdx, 
                          kernel="Gaussian", adapt=FALSE,
                          by=1, aicType="c") {
    "
    Bandwidth selection by corrected AIC or AIC. Loop over all possible bandwidths and return the optimal one. If the hat matrix is singular, return NA.
    [Return]
        (float): optimal bandwidth
        (.csv): a matrix of kernel, bandwidth, and scores
    "
    # if 'aicType' is invalidly specified, make it type 'c'
    aicType <- ifelse(aicType %in% c("c", "nc"), aicType, "c")

    # get the maximum distance between any 2 locations
    coord <- data[, coordIdx]
    D <- distanceMtx(coord)
    n <- nrow(data)
    minBW <- min(D[D > 0])
    maxBW <- max(D)
    if (adapt) { minBW <- 1; maxBW <- n }

    # AICc score for bandwidth 'bw'
    aicScore <- function (bw) {
        # returns the hat matrix and an additional error column
        # NOTE: each column in 'S_' is a sapply() returned value
        S_ <- sapply(1:n, FUN = function (i) {
            # calculate weights for all points
            data$Wi <- weighting(D[, i], kernel, bandwidth=bw, adapt)
            data$Wi[data$Wi < 1e-10] <- 0

            # fit the WLS
            model.i <- lm(formula=formula, weights=Wi, data=data)
            yName <- colnames(model.i$model)[1]

            X <- model.matrix(formula, data)  # design matrix X
            W.i <- diag(data$Wi)

            # if the psudo inverse is singular, return NA
            psuInv <- t(X) %*% W.i %*% X
            if (isSingular(psuInv)) return(rep(NA, (n+1)))

            # get row i for the hat matrix
            Hi <- X[i, ] %*% solve(psuInv) %*% t(X) %*% W.i

            # predict for location i, calculate squared error
            error2 <- (data[i, yName] - predict(model.i, data[i, ]))^2
            return(c(Hi, error2))
        })
        S_ <- t(S_) 
        S <- S_[, 1:n]; error2s <- S_[, (n+1)]
        RSS <- sum(error2s); sigmaHat <- sqrt(RSS/n)

        if (aicType == "c") AIC <- ms.AICc(n, sigmaHat, S)
        else                AIC <- ms.AIC(n, sigmaHat, S)

        return(AIC)
    }

    m <- paste0("AIC", aicType)
    gss <- goldenSectionSearch(minBW, maxBW, aicScore, m, discrete=adapt)
    bw <- gss$result
    return(bw)
}


mcTest <- function (nIter, data, coordIdx, formula, 
                    silent=FALSE, seed=123, coreNum=NULL, ...) {
    "
    Monte Carlo based tests for significance of the local coefficients.
    [Args]
        nIter (int): number of iterations
        seed (float): random seed
        coreNUm (int): amount of core to utilize 
        ... : other arguments passed to GWR()
    [Return]
        (numeric): p-Values for each coefficient
    "
    if (is.null(coreNum)) 
        coreNum <- ifelse(detectCores()-5 > 1, detectCores()-5, 1)
    c1 <- makeCluster(coreNum)  # specify number of nodes
    # 'clusterExport' export vars. from the global environment by default
    # to export vars. from a function, must specify the environment 
    # 'environment()' gets the current environment
    # ensure reproducibility when using parallel computing
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    c1 <- makeCluster(coreNum)  # specify number of nodes
    clusterSetRNGStream(c1, seed)
    clusterExport(c1, ls(), envir=environment()) 
    clusterExport(c1, ls(envir=.GlobalEnv), envir=.GlobalEnv) 
    
    Vs <- parLapply(c1, 1:nIter, fun = function(i_, ...) {
        library(parallel); library(matrixcalc)  # for 'isSingular()'
    #Vs <- lapply(1:nIter, FUN = function(i_, ...) {
        cat("Conducting Monte Carlo test, iteration: ", i_, 
            "\r", sep="")
        shuf <- sample(1:nrow(data))
        data.i <- data
        data.i[, coordIdx] <- data.i[shuf, coordIdx]
        
        betaMtx <- GWR(formula, data.i, coordIdx, ...)$allBeta
        vs <- apply(betaMtx, 2, popVar)
        # 'vs' is a vector of SEs for all coefficients
        return(vs) }, 
        ...)
    stopCluster(c1)  # stop the cluster, free up space (important!!!)
    cat("\n")
    Vs <- t(do.call(cbind, Vs))  # turn a list of vectors into a matrix
    #print(Vs)
    #Vs <- t(Vs)  # each column has variances for a coefficient
    # variance of coefficient under correct location
    gBeta <- GWR(formula, data, coordIdx, ...)$allBeta
    Vks <- apply(gBeta, 2, popVar)

    pValues <- sapply(1:length(Vks), FUN = function(i) {
       return(mean(Vs[, i] >= Vks[i]))
    })

    #print(pValues)

    # extract RHS variables
    names(pValues) <- c("(Intercept)", all.vars(formula)[-1])  

    if (silent == FALSE) {
        cat("Monte Carlo p-Values\n", sep="")
        print(pValues)
    }

    return(pValues)
}

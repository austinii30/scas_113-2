
source("../../myPKG/R/figure.R")
source("code.R")
source("../../myPKG/R/EDA.R")
source("../../myPKG/R/dfManipulation.R")


# 1.(a)
set.seed(2025)

k <- 1000
#n <- 1000
#n <- 100
#n <- 200
n <- 500
theta <- c(0, 0.05, 0.1, 0.15, 0.2)

simDat <- list()
# for each theta, draw k random samples, each with n=100
for (th in theta) {
    dat <- c()
    for (i in 1:k)
        dat <- c(dat, arima.sim(n=n, list(ar=c(th, th))))
    # a column is a sample (k columns in total)
    simDat[[length(simDat)+1]] <- matrix(dat, ncol=k)
}

# calculate correlation
simDatCor <- list()
for (dat in simDat) {
    # returns a matirx with 2 rows, each is a correlation
    simCor <- sapply(1:ncol(dat), FUN = function(sIdx) {
            sample <- dat[, sIdx]
            cor1 <- cor(sample[-1], sample[-n])
            cor2 <- cor(sample[-c(1, 2)], sample[-c(n-1,n)])
            return(c(cor1, cor2))
        }) 
    # transpose so row: a sample, col: cor1 & cor2
    simDatCor[[length(simDatCor)+1]] <- t(simCor)
}
names(simDatCor) <- c("0", "0.05", "0.1", "0.15", "0.2")

HDB <- function(dat, varName, ct="Set1") {
    "
    Histogram, density, and box plots for numeric and integer variables.
    "
    #devtools::install_github("psyteachr/introdataviz")
    #install.packages("Hmisc")
    cat("Drawing univariate raincloud plots...\n")

    if (! class(dat) %in% c("integer", "numeric"))
        stop("'dat' can only be 'integer' or 'numeric'!")

    df <- data.frame(dat)
    colnames(df) <- c(varName)
    
    #colortheme <- "Set1"
    #colortheme <- "Dark2"
    colortheme <- ct
    legendPos <- "none"
    titleText <- 16
    axisText  <- 21
    rain_height <- 0.1
    
    # add a variable 'dummy' which all rows have equal value
    df$dummy <- factor("same")

    # draw each plot separately
    for (varName in colnames(df)) {
        if (varName == "dummy") next

        dynamicCode <- "
var <- df$`__var__`
varName <- \"__var__\"

fig <- 
ggplot(df, aes(x = \"\", y = __var__, fill = dummy)) +

# clouds
# 'alpha' controls the transparancy of the graph (1: opague, 0: transparent)
introdataviz::geom_flat_violin(trim=TRUE, alpha = 0.4,
                               position = position_nudge(x = rain_height+.05)) +
# rain
geom_point(aes(colour = dummy), size = 2, alpha = .5, show.legend = FALSE, 
           position = position_jitter(width = rain_height, height = 0)) +
# boxplots
geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
             outlier.shape = NA,
             position = position_nudge(x = -rain_height*2)) +
# mean and SE point in the cloud
#stat_summary(fun.data = mean_cl_normal, mapping = aes(color = dummy), show.legend = FALSE,
stat_summary(fun = \"mean\", mapping = aes(color = dummy), show.legend = FALSE,
             position = position_nudge(x = rain_height * 3)) +
# adjust layout
#ggtitle(varName) + 
    scale_x_discrete(name = \"\", expand = c(rain_height*3, 0, 0, 0.7)) +
    scale_y_continuous(name = \"\", limits = c(-0.5, 0.5)) + 
    coord_flip() +
    # custom colours and theme
    scale_fill_brewer(palette = colortheme, name = \"Location\") +
    scale_colour_brewer(palette = colortheme) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = legendPos,
          legend.background = element_rect(fill = \"white\", color = \"white\"),
          text = element_text(size = 14), 
          plot.title = element_text(size = titleText,  hjust = 0.5),
          axis.text = element_text(size = axisText))
"
        dynamicCode <- gsub("__var__", varName, dynamicCode)
        return(eval(parse(text=dynamicCode)))
    }
}


# calculate summary statistics
d <- simDatCor[[1]]
for (i in 2:length(simDatCor))
    d <- cbind(d, simDatCor[[i]])
colnames(d) <- c("0_1", "0_2", "0.05_1", "0.05_2", 
                 "0.1_1", "0.1_2", "0.15_1", "0.15_2", 
                 "0.2_1", "0.2_2")
d <- data.frame(d)
ss <- summaryStatistics(d)
for (i in 1:ncol(ss)) ss[, i] <- round(ss[, i], 3)
write.csv(ss, paste0("1(a)_summaryStatistics_n", n, "_k", k, ".csv"))


# draw HBDs 
allHBD <- list()
for (i in 1:ncol(d)) {
    dat <- d[, i]
    varName <- colnames(d)[i]

    if (i %% 2 == 0) p <- HDB(dat, varName, "Dark2")
    else             p <- HDB(dat, varName)
    allHBD[[length(allHBD)+1]] <- p

    #expPath <- paste0(varName, ".pdf")
    #pdf(expPath, width = 9, height = 6)
    #print(p)
    #dev.off()
}
pdf(paste0("1(a)_corHBD_n", n, "_k", k, ".pdf"), width = 18, height = 25)
grid.arrange(grobs = allHBD, ncol = 2, nrow = 5)
dev.off()


warnings()



# 1.(b)
a <- 0.05
pRes <- matrix(nrow=3, ncol=0)
pNAs <- matrix(nrow=3, ncol=0)
for (i in 1:length(simDat)) {
    dat <- simDat[[i]]
    p.gap <- apply(dat, 2, FUN = gapTest)
    p.uad <- apply(dat, 2, FUN = uadTest)
    p.pmt <- apply(dat, 2, FUN = pmtTest)

    pRes <- cbind(pRes, c(sum(p.gap<a, na.rm=TRUE), 
                          sum(p.uad<a, na.rm=TRUE), 
                          sum(p.pmt<a, na.rm=TRUE)))
    pNAs <- cbind(pNAs, c(sum(is.na(p.gap)), sum(is.na(p.uad)), sum(is.na(p.pmt))))
}
colnames(pRes) <- colnames(pNAs) <- names(simDatCor)
rownames(pRes) <- rownames(pNAs) <- c("Gap", "UaD", "Pmt")

write.csv(pRes, paste0("1(b)_rejH0Times_n", n, "_k", k, "_a", a, ".csv"))
write.csv(pNAs, paste0("1(b)_cannotTest_n", n, "_k", k, "_a", a, ".csv"))


warnings()



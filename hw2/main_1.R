
source("../../myPKG/R/figure.R")
source("code.R")
source("../../myPKG/R/EDA.R")
source("../../myPKG/R/dfManipulation.R")


# 1.(a)
set.seed(2025)

k <- 1000
#k <- 100
n <- 100
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



HDB <- function(dat, varName, expPath="") {
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
    
    colortheme <- "Set1"
    legendPos <- "none"
    titleText <- 16
    axisText  <- 11
    rain_height <- 0.1
    
    # add a variable 'dummy' which all rows have equal value
    df$dummy <- factor("same")

    # draw each plot separately
    for (varName in colnames(df)) {
        if (varName == "dummy") next

        dynamicCode <- "
if (expPath != \"\")
    pdf(expPath, width = 8, height = 6)

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
ggtitle(varName) + 
    scale_x_discrete(name = \"\", expand = c(rain_height*3, 0, 0, 0.7)) +
    scale_y_continuous(name = \"\") + 
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

if (expPath != \"\")
    dev.off()
"
        dynamicCode <- gsub("__var__", varName, dynamicCode)
        eval(parse(text=dynamicCode))
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
write.table(ss, "ss_1.csv", fileEncoding="UTF-8")

stop()

# draw HBDs 
for (i in 1:length(simDatCor)) {
    for (j in 1:2) {
        dat <- simDatCor[[i]][, j]
        varName <- names(simDatCor)[i]
        #varName <- expression(paste(theta,"=",names(simDatCor)[i]))
        HDB(dat, varName, paste0(varName, "_", j, ".pdf"))
    }
}


# 1.(b)





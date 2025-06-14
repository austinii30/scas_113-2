###########################################################
# Filename: q2.R
# Course  : 113-2 Statistical Computation and Simulation, HW5, 
#           Question 2 
# Author  : Potumas Liu 
# Date    : 2025/05/26
###########################################################

library(MASS)         # for kde2d()
library(fields)       # For image.plot (better legends)
library(RColorBrewer) # For color palette
library(grDevices)    # For colorRampPalette
library(ks)
library(np)
#library(ggtern)

setwd(this.path::this.dir())
dat <- read.csv("data-q2.csv")

# missing values (no missing values)
print(nrow(dat[!complete.cases(dat), ]))

dat88 <- dat[dat$year == 88, ]
dat89 <- dat[dat$year == 89, ]
write.csv(dat88, "dat88-q2.csv")
write.csv(dat89, "dat89-q2.csv")

jet.colors <- colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
      "#FF7F00", "red", "#7F0000"))

# maximum and minimum of longitude and latitude across 1999~2000
#print(min(dat$longitude))
#print(max(dat$longitude))
#print(min(dat$latitude))
#print(max(dat$latitude))
minlong <- 120; maxlong <- 124
minlat <- 22; maxlat <- 25

min(dat$depth)
max(dat$depth)
min(dat$intensity)
max(dat$intensity)

xlimcoord <- c(120, 124)
ylimcoord <- c(22, 25)
ilimcoord <- c(4, 7)
dlimcoord <- c(0, 120)

xminkde <- c(120, 22)
xmaxkde <- c(124, 25)


kernel_smooth_grid <- function(xo1, xo2, y, xn1, xn2, h1, h2) {
  z_mat <- matrix(NA, nrow=length(xn1), ncol=length(xn2))

  for (i in seq_along(xn1)) {
    for (j in seq_along(xn2)) {
      dx1 <- dnorm((xn1[i] - xo1) / h1)
      dx2 <- dnorm((xn2[j] - xo2) / h2)
      w <- dx1 * dx2
      w <- w / sum(w)

      z_mat[i, j] <- sum(w * y)
    }
  }
  return(z_mat)
}

plotlldi <- function (plotname, x, y, z, 
                      grid.size=200, xlcd=xlimcoord, ylcd=ylimcoord, zlcd=zlimcoord) {
# Extract coordinates
dat <- data.frame(x = x, y = y, z = z)

# Step 3: Prediction grid
grid.x <- seq(minlong, maxlong, length.out=grid.size)
grid.y <- seq(minlat, maxlat, length.out=grid.size)
grid.coords <- expand.grid(grid.x, grid.y)
colnames(grid.coords) <- c("x", "y")

z.matrix <- kernel_smooth_grid(x, y, z, grid.x, grid.y, 0.2, 0.2)

dens <- list(
    x = grid.x,
    y = grid.y,
    z = z.matrix
)

# Plot with a square aspect ratio and a tall legend
pdf(plotname, width=8, height=6.5)
par(mar = c(4.5, 5.5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Longitude", ylab = "", main = "",
      las = 1, yaxt="n",
      cex.axis=1.7, cex.lab=1.9,
      xlim=xlcd, ylim=ylcd, zlim=zlcd)

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Longitude", ylab = "", main = "", 
           las = 1, yaxt="n",
           cex.axis=1.7, cex.lab=1.9,
           xlim=xlcd, ylim=ylcd, zlim=zlcd)

# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.7)
mtext("Latitude", side=2, line=4, cex=1.9) 

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()
}


# -----------------------------------------------------------------------------
# 1999 depth 
# -----------------------------------------------------------------------------
x <- dat88$longitude
y <- dat88$latitude
z <- dat88$depth
plotlldi("q2-1999-depth.pdf", x, y, z, zlcd=dlimcoord)

x <- dat88$longitude
y <- dat88$latitude
z <- dat88$intensity
plotlldi("q2-1999-intensity.pdf", x, y, z, zlcd=ilimcoord)

x <- dat89$longitude
y <- dat89$latitude
z <- dat89$depth
plotlldi("q2-2000-depth.pdf", x, y, z, zlcd=dlimcoord)

x <- dat89$longitude
y <- dat89$latitude
z <- dat89$intensity
plotlldi("q2-2000-intensity.pdf", x, y, z, zlcd=ilimcoord)

stop()

# Extract coordinates
x <- dat88$longitude
y <- dat88$latitude
z <- dat88$depth
dat <- data.frame(x = x, y = y, z = z)

# Step 1: Bandwidth selection via least-squares cross-validation (CV)
bw <- npregbw(
  formula = z ~ x + y,
  regtype = "lc",         # local constant (Nadaraya-Watson)
  bwmethod = "cv.ls",     # cross-validation for least squares
  ckertype = "gaussian",  # Gaussian kernel
  data = dat
)
print(bw)

# Step 2: Fit kernel regression model
model <- npreg(bws = bw, data = dat)

# Step 3: Prediction grid
grid.size <- 200
grid.x <- seq(minlong, maxlong, length.out=grid.size)
grid.y <- seq(minlat, maxlat, length.out=grid.size)
grid.coords <- expand.grid(grid.x, grid.y)
class(grid.coords)
colnames(grid.coords) <- c("x", "y")

# Step 4: Predict smoothed z on grid
z_hat <- predict(model, newdata = grid.coords)
z.matrix <- matrix(z_hat, nrow = grid.size, byrow = FALSE)

## Combine into a matrix of coordinates
#coords <- cbind(x, y)

## Perform smoothing
#fit <- Tps(coords, z, lon.lat=TRUE)  # Thin Plate Spline smooth

#z.smooth <- predict(fit, grid.coords)

## Plot the smoothed surface
#z.matrix <- matrix(z.smooth, nrow = grid.size, ncol = grid.size)

# Compute 2D density
#res <- kde(x=cbind(x, y), xmin=xminkde, xmax=xmaxkde, w=z, gridsize=c(200, 200))
dens <- list(
    x = grid.x,
    y = grid.y,
    z = z.matrix
    #x = res$eval.points$x,
    #y = res$eval.points$y, 
    #z = res$estimate
)

# Plot with a square aspect ratio and a tall legend
pdf("q2-1999-depth.pdf", width=8, height=6.5)
par(mar = c(4.5, 5.5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Longitude", ylab = "", main = "",
      las = 1, yaxt="n",
      cex.axis=1.7, cex.lab=1.9,
      xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Longitude", ylab = "", main = "", 
           las = 1, yaxt="n",
           cex.axis=1.7, cex.lab=1.9,
           xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.7)
mtext("Latitude", side=2, line=4, cex=1.9) 

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()
stop()

# -----------------------------------------------------------------------------
# 1999 intensity
# -----------------------------------------------------------------------------
# Extract coordinates
x <- dat88$longitude
y <- dat88$latitude
z <- dat88$intensity

# Compute 2D density
res <- kde(x=cbind(x, y), xmin=xminkde, xmax=xmaxkde, w=z, gridsize=c(200, 200))
dens <- list(
    x = res$eval.points$x,
    y = res$eval.points$y, 
    z = res$estimate
)

# Plot with a square aspect ratio and a tall legend
pdf("q2-1999-intensity.pdf", width=8, height=6.5)
par(mar = c(4.5, 5.5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Longitude", ylab = "", main = "",
      las = 1, 
      cex.axis=1.7, cex.lab=1.9, yaxt="n",
      xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Longitude", ylab = "", main = "", 
           las = 1, yaxt="n",
           cex.axis=1.7, cex.lab=1.9,
           xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.7)
mtext("Latitude", side=2, line=4, cex=1.9) 

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()


# -----------------------------------------------------------------------------
# 2000 depth
# -----------------------------------------------------------------------------
# Extract coordinates
x <- dat89$longitude
y <- dat89$latitude
z <- dat89$depth

# Compute 2D density
res <- kde(x=cbind(x, y), xmin=xminkde, xmax=xmaxkde, w=z, gridsize=c(200, 200))
dens <- list(
    x = res$eval.points$x,
    y = res$eval.points$y, 
    z = res$estimate
)

# Plot with a square aspect ratio and a tall legend
pdf("q2-2000-depth.pdf", width=8, height=6.5)
par(mar = c(4.5, 5.5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Longitude", ylab = "", main = "",
      las = 1, yaxt="n", 
      cex.axis=1.7, cex.lab=1.9,
      xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Longitude", ylab = "", main = "", 
           las = 1, yaxt="n", 
           cex.axis=1.7, cex.lab=1.9,
           xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.7)
mtext("Latitude", side=2, line=4, cex=1.9) 

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()


# -----------------------------------------------------------------------------
# 2000 intensity
# -----------------------------------------------------------------------------
# Extract coordinates
x <- dat89$longitude
y <- dat89$latitude
z <- dat89$intensity

# Compute 2D density
res <- kde(x=cbind(x, y), xmin=xminkde, xmax=xmaxkde, w=z, gridsize=c(200, 200))
dens <- list(
    x = res$eval.points$x,
    y = res$eval.points$y, 
    z = res$estimate
)

# Plot with a square aspect ratio and a tall legend
pdf("q2-2000-intensity.pdf", width=8, height=6.5)
par(mar = c(4.5, 5.5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Longitude", ylab = "", main = "",
      las = 1, yaxt="n", 
      cex.axis=1.7, cex.lab=1.9,
      xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Longitude", ylab = "", main = "", 
           las = 1, yaxt="n", 
           cex.axis=1.7, cex.lab=1.9,
           xlim=xlimcoord, ylim=ylimcoord)#, zlim=zlimcoord)

# custom y-axis labels and axis
axis(2, las=1, cex.axis=1.7)
mtext("Latitude", side=2, line=4, cex=1.9) 

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()

#stop()

# -----------------------------------------------------------------------------
# 1999 intensity vs. depth
# -----------------------------------------------------------------------------
# Extract coordinates
x <- dat88$intensity
y <- dat88$depth

# Compute 2D density
dens <- kde2d(x, y, n = 200,
              lims=c(2.5, 7.5, 0, 120))  # Increase 'n' for better resolution

# Open a square plotting window (works best in RStudio or X11)
#asp_ratio <- 1
#plot_width <- 7
#plot_height <- 7
# For standalone script:
# windows(width = plot_width, height = plot_height)  # Windows
# quartz(width = plot_width, height = plot_height)   # macOS
# X11(width = plot_width, height = plot_height)      # Linux

# Plot with a square aspect ratio and a tall legend
pdf("q2-1999_intensity-depth.pdf", width=8, height=6.5)
par(mar = c(4.5, 5, 0.5, 6))  # Extra space on right for the legend

# draw the image without the legend
image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Intensity", ylab = "Depth", main = "",
      las = 1, 
      cex.axis=1.7, cex.lab=1.9,
      xlim=c(2.5, 7.5), ylim=c(0, 120), zlim=c(0, 0.045))

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Intensity", ylab = "Depth", main = "", 
           las = 1, 
           cex.axis=1.7, cex.lab=1.9,
           xlim=c(2.5, 7.5), ylim=c(0, 120), zlim=c(0, 0.045))

# Add contour lines if desired
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()


# -----------------------------------------------------------------------------
# 2000 intensity vs. depth
# -----------------------------------------------------------------------------
# Extract coordinates
x <- dat89$intensity
y <- dat89$depth

# Compute 2D density
dens <- kde2d(x, y, n = 200,
              lims=c(2.5, 7.5, 0, 120))  # Increase 'n' for better resolution
pdf("q2-2000_intensity-depth.pdf", width=8, height=6.5)
par(mar = c(4.5, 5, 0.5, 6))  # Extra space on right for the legend

image(dens$x, dens$y, dens$z,
      col = jet.colors(100),
      xlab = "Intensity", ylab = "Depth", main = "",
      las = 1, 
      cex.axis=1.7, cex.lab=1.9,
      xlim=c(2.5, 7.5), ylim=c(0, 120), zlim=c(0, 0.045))

image.plot(dens,
           legend.only=TRUE,
           col = jet.colors(100),
           legend.width = 2,        # Make legend longer
           legend.mar = 5,
           legend.shrink = 1.0,       # Don't shrink legend
           axis.args=list(cex.axis=1.7),  # change legend number size
           xlab = "Intensity", ylab = "Depth", main = "", 
           las = 1, 
           cex.axis=1.7, cex.lab=1.9,
           xlim=c(2.5, 7.5), ylim=c(0, 120), zlim=c(0, 0.045))

# contour lines
contour(dens, add = TRUE, drawlabels = FALSE, col = "black", lwd = 0.5)

dev.off()

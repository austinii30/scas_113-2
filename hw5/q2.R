###########################################################
# Filename: q2.R
# Course  : 113-2 Statistical Computation and Simulation, HW5, 
#           Question 2 
# Author  : Potumas Liu 
# Date    : 2025/05/26
###########################################################


library(MASS)
library(plotly)
library(lattice)
library(ggplot2)

setwd(this.path::this.dir())

dat <- read.csv("data-q2.csv")

str(dat)

# missing values (no missing values)
print(nrow(dat[!complete.cases(dat), ]))

dat88 <- dat[dat$year == 88, ]
dat89 <- dat[dat$year == 89, ]

write.csv(dat88, "dat88-q2.csv")
write.csv(dat89, "dat89-q2.csv")

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


pdf("test.pdf")

#ggplot(dat88, aes(x = longitude, y = latitude)) +
#  geom_density_2d(color="red") +
#  labs(
#    title = "2D Filled Density Plot of Old Faithful Geyser Data",
#    x = "Eruption Duration (minutes)",
#    y = "Waiting Time (minutes)"
#  ) +
#  theme_minimal()

#my_colors <- colorRampPalette(c("blue", "red"))(14)
#ggplot(dat88, aes(x = intensity, y = depth)) +
#  geom_density_2d_filled() +
#  scale_fill_manual(values = my_colors) +
#  labs(
#    title = "2D Filled Density Plot of Location Data",
#    x = "Longitude",
#    y = "Latitude",
#    fill = "Density Level"
#  ) +
#  theme_minimal()



#ggplot(dat88, aes(x = longitude, y = latitude)) +
#  stat_density_2d(aes(color = after_stat(level)), geom = "density_2d") +
#  scale_color_gradient(low = "blue", high = "red") +
#  labs(
#    title = "2D Density Contour Plot of Location Data",
#    x = "Longitude",
#    y = "Latitude",
#    color = "Density Level"
#  ) +
#  theme_minimal()

#ggplot(dat88, aes(x = longitude, y = latitude)) +
#  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#  scale_fill_gradientn(
#    colors = c("blue", "cyan", "green", "yellow", "orange", "red"),
#    name = "Density Level"
#  ) +
#  labs(
#    title = "2D Density Plot with Jet Color Scale",
#    x = "Longitude",
#    y = "Latitude"
#  ) +
#  theme_minimal()


#ggplot(dat88, aes(x = longitude, y = latitude)) +
#  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#  scale_fill_gradientn(
#    colours = c("blue", "cyan", "green", "yellow", "orange", "red"),
#    name = "Density"
#  ) +
#  labs(
#    title = "Heatmap-Style 2D Density",
#    x = "Longitude",
#    y = "Latitude"
#  ) +
#  theme_minimal() +
#  coord_fixed()  # keep aspect ratio square if appropriate

#ggplot(dat88, aes(x = longitude, y = latitude, fill = depth)) +
#  geom_raster() +
#  scale_fill_gradientn(
#    colours = c("blue", "cyan", "green", "yellow", "orange", "red"),
#    name = "Depth"
#  ) +
#  labs(
#    title = "Depth Heatmap",
#    x = "Longitude",
#    y = "Latitude"
#  ) +
#  theme_minimal() +
#  coord_fixed()  # keep square aspect ratio

#ggplot(dat88, aes(x = longitude, y = latitude)) +
#  geom_tile(aes(fill = depth)) +
#  scale_fill_gradientn(colors = jet.colors(100)) +
#  coord_fixed() +
#  theme_minimal()



#levelplot(depth ~ longitude * latitude, data = dat88,
#          cuts = 100, col.regions = jet.colors(100),
#          xlab = "Longitude", ylab = "Latitude")


dens <- kde2d(dat88$intensity, dat88$depth, n = 100)

plot_ly(
  x = dens$x, y = dens$y, z = dens$z,
  type = "heatmap",
  colors = jet.colors(100)
)



dev.off()




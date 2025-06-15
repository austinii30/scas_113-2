library(sp)
library(spgwr)
library(spData)
library(GWmodel)

source("./script/myGWRfunc.R")

data(meuse)

a <- meuse
coordinates(a) <- c("x", "y")
str(a)


formula <- formula(cadmium ~ dist)

ab <- gwr.basic(formula, a, bw=228, kernel="gaussian")
print(ab)

bw.gwr(formula, data=a, approach="CV", kernel="gaussian")
aa <- gwr.cv(bw=228, X=cbind(1, a$dist), Y=a$copper, kernel="gaussian", dp.locat=a@coords, dMat=gw.dist(a@coords))
print(aa/33.66762)

print(GWR(formula, meuse, c("x", "y"), 228, AICs=TRUE))



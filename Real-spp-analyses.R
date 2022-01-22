##Load climate data and run PCA
library(spatstat)
library(raster)

clim <- stack(paste0("Real-spp/Climate/bio", 1:19, ".tif"))
clim.points <- na.omit(data.frame(rasterToPoints(clim)))
clim.vars <- clim.points[, 3:ncol(clim.points)]

clim.pca <- princomp(clim.vars)
clim.pca.vals <- predict(clim.pca)
clim.pca.points <- data.frame(clim.points[, 1:2], clim.pca.vals)

clim.pca.r <- rasterFromXYZ(clim.pca.points)

dir.create("Real-spp/Climate/PCA")
for(i in 1:3)writeRaster(clim.pca.r[[i]], paste0("Real-spp/Climate/PCA/Comp-", i), "GTiff")

## Load spp data and format

spp <- lapply(paste0("Real-spp/Presence/", c("cal_cal_", "cal_mel_"), "train.csv"), read.csv)

#spp.lay.df <- lapply(spp.layers, function(x){data.frame(rasterToPoints(x))}) #Transform layers to dataframes
#spp.point.vals <- lapply(seq_along(spp.layers), function(x){extract(spp.layers[[x]], spp.points[[x]])}) #Extracting values at points


## Coercing raster package data to spatstat
ux = sort(unique(clim.pca.points$x)) #Extracting unique coordinates
uy = sort(unique(clim.pca.points$y))
nx = length(ux) #length of unique coordinates
ny = length(uy)
ref.cols = match(clim.pca.points$x, ux) #position of every data point
ref.lines = match(clim.pca.points$y, uy)
vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
vec[ref.vec] = 1
data.mask = matrix(vec, max(ref.lines), max(ref.cols), dimnames = list(uy, ux))
win = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window

## transforming species points to a planar point pattern
spp.ppp <- lapply(spp, function(x){
      ppp(x[, 'Latitud'], x[, 'Longitud'], window = win, check = F)
})

spp.lay <- list(clim.pca.points[, c("x", "y", "Comp.1")],
                clim.pca.points[, c("x", "y", "Comp.2")],
                clim.pca.points[, c("x", "y", "Comp.3")])

## Transforming species' layers to spatstat images
library(foreach)
lay.im.list <- foreach(i = 3:ncol(clim.pca.points)) %do% {
            vec.all = rep(NA, max(ref.lines)*max(ref.cols))
            vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
            vec.all[ref.vec] = clim.pca.points[,i]
            lay <- im(matrix(vec.all, max(ref.lines), max(ref.cols),
                             dimnames = list(uy, ux)), xcol = ux, yrow = uy)
            return(lay)
      }
names(lay.im.list) <- c("a", "b", "c")

##Analisis exploratorio de las respuestas para seleccionar variables bioclimáticas


## Ajuste de modelo
spp.ppms <- lapply(seq_along(spp.ppp), function(i){
      mod <-ppm(spp.ppp[[i]],
                trend = ~ a + b + c +
                      I(a^2) + I(b^2) + I(c^2),
                covariates = lay.im.list)
      s.mod <- step(mod)
      coef <- coefficients(mod)
      rast <- raster(predict(mod, type = "trend", ngrid = c(209, 259)))
      return(list(mod = mod,
                  s.mod = s.mod,
                  pred = rast,
                  coef = coef))
})


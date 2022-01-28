##Load climate data and run PCA
library(spatstat)
library(raster)

clim <- stack(paste0("Real-spp/Climate/bio", 1:19, ".tif"))
clim <- projectRaster(clim, crs = CRS("+init=epsg:4269"))

clim.points <- na.omit(data.frame(rasterToPoints(clim)))
clim.vars <- clim.points[, 3:ncol(clim.points)]

clim.pca <- princomp(clim.vars)
clim.pca.vals <- predict(clim.pca)
clim.pca.points <- data.frame(clim.points[, 1:2], clim.pca.vals)

clim.pca.r <- rasterFromXYZ(clim.pca.points[, 1:5])
names(clim.pca.r) <-  c("Comp.1", "Comp.2", "Comp.3")
proj4string(clim.pca.r) <- CRS("+init=epsg:4269")

dir.create("Real-spp/Climate/PCA")
for(i in 1:3)writeRaster(clim.pca.r[[i]], paste0("Real-spp/Climate/PCA/Comp-", i, "-NAD83"), "GTiff")

## Load spp data and format
spp <- lapply(paste0("Real-spp/Presence/", c("cal_cal_", "cal_mel_"), "train.csv"), function(x){
   require(rgdal)
   x1 <- read.csv(x)[, 2:3]
   coordinates(x1) <- c("Longitud", "Latitud")
   proj4string(x1) <- CRS("+init=epsg:4326")
   x2 <- spTransform(x1, CRS = CRS("+init=epsg:4269"))
   x3 <- data.frame(coordinates(x2))
   return(x3)
})
names(spp) <- c("Cal-cal", "Cal-mel")

spp.buffers <- lapply(spp, function(x){
   require(rgdal); require(rgeos)
   coordinates(x) <- c("Longitud", "Latitud")
   buf <- rgeos::gBuffer(x, width = 5, byid = F)
   return(buf)
})

## Cropping climate
clim.spp <- lapply(spp.buffers, function(x1){
   m <- mask(x = clim, mask = x1)
   m.df <- na.omit(rasterToPoints(m))
   m <- rasterFromXYZ(m.df)
   return(m)
})

clim.spp.points <- lapply(clim.spp, function(x){
   data.frame(rasterToPoints(x))
})

pca.spp <- lapply(spp.buffers, function(x1){
   m <- mask(x = clim.pca.r, mask = x1)
   m.df <- na.omit(rasterToPoints(m))
   m <- rasterFromXYZ(m.df)
   return(m)
})

##Working windows
win <- list(
   winFromRaster(raster(clim.spp[[1]][[1]])),
   winFromRaster(raster(clim.spp[[2]][[1]]))
) 

## transforming species points to a planar point pattern
spp.ppp <- lapply(1:2, function(x){
      ppp(spp[[x]][, 'Longitud'], spp[[x]][, 'Latitud'], window = win[[x]], check = F)
})

## Transforming species' layers to spatstat images
pca.im.list <- lapply(pca.spp, imFromStack)
bio.im.list <- lapply(clim.spp, imFromStack)

##Computing quadrat counts to smooth
quads <- list(
   ppp(clim.spp.points[[1]]$x, clim.spp.points[[1]]$y, window = win[[1]]),
   ppp(clim.spp.points[[2]]$x, clim.spp.points[[2]]$y, window = win[[2]])
   )

Q <- list(
   quadscheme(data = spp.ppp[[1]], dummy = quads[[1]], method = "grid",
              ntile = c(116, 168), npix = c(116, 168)),
   quadscheme(data = spp.ppp[[2]], dummy = quads[[2]], method = "grid",
              ntile = c(178, 194), npix = c(178, 194))
   )

####Exploratory analysis

##Analyses of dependence on bioclimatic variables

for(j in 1:2){
   Z.bio <- bio.im.list[[j]]
   
   cuts <- seq(0, 1, by = 0.1)
   
   quants <- lapply(Z.bio, function(x){quantile(x, probs = cuts, labels = 1:(length(cuts)-1))})
   
   bio.cut <- foreach(k = seq_along(Z.bio)) %do% {
      cut(Z.bio[[k]], breaks = quants[[k]], labels = 1:(length(cuts)-1))
   }
   
   V <- lapply(bio.cut, function(x)tess(image = x))
   
   counts <- foreach(k = seq_along(V)) %do% {
         quadratcount(Q[[j]]$data, tess = V[[k]])
   }
   
   pdf(paste0("Real-spp/",c("Cal-cal.pdf", "Cal-mel.pdf")[j]), width = 15, height = 5)
   par(mfrow = c(1,3))
   for(i in 1:19){
      plot(counts[[i]], main = paste(c("Callipepla californica",
                                           "Calamospiza melanocorys")[j], names(clim)[i], sep = ", "))
      plot(V[[i]], main = "")
      points(spp.ppp[[j]], pch = "+", col = "green", cex = 2)
      plot(rhohat(Q[[j]]$data, Z.bio[[i]]) , main = "")
   }
   dev.off()
}

png("Real-spp/Pairs-Cal-cal.png", width = 2000, height = 2000)
pairs(clim.spp[[1]])
dev.off()

png("Real-spp/Pairs-Cal-mel.png", width = 2000, height = 2000)
pairs(clim.spp[[2]])
dev.off()

cal.cal.formulas <- c("~ bio1 + bio8 + bio19 + I(bio1^2) + I(bio8^2) + I(bio19^2)",
                      "~ bio1 + bio8 + bio12 + I(bio1^2) + I(bio8^2) + I(bio12^2)",
                      "~ bio5 + bio8 + bio12 + I(bio2^2) + I(bio8^2) + I(bio12^2)",
                      "~ bio8 + bio11 + bio12 + I(bio8^2) + I(bio11^2) + I(bio12^2)")

cal.mel.formulas <- c("~ bio4 + bio5 + bio10 + bio16 + I(bio4^2) + I(bio5^2) + I(bio10^2) + I(bio16^2)",
                      "~ bio4 + bio5 + bio10 + I(bio4^2) + I(bio5^2) + I(bio10^2)",
                      "~ bio4 + bio5 + bio16 + I(bio4^2) + I(bio5^2) + I(bio16^2)",
                      "~ bio4 + bio10 + bio16 + I(bio4^2) + I(bio10^2) + I(bio16^2)",
                      "~ bio5 + bio10 + bio16 + I(bio5^2) + I(bio10^2) + I(bio16^2)",
                      "~ bio5 + bio7 + bio16 + I(bio5^2) + I(bio7^2) + I(bio16^2)",
                      "~ bio1 + bio4 + bio16 + I(bio1^2) + I(bio4^2) + I(bio16^2)",
                      "~ bio1 + bio4 + bio13 + I(bio1^2) + I(bio4^2) + I(bio13^2)",
                      "~ bio1 + bio4 + bio5 + bio13 + I(bio1^2) + I(bio4^2) + I(bio5^2) + I(bio13^2)",
                      "~ bio4 + bio5 + bio13 + I(bio4^2) + I(bio5^2) + I(bio13^2)"
                      )


## Ajuste de modelo
cal.cal.ppms <- foreach(i = seq_along(cal.cal.formulas))%do% {
      ppm(spp.ppp[[1]],
                trend = formula(cal.cal.formulas[i]),
                covariates = bio.im.list[[1]])
}

cal.mel.ppms <- foreach(i = seq_along(cal.mel.formulas))%do% {
   ppm(spp.ppp[[2]],
       trend = formula(cal.mel.formulas[i]),
       covariates = bio.im.list[[2]])
}

#Diagnóstico de los modelos

cal.cal.aic <- sapply(cal.cal.ppms, AIC)
cal.mel.aic <- sapply(cal.mel.ppms, AIC)

summary(cal.cal.ppms[[which.min(cal.cal.aic)]])
summary(cal.mel.ppms[[6]])

cal.cal.best <- cal.cal.ppms[[which.min(cal.cal.aic)]]
cal.mel.best <- cal.mel.ppms[[6]]

png("Real-spp/Diag-cal-cal.png", width = 1000, height = 1000)
diagnose.ppm(cal.cal.best)
dev.off()

png("Real-spp/Diag-cal-mel.png", width = 1000, height = 1000)
diagnose.ppm(cal.mel.best)
dev.off()

cal.cal.pred <- predict(cal.cal.best, type = "intensity", dimyx = c(116, 168))
cal.mel.pred <- predict(cal.mel.best, type = "intensity", dimyx = c(178, 194))

cal.cal.pred.r <- raster(cal.cal.pred)
cal.mel.pred.r <- raster(cal.mel.pred)

dir.create("Real-spp/Predictions")
writeRaster(cal.cal.pred.r, "Real-spp/Predictions/Callipepla-californica-PPM-NAD83", "GTiff")
writeRaster(cal.mel.pred.r, "Real-spp/Predictions/Calamospiza-melanocorys-PPM-NAD83", "GTiff")

########################
####Fitting ellipses####
########################



###########################
####Comparing centroids####
###########################
s.mod <- step(mod)
coef <- coefficients(mod)
rast <- raster(predict(mod, type = "trend", ngrid = c(209, 259)))


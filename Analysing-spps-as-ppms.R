library(raster); library(doParallel); library(spatstat)

l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
config <- readRDS("Simulated-species/Sim-config-species-list.rds")
spp.layers <- lapply(1:ncol(config$layers), function(x){dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config$layers[, x])))})
spp.points <- readRDS("Simulated-species/Species-presences.rds")
p.spp <- readRDS("Simulated-species/P-presence.rds")
spp.cent.cov <- readRDS("Simulated-species/Spp-cent-covs.rds")


### Formatting data for spatstat

spp.lay.df <- lapply(spp.layers, function(x){data.frame(rasterToPoints(x))}) #Transform layers to dataframes
spp.point.vals <- lapply(seq_along(spp.layers), function(x){extract(spp.layers[[x]], spp.points[[x]])}) #Extracting values at points


## Coercing raster package data to spatstat
ux = sort(unique(spp.lay.df[[1]]$x)) #Extracting unique coordinates
uy = sort(unique(spp.lay.df[[1]]$y))
nx = length(ux) #length of unique coordinates
ny = length(uy)
ref.cols = match(spp.lay.df[[1]]$x, ux) #position of every data point
ref.lines = match(spp.lay.df[[1]]$y, uy)
vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
vec[ref.vec] = 1
data.mask = matrix(vec, max(ref.lines), max(ref.cols), dimnames = list(uy, ux))
win = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window

## transforming species points to a planar point pattern
spp.ppp <- lapply(spp.points, function(x){
            ppp(x[, 'x'], x[, 'y'], window = win, check = F)
      }) 

## Transforming species' layers to spatstat images
spp.lay.im <- lapply(spp.lay.df, function(x){
      require(foreach)
      names(x) <- c("x", "y", "a", "b", "c")
      X <- with(x, cbind(a, b, c))
      lay.im.list <- foreach(i = 1:ncol(X)) %do% {
            vec.all = rep(NA, max(ref.lines)*max(ref.cols))
            vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
            vec.all[ref.vec] = X[,i]
            lay <- im(matrix(vec.all, max(ref.lines), max(ref.cols),
                              dimnames = list(uy, ux)), xcol = ux, yrow = uy)
            return(lay)
      }
      names(lay.im.list) <- c("a", "b", "c")
      return(lay.im.list)
})

## Fitting ppms
spp.ppms <- lapply(seq_along(spp.ppp), function(i){
      ppm.mod <-ppm(spp.ppp[[i]],
                    trend = ~ a + b + c + 
                          I(a^2) + I(b^2) + I(c^2), 
                             covariates = spp.lay.im[[i]])
      return(ppm.mod)
})

dir.create("../Resultados/Analysis-centroids/Fitted-PPMs")

for(i in seq_along(spp.ppms)){ 
   saveRDS(spp.ppms[[i]], paste0("../Resultados/Analysis-centroids/Fitted-PPMs/PPM-", i, ".rds"))
}

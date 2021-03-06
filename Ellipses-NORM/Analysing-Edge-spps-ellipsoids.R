library(raster); library(ntbox); library(doParallel)

all.layers <- readRDS("Simulated-layers/All-simulated-layers-NORM.rds")
config <- readRDS("Simulated-species/Sim-config-species-list-EdgeCentroids.rds")
l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
spp.points <- readRDS("Simulated-species/Species-presences-EdgeCentroids.rds")
mahal.dists <- readRDS("Simulated-species/Mahal-dists-EdgeCentroids.rds")
p.spp <- readRDS("Simulated-species/P-presence-EdgeCentroids.rds")

spp.layers <- lapply(1:ncol(config$layers), function(x){dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config$layers[, x])))})


registerDoParallel(cores = 2)
all.ellips <- foreach(i = seq_along(spp.points)) %dopar% {
      data <- data.frame(extract(spp.layers[[i]], spp.points[[i]]))
      data <- na.omit(data)
      names(data) <- c("a", "b", "c")
      cent <- cov_center(data, vars = c("a", "b", "c"), level = 0.99)
      ellip <- ellipsoidfit(envlayers = spp.layers[[i]], 
                            centroid = cent$centroid,
                            covar = cent$covariance,
                            level = 0.99,
                            size = 1, plot = F)
      return(list(cent, ellip))
}


dir.create("../Resultados/Analysis-centroids/Fitted-Edge-Ellipses-NORM")

for(i in seq_along(all.ellips)){ 
      saveRDS(all.ellips[[i]], paste0("../Resultados/Analysis-centroids/Fitted-Edge-Ellipses-NORM/Ellips-", i, ".rds"))
}


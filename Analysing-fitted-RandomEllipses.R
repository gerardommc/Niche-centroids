library(raster); library(ntbox); library(doParallel)

all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
spp.points <- readRDS("Simulated-species/Species-presences-RandomCentroids.rds")
p.spp <- readRDS("Simulated-species/P-presence-RandomCentroids.rds")
config <- readRDS("Simulated-species/Sim-config-species-list-RandomCentroids.rds")
ellips <- lapply(1:1000, function(x){
      readRDS(paste0("../Resultados/Analysis-centroids/Fitted-Random-Ellipses/Ellips-", x, ".rds"))
})
spp.layers <- lapply(1:ncol(config$layers), function(x){dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config$layers[, x])))})

spp.cent.cov <- readRDS("Simulated-species/Spp-cent-covs-RandomCentroids.rds")

registerDoParallel(cores = 4)
corr.p.pres <- foreach(i = seq_along(ellips)) %dopar% {
      cor.test(p.spp[[i]][], ellips[[i]][[2]]$suitRaster[])
}

corr.estimates <- sapply(corr.p.pres, function(x){x$estimate})

dist.true.cents <- sapply(seq_along(ellips), function(x){
      mahalanobis(ellips[[x]][[1]]$centroid, spp.cent.cov$centroids[[x]], spp.cent.cov$covariances[[x]])
})

#Saving results

df.centroids <- data.frame(Corr.surf = corr.estimates, Dist.true.cent = dist.true.cents)

vars.spp <- foreach(i = seq_along(config$layer.names), .combine = rbind) %do% {
      vars <- sapply(c("mean", "Log", "Beta", "Gamma"), function(x){
            grep(pattern = x, config$layer.names[[i]])
      })
      return(sapply(vars, length))
}
vars.spp <- data.frame(vars.spp)
names(vars.spp) <- c("Normal", "Log.norm", "Beta", "Gamma")

df.results <- data.frame(df.centroids, vars.spp, approach = "Ellipses", centr.conf = "random")

write.csv(df.results, "Simulated-species/Results-RandomEllipses.csv", row.names = F)



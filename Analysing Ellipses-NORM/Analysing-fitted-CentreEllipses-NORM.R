library(raster); library(ntbox); library(doParallel); library(bestNormalize)

all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
spp.points <- readRDS("Simulated-species/Species-presences.rds")
p.spp <- readRDS("Simulated-species/P-presence.rds")
config <- readRDS("Simulated-species/Sim-config-species-list.rds")
ellips <- lapply(1:2500, function(x){
      readRDS(paste0("../Resultados/Analysis-centroids/Fitted-ellipses-NORM/Ellips-", x, ".rds"))
      })

#Normalising the layers in order to transform the original centroids
registerDoParallel(cores = 8)
all.layers.norm <- foreach(i = 1:nlayers(all.layers)) %dopar% {
   df <- all.layers[[i]][]
   
   n <- bestNormalize(df)
   df.t <- predict(n)
   
   r <- all.layers[[i]]
   r[] <- matrix(df.t, 100, 100, byrow = T)
   return(list(normLayer = r, normMod = n))
}
normMods <- lapply(all.layers.norm, function(x){x$normMod})

all.layers.norm <- lapply(all.layers.norm, function(x){x$normLayer})
all.layers.norm <- stack(all.layers.norm)
names(all.layers.norm) <- names(all.layers)

#####
spp.layers <- lapply(1:ncol(config$layers), function(x){dropLayer(all.layers.norm, i = c(which(! 1:nlayers(all.layers.norm) %in% config$layers[, x])))})

### Retrieving true centroids
spp.cent.cov <- readRDS("Simulated-species/Spp-cent-covs.rds")
centroids <- covariances <- lapply(ellips, function(x){x[[1]]$centroid})

### Predicting true centroids on the normalised scale

centroids.norm.inv <- foreach(i = seq_along(centroids)) %dopar% {
   ids <- which(1:nlayers(all.layers.norm) %in% config$layers[, i])
   mods <- normMods[ids]
   cent <- centroids[[i]]
   new.cents <- sapply(1:3, function(x){predict(mods[[x]], newdata = cent[x], inverse = T)})
   return(new.cents)
}

### Estimating similarities

registerDoParallel(cores = 4)
corr.p.pres <- foreach(i = seq_along(ellips)) %dopar% {
      cor.test(p.spp[[i]][], ellips[[i]][[2]]$suitRaster[])
}

corr.estimates <- sapply(corr.p.pres, function(x){x$estimate})

dist.true.cents <- sapply(seq_along(ellips), function(x){
   mahalanobis(centroids.norm.inv[[x]], spp.cent.cov$centroids[[x]], spp.cent.cov$covariances[[x]])
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

df.results <- data.frame(df.centroids, vars.spp, approach = "Ellipses-NORM", centr.conf = "centre")

write.csv(df.results, "Simulated-species/Results-CentreEllipses-NORM.csv", row.names = F)

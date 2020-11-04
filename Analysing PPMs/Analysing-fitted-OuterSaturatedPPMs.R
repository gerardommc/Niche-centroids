library(raster); library(doParallel)

l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
config <- readRDS("Simulated-species/Sim-config-species-list-OuterCentroids.rds")
spp.layers <- lapply(1:ncol(config$layers), function(x){dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config$layers[, x])))})
spp.points <- readRDS("Simulated-species/Species-presences-OuterCentroids.rds")
p.spp <- readRDS("Simulated-species/P-presence-OuterCentroids.rds")
spp.cent.cov <- readRDS("Simulated-species/Spp-cent-covs-OuterCentroids.rds")

spp.ppms <- lapply(1:2500, function(x){
      readRDS(paste0("../Resultados/Analysis-centroids/Fitted-Outer-Saturated-PPMs/PPM-", x, ".rds"))
})

ppm.preds <- lapply(spp.ppms, function(x){x$pred})

cor.ppm.preds <- lapply(1:2500, function(x){cor.test(ppm.preds[[x]][], p.spp[[x]][])})

corr.estimates <- sapply(cor.ppm.preds, function(x){x$estimate})

## Computing centroids

centroids <- lapply(spp.ppms, function(x){
      effects <- x$coef[2:7]
      cent.a <- - effects[1]/(2 * effects[4])
      cent.b <- - effects[2]/(2 * effects[5])
      cent.c <- - effects[3]/(2 * effects[6])
      centroid <- c(a = cent.a,
                    b = cent.b,
                    c = cent.c)
      return(centroid)
})

#Calculating the distance to true centroids

dist.true.cents <- sapply(seq_along(centroids), function(x){
      mahalanobis(centroids[[x]], spp.cent.cov$centroids[[x]], spp.cent.cov$covariances[[x]])
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

df.results <- data.frame(df.centroids, vars.spp, approach = "PPM-sat", centr.conf = "outer")

write.csv(df.results, "Simulated-species/Results-OuterSaturatedPPMs.csv", row.names = F)

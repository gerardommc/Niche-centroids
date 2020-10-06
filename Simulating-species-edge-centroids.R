library(raster)

all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
l.sum <- read.csv("Simulated-layers/Layer-summaries.csv", stringsAsFactors = F)

total.combinations <- combn(x = 100, m = 3)
s <- sample(1:ncol(total.combinations), 1000) 
sample.layers <- total.combinations[, s]

layer.names <- lapply(1:ncol(sample.layers), function(x){l.sum$Layer.name[sample.layers[, x]]})

species.config.list <- list(Readme = "This list containes three other objects \n
                            1- The total possible number of three layer combinations out of
                            the 100 simulated layers. \n
                            2- An index of the random sample of the entire number of combinations and \n
                            3- The matrix of combinations \n
                            4- The names of each layer in the combination",
                            combinations = total.combinations,
                            sample = s,
                            layers = sample.layers,
                            layer.names = layer.names)

dir.create("Simulated-species")
saveRDS(species.config.list, "Simulated-species/Sim-config-species-list-EdgeCentroids.rds")

####
species.layers <- lapply(1:ncol(sample.layers), function(x){dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% sample.layers[, x])))})


species.centroids <- lapply(species.layers, 
                            function(x){
                                  q <- sample(c(0.05, 0.95), 3, replace = T)
                                  cent <- sapply(1:3, function(x1){quantile(x[[x1]], q[x1])})
                                  names(cent) <- paste0(names(x), "-", names(cent))
                                  return(cent)
                            })

species.covariances <- lapply(1:ncol(sample.layers), function(x){
      df <- data.frame(rasterToPoints(species.layers[[x]]))[, 3:5]
      return(MASS :: cov.rob(df, method = "mve", nsamp = "sample", seed = 1234)$cov)})

spp.cent.cov <- list(centroids = species.centroids,
                     covariances = species.covariances)
saveRDS(spp.cent.cov, "Simulated-species/Spp-cent-covs-EdgeCentroids.rds")

mahal.dists <- lapply(1:ncol(sample.layers), function(x){
      df <- data.frame(rasterToPoints(species.layers[[x]]))
      dist <- mahalanobis(df[, 3:5], center = species.centroids[[x]], cov = species.covariances[[x]])
      dist <- round(dist, 3)
      r <- rasterFromXYZ(data.frame(df[, c("x", "y")], dist))
      return(r)
})

saveRDS(mahal.dists, "Simulated-species/Mahal-dists-EdgeCentroids.rds")

p.presence <- lapply(mahal.dists, function(x){x1 <- log(1/x)
x1 <- x1 - cellStats(x1, min)
x1 <- x1/cellStats(x1, max)
return(x1)
})

p.presence <- lapply(mahal.dists, function(x){
   r <- 10/(cellStats(x, max) - cellStats(x, min))
   return(exp(-r * x))}
   )
saveRDS(p.presence, "Simulated-species/P-presence-EdgeCentroids.rds")

library(dismo)

spp.points <- lapply(p.presence, function(x){
   points <- randomPoints(mask = x, n = 500, prob = T)
   points <- t(apply(points, 1, function(x){ x + rnorm(2)}))
   return(points)
})
saveRDS(spp.points, "Simulated-species/Species-presences-EdgeCentroids.rds")


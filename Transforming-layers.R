library(MVN); library(doParallel); library(raster); library(bestNormalize)

all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
config <- lapply(paste0("Simulated-species/Sim-config-species-list", 
                        c("", "-EdgeCentroids",
                          "-OuterCentroids",
                          "-RandomCentroids"), ".rds"), readRDS)

spp.layers <- foreach(k = 1:ncol(config[[1]]$layers)) %do% {
      dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config[[1]]$layers[, k])))
}

#source("Normalise-functions.R")

registerDoParallel(cores = 8)
all.layers.norm <- foreach(i = 1:nlayers(all.layers)) %dopar% {
      df <- all.layers[[i]][]

      n <- bestNormalize(df)
      df.t <- predict(n)
      
      r <- all.layers[[i]]
      r[] <- matrix(df.t, 100, 100, byrow = T)
      return(list(normLayer = r, normMod = n))
}
normMods <- sapply(all.layers.norm, function(x){x$normMod})

all.layers.norm <- lapply(all.layers.norm, function(x){x$normLayer})
all.layers.norm <- stack(all.layers.norm)
names(all.layers.norm) <- names(all.layers)

saveRDS(all.layers.norm, "Simulated-layers/All-simulated-layers-NORM.rds")

spp.layers.norm <- foreach(k = 1:ncol(config[[1]]$layers)) %do% {
   dropLayer(all.layers.norm, i = c(which(! 1:nlayers(all.layers.norm) %in% config[[1]]$layers[, k])))
}

## Measuring skewness of transformed variables
s <- seq(1, 10000, by = 10)
registerDoParallel(cores = 4)
env.skew.norm <- foreach(i = seq_along(spp.layers.norm), .combine = c) %dopar% {
      df <- as.matrix(spp.layers.norm[[i]])[s, ]
      res <- MVN::mvn(df, mvnTest = "hz", multivariatePlot = F)$multivariateNormality$HZ
      return(res)
}

env.skew <- foreach(i = seq_along(spp.layers), .combine = c) %dopar% {
   df <- as.matrix(spp.layers[[i]])[s, ]
   res <- MVN::mvn(df, mvnTest = "hz", multivariatePlot = F)$multivariateNormality$HZ
   return(res)
}

df.sk <- data.frame(set = rep(1:2500, 2), sk = c(env.skew.norm, env.skew), type = rep(c("Corrected", "Uncorrected"), each = 2500))
write.csv(df.sk, "Simulated-layers/Skewness-difference-NORM-DATA.csv")

library(ggplot2)

pdf("../Graphs/Skewness-Corr-Uncorr.pdf", width = 4, height = 4)
ggplot(df.sk) + geom_boxplot(aes(x = type, y = sk)) + 
   labs(x = "Variables", y = "Multivariate skewness")
dev.off()

t.test(env.skew, env.skew.norm, paired = T)

length(which((env.skew.norm - env.skew) < 0))



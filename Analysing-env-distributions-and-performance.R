library(data.table); library(doParallel); library(reshape)
library(raster); library(ggplot2)

#Analysis of models
ppms <- rbindlist(lapply(list.files("Simulated-species", "PPMs", full.names = T), read.csv))
ellips <-  rbindlist(lapply(list.files("Simulated-species", "llipses", full.names = T), read.csv))
all.results <- rbind(ppms, ellips)

#Model subsets
ppms.sat <- subset(ppms, approach == "PPM-sat")
ppms.step <- subset(ppms, approach == "PPM-step")

corr.dif.sat <- ppms.sat$Corr.surf - ellips$Corr.surf
dist.dif.sat <- ellips$Dist.true.cent - ppms.sat$Dist.true.cent

#Environmental layers
all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
config <- lapply(paste0("Simulated-species/Sim-config-species-list", 
                        c("", "-EdgeCentroids",
                          "-OuterCentroids",
                          "-RandomCentroids"), ".rds"), readRDS)

spp.layers <- foreach(k = 1:ncol(config[[1]]$layers)) %do% {
      dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config[[1]]$layers[, k])))
}

pres <- lapply(paste0("Simulated-species/Species-presences",
                      c(".rds",
                        "-EdgeCentroids.rds",
                        "-OuterCentroids.rds",
                        "-RandomCentroids.rds")), readRDS)

#Skewness analyses
registerDoParallel(cores = 4)

s <- seq(1, 10000, by = 10) #Sample of the environment

env.skew <- foreach(i = seq_along(spp.layers), .combine = c) %dopar% {
      df <- as.matrix(spp.layers[[i]])[s, ]
      res <- MVN::mvn(df, mvnTest = "hz", multivariatePlot = F)$multivariateNormality$HZ
      return(res)
}

#Adding the analyses to model data
ppms.sat$Skewness.env <- rep(env.skew, 4)
ppms.step$Skewness.env <- rep(env.skew, 4)
ellips$Skewness.env <-  rep(env.skew, 4)

dat.skew <- rbind(ppms.sat, ellips)

dir.create("../Graphs")
#Plots of th relationships
pdf("../Graphs/Env-skewness-Distance-cent.pdf", width = 9, height = 4)
ggplot(dat.skew) + geom_hex(aes(x = Skewness.env, y = log10(Dist.true.cent))) + 
      scale_fill_continuous(type = "viridis", trans = "log10") +
      geom_smooth(aes(x = Skewness.env, y = log10(Dist.true.cent)), colour = "red", alpha = 0.5) +
      labs(fill = expression(paste(log[10], " count")), colour = "", 
           x = "Multivariate environmental skewness",
           y = expression(paste(log[10], " Distance to true centroid"))) +
      facet_grid(rows = vars(approach), cols = vars(centr.conf))
dev.off()

pdf("../Graphs/Env-skewness-Corr-surf.pdf", width = 9, height = 4)
ggplot(dat.skew) + geom_hex(aes(x = Skewness.env, y = Corr.surf)) + 
      scale_fill_continuous(type = "viridis", trans = "log10") +
      geom_smooth(aes(x = Skewness.env, y = Corr.surf),  colour = "red", alpha = 0.5) +
      labs(fill = expression(paste(log[10], " count")), colour = "", 
           x = "Multivariate environmental skewness",
           y = "Correlation with generating surface") +
      facet_grid(rows = vars(approach), cols = vars(centr.conf))
dev.off()

#Analysis of the frequency of positive squared terms B'

ppm.files <- paste0("../Resultados/Analysis-centroids/Fitted-", 
                    rep(c("Centre", "Edge", "Outer", "Random"), each = 2500), 
                    "-Saturated-PPMs/PPM-", rep(1:2500, 4), ".rds")
ppm.mods <- lapply(ppm.files, readRDS)


#ellip.files <- paste0("../Resultados/Analysis-centroids/Fitted-", 
                      #rep(c("Centre", "Edge", "Outer", "Random"), each = 2500), 
                      #"-Ellipses/Ellips-", rep(1:2500, 4), ".rds")
#ellip.mods <- lapply(ellip.files, readRDS)

   
neg.coefs <- sapply(ppm.mods, function(x){(any(x$coef[5:7] > 0))})

ppms.sat$pos.neg <- neg.coefs

ggplot(ppms.sat) + geom_boxplot(aes(x = pos.neg, y = log10(Dist.true.cent))) +
   facet_grid(~ centr.conf)

ggplot(ppms.sat) + geom_bar(aes(x = pos.neg), position = "dodge", alpha = 0.5) +
   facet_grid(~ centr.conf)

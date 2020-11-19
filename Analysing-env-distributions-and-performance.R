library(data.table); library(doParallel); library(reshape)
library(raster); library(ggplot2)

#Analysis of models
ppms <- rbindlist(lapply(list.files("Simulated-species", "PPMs", full.names = T), read.csv))
ellips <-  rbindlist(lapply(list.files("Simulated-species", "llipses", full.names = T), read.csv))
all.results <- rbind(ppms, ellips)

#Model subsets
ppms.sat <- subset(ppms, approach == "PPM-sat")
ppms.sat.norm <- subset(ppms, approach == "PPM-sat-NORM")
ppms.step <- subset(ppms, approach == "PPM-step")

ellips.norm <- subset(ellips, approach == "Ellipses-NORM")
ellips <- subset(ellips, approach == "Ellipses")

corr.dif.sat <- ppms.sat$Corr.surf - ellips$Corr.surf
dist.dif.sat <- ellips$Dist.true.cent - ppms.sat$Dist.true.cent

#Environmental layers
all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
all.layers.norm <- readRDS("Simulated-layers/All-simulated-layers-NORM.rds")

config <- lapply(paste0("Simulated-species/Sim-config-species-list", 
                        c("", "-EdgeCentroids",
                          "-OuterCentroids",
                          "-RandomCentroids"), ".rds"), readRDS)

library(foreach)

bimod <- lapply(1:nlayers(all.layers), function(x){diptest::dip.test(all.layers[[x]][])})
bimod.norm <- lapply(1:nlayers(all.layers.norm), function(x){diptest::dip.test(all.layers.norm[[x]][])})

spp.layers <- foreach(k = 1:ncol(config[[1]]$layers)) %do% {
      dropLayer(all.layers, i = c(which(! 1:nlayers(all.layers) %in% config[[1]]$layers[, k])))
}

spp.layers.norm <- foreach(k = 1:ncol(config[[1]]$layers)) %do% {
   dropLayer(all.layers.norm, i = c(which(! 1:nlayers(all.layers.norm) %in% config[[1]]$layers[, k])))
}

comp.bimod <- foreach(k = 1:ncol(config[[1]]$layers), .combine = c) %do% {
   b <- bimod[c(which(1:nlayers(all.layers) %in% config[[1]]$layers[, k]))]
   sum(sapply(b, function(x){x$statistic}))
}

comp.bimod.norm <- foreach(k = 1:ncol(config[[1]]$layers), .combine = c) %do% {
   b <- bimod.norm[c(which(1:nlayers(all.layers.norm) %in% config[[1]]$layers[, k]))]
   sum(sapply(b, function(x){x$statistic}))
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

env.skew.norm <- foreach(i = seq_along(spp.layers.norm), .combine = c) %dopar% {
   df <- as.matrix(spp.layers.norm[[i]])[s, ]
   res <- MVN::mvn(df, mvnTest = "hz", multivariatePlot = F)$multivariateNormality$HZ
   return(res)
}

#Adding the analyses to model data
ppms.sat$Skewness.env <- rep(env.skew, 4)
ppms.sat.norm$Skewness.env <- rep(env.skew.norm, 4)
ppms.step$Skewness.env <- rep(env.skew, 4)

ellips$Skewness.env <-  rep(env.skew, 4)
ellips.norm$Skewness.env <-  rep(env.skew.norm, 4)

ppms.sat$bimod <-rep(comp.bimod, 4)
ppms.sat.norm$bimod <- rep(comp.bimod.norm, 4)
ppms.step$bimod <- rep(comp.bimod, 4)

ellips$bimod <- rep(comp.bimod, 4)
ellips.norm$bimod <- rep(comp.bimod.norm, 4)

dat.skew <- rbind(ppms.sat, ellips)
dat.skew.norm <- rbind(ppms.sat.norm, ellips.norm)

dat.skew.all <- rbind(dat.skew, dat.skew.norm)

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

ggplot(dat.skew.norm) + geom_hex(aes(x = Skewness.env, y = log10(Dist.true.cent))) + 
   scale_fill_continuous(type = "viridis", trans = "log10") +
   geom_smooth(aes(x = Skewness.env, y = log10(Dist.true.cent)), colour = "red", alpha = 0.5) +
   labs(fill = expression(paste(log[10], " count")), colour = "", 
        x = "Multivariate environmental skewness",
        y = expression(paste(log[10], " Distance to true centroid"))) +
   facet_grid(rows = vars(approach), cols = vars(centr.conf))



pdf("../Graphs/Env-skewness-Corr-surf.pdf", width = 9, height = 4)
ggplot(dat.skew) + geom_hex(aes(x = Skewness.env, y = Corr.surf)) + 
      scale_fill_continuous(type = "viridis", trans = "log10") +
      geom_smooth(aes(x = Skewness.env, y = Corr.surf),  colour = "red", alpha = 0.5) +
      labs(fill = expression(paste(log[10], " count")), colour = "", 
           x = "Multivariate environmental skewness",
           y = "Correlation with generating surface") +
      facet_grid(rows = vars(approach), cols = vars(centr.conf))
dev.off()

ggplot(dat.skew.all) + geom_hex(aes(x = Skewness.env, y = Corr.surf)) + 
   scale_fill_continuous(type = "viridis", trans = "log10") +
   geom_smooth(aes(x = Skewness.env, y = Corr.surf),  colour = "red", alpha = 0.5) +
   labs(fill = expression(paste(log[10], " count")), colour = "", 
        x = "Multivariate environmental skewness",
        y = "Correlation with generating surface") +
   facet_grid(rows = vars(approach), cols = vars(centr.conf))

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

pdf("../Graphs/NegCoef-DistCent.pdf", width = 9, height = 4.5)
ggplot(ppms.sat) + geom_point(aes(x = pos.neg, y = log10(Dist.true.cent), 
                                  colour = Skewness.env),
                              position = position_jitter(width = 0.2), alpha = 0.3) +
   scale_colour_continuous(type = "viridis") + 
   geom_boxplot(aes(x = pos.neg, y = log10(Dist.true.cent)), alpha = 0.5) +
   facet_grid(~ centr.conf) +
   labs(x = "Presence of positive squared terms",
        y = "Distance to true centroid", 
        colour = "Skewness")
dev.off()

pdf("../Graphs/NegCoef-NumModels.pdf", width = 9, height = 4.5)
ggplot(ppms.sat) + geom_bar(aes(x = pos.neg), position = "dodge", alpha = 0.5) +
   facet_grid(~ centr.conf) +
   labs(x = "Presence of positive squared terms",
        y = "Number of models")
dev.off()

pdf("../Graphs/NegCoef-Skewness.pdf", width = 9, height = 4.5)
ggplot(ppms.sat) + geom_boxplot(aes(x = pos.neg, y = Skewness.env)) + 
   geom_violin(aes(x = pos.neg, y = Skewness.env), alpha = 0.3) +
   facet_grid(~ centr.conf) +
   labs(x = "Presence of positive squared terms", 
        y = "Multivariate environmental skewness")
dev.off()

#Statistics for skewness data & vis
skew.stats <- reshape :: melt(ppms.sat[1:2500,], 
                             measure.vars = c("Normal", "Log.norm", "Gamma", "Beta"),
                             id.vars = c("Skewness.env", "Corr.surf", "Dist.true.cent"))

library(dplyr)

skew.summ <- skew.stats %>% group_by(variable, value) %>% 
   summarise(Skewness.env = mean(Skewness.env))

pdf("../Graphs/Skewness-distributions.pdf", width = 6, height = 4.5)
ggplot(skew.summ) + geom_tile(aes(x = value, y = variable, fill = Skewness.env)) + 
   scale_fill_viridis_c() + 
   labs(x = "Number of dimensions",
        y = "Distribution", fill = "Skewness")
dev.off()


#Statistics for bimodality data & vis
bimod.stats <- reshape :: melt(ppms.sat[1:2500,], 
                              measure.vars = c("Normal", "Log.norm", "Gamma", "Beta"),
                              id.vars = c("bimod", "Corr.surf", "Dist.true.cent"))

bimod.summ <- bimod.stats %>% group_by(variable, value) %>% 
   summarise(bimod = mean(bimod))

pdf("../Graphs/Bimodality-distributions.pdf", width = 6, height = 4.5)
ggplot(bimod.summ) + geom_tile(aes(x = value, y = variable, fill = bimod)) + 
   scale_fill_viridis_c() + 
   labs(x = "Number of dimensions",
        y = "Distribution",
        fill = "Bimodality")
dev.off()

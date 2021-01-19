library(raster); library(spatstat); library(foreach)

##Reading data
# Simulation config
spp.dat <- readRDS("../Resultados/Interacting-spp-suitability/Simulated-interacting-spp.rds")

#Environmental data
env <- stack(paste0("../Chelsa/CHELSA_bio10_", c("01", "02", "12", "14"), ".tif"))
env <- crop(env, y = extent(c(-105, -85, 7.5, 22.5)))
env <- aggregate(env, 5)
names(env) <- paste0("bio", c(1, 2, 12, 14))

env.ranges <- data.frame(cellStats(env, range))

t.mean <- with(env.ranges, seq(bio1[1]/10, bio1[2]/10, len = 1000))
t.ran <- with(env.ranges, seq(bio2[1]/10, bio2[2]/10, len = 1000))
rain <- with(env.ranges, seq(bio12[1], bio12[2], len = 1000))
rain.dry <- with(env.ranges, seq(bio14[1], bio14[2], len = 1000))

##Extracting true centroids
source("Functions-sim-true-spp/Sim-spp-funcitons.R")
#Species 1
t.perf <- temp.optim(x = t.mean,
                     k = list(k1 = 1, k2 =0, k5 = 40,
                              k3 = 2, k4 = 28.3,
                              k6 = 30, k7 = -0.02))
t.ran.perf <- 0.3 * t.ran * exp(-0.3 * t.ran)
r.perf <- rain.resp(x = rain, k = list(k1 = 1, k2 = 0.001))
r.dry.perf <- 0.05 * rain.dry * exp(-0.05 * rain.dry)

#Species 2
t.perf.2 <- temp.optim(x = t.mean, k = list(k1 = 1, k2 = 0, k5 = 40,
                                                k3 = 2, k4 = 28.3, 
                                                k6 = 30, k7 = -0.02))
t.ran.perf.2 <-  0.2 * t.ran *  exp(-0.2 * t.ran)
r.perf.2 <- 0.0004 * rain * exp(-0.0004* rain)
r.dry.perf.2 <- rain.resp(x = rain.dry, k = list(k1 = 1, k2 = 0.02))

## True centroids
#Species 1
temp.1 <-  t.mean[which.max(t.perf)]
t.ran.1 <- t.ran[which.max(t.ran.perf)]
rain.1 <- rain[which.max(r.perf)]
rain.dry.1 <- rain.dry[which.max(r.dry.perf)]

#Species 2
temp.2 <-  t.mean[which.max(t.perf.2)]
t.ran.2 <- t.ran[which.max(t.ran.perf.2)]
rain.2 <- rain[which.max(r.perf.2)]
rain.dry.2 <- rain.dry[which.max(r.dry.perf.2)]


par(mfrow = c(2,2))
plot(t.mean, t.perf); abline(v = temp.1, col = "red")
plot(t.ran, t.ran.perf); abline(v = t.ran.1, col = "red")
plot(rain, r.perf); abline(v = rain.1, col = "red")
plot(rain.dry, r.dry.perf); abline(v = rain.dry.1, col = "red")

plot(t.mean, t.perf.2); abline(v = temp.2, col = "red")
plot(t.ran, t.ran.perf.2); abline(v = t.ran.2, col = "red")
plot(rain, r.perf.2); abline(v = rain.2, col = "red")
plot(rain.dry, r.dry.perf.2); abline(v = rain.dry.2, col = "red")

## Loading Fitted models
hard.mod <- readRDS("../Resultados/Interacting-spp-suitability/Harcore-Models.rds")
pois.mod <- readRDS("../Resultados/Interacting-spp-suitability/Poisson-Models.rds")
ellip.mod <- readRDS("../Resultados/Interacting-spp-suitability/Ellisoid-models.rds")

#Extracting centroids
hard.cents <- foreach(i = seq_along(hard.mod), .combine = rbind) %do% {
      betas <- hard.mod[[i]]$coef
      betas.1 <- betas[seq(3, 17, by = 2)]
      betas.2 <- betas[seq(4, 18, by = 2)]
      cents.1 <- - betas.1[seq(1, 7, by = 2)]/(2 * betas.1[seq(2, 8, by = 2)])
      cents.2 <- - betas.2[seq(1, 7, by = 2)]/(2 * betas.2[seq(2, 8, by = 2)])
      cent <- data.frame(rbind(cents.1, cents.2))
      names(cent) <- c("bio1", "bio2", "bio12", "bio14")
      cent$spp <- c("spp.1", "spp.2")
      cent$method <- "hard"
      cent$rad <- c(spp.dat$hradii$spp.1[i], spp.dat$hradii$spp.2[i])
      return(cent)
}

pois.cents <- foreach(i = seq_along(pois.mod), .combine = rbind) %do% {
      betas <- pois.mod[[i]]$coef
      betas.1 <- betas[seq(3, 17, by = 2)]
      betas.2 <- betas[seq(4, 18, by = 2)]
      cents.1 <- - betas.1[seq(1, 7, by = 2)]/(2 * betas.1[seq(2, 8, by = 2)])
      cents.2 <- - betas.2[seq(1, 7, by = 2)]/(2 * betas.2[seq(2, 8, by = 2)])
      cent <- data.frame(rbind(cents.1, cents.2))
      names(cent) <- c("bio1", "bio2", "bio12", "bio14")
      cent$spp <- c("spp.1", "spp.2")
      cent$method <- "pois"
      cent$rad <- c(spp.dat$hradii$spp.1[i], spp.dat$hradii$spp.2[i])
      return(cent)
}

ellip.cents <- foreach(i = seq_along(ellip.mod$spp.1$cent.cov), .combine = rbind) %do% {
   ellip.cents.1 <- ellip.mod$spp.1$cent.cov[[i]]$centroid
   ellip.cents.2 <- ellip.mod$spp.2$cent.cov[[i]]$centroid
   cent <- data.frame(rbind(ellip.cents.1, ellip.cents.2))
   names(cent) <- c("bio1", "bio2", "bio12", "bio14")
   cent$spp <- c("spp.1", "spp.2")
   cent$method <- "ellip"
   cent$rad <- c(spp.dat$hradii$spp.1[i], spp.dat$hradii$spp.2[i])
   return(cent)
}

centroids <- rbind(hard.cents, pois.cents, ellip.cents)

# Calculating distance to true centroids and correlation with generating surface
cents.spp.1 <- subset(centroids, spp == "spp.1")
cents.spp.2 <- subset(centroids, spp == "spp.2")

samps <- data.frame(rasterToPoints(env))[, c("bio1", "bio2", "bio12", "bio14")]
cov.mat <- MASS :: cov.mve(samps[sample(1:nrow(samps), 5000),])

dists.1 <- foreach(i = 1:nrow(cents.spp.1), .combine = c) %do% {
   mahalanobis(with(cents.spp.1, c(bio1[i], bio2[i], bio12[i], bio14[i])), 
               c(temp.1 * 10, t.ran.1 * 10, rain.1, rain.dry.1),
               cov = cov.mat$cov)
}

dists.2 <- foreach(i = 1:nrow(cents.spp.2), .combine = c) %do% {
   mahalanobis(with(cents.spp.2, c(bio1[i], bio2[i], bio12[i], bio14[i])), 
               c(temp.2 * 10, t.ran.2 * 10, rain.2, rain.dry.2),
               cov = cov.mat$cov)
}

cents.spp.1$dist.t.cent <- dists.1
cents.spp.2$dist.t.cent <- dists.2

## Correlation with generating surface

hard.surfs <- lapply(hard.mod, function(x){
   r <- predict(x, dimyx = c(360, 480))
   r <- lapply(r, raster)
   return(stack(r))
   })
pois.surfs <- lapply(pois.mod, function(x){
   r <- predict(x, dimyx = c(360, 480))
   r <- lapply(r, raster)
   return(stack(r))
   })
ellip.surfs <- foreach(i = 1:25) %do% {
   r1 <- ellip.mod$spp.1$ellip[[i]]$suitRaster
   r2 <- ellip.mod$spp.2$ellip[[i]]$suitRaster
   return(stack(r1, r2))
}

suit.sim <- stack(raster(spp.dat$suit[[1]]), raster(spp.dat$suit[[2]]))
suit.sim <- data.frame(rasterToPoints(suit.sim))
names(suit.sim) <- c("x", "y", "spp.1", "spp.2")

hard.cors <- foreach(i = seq_along(hard.mod), .combine = rbind) %do% {
   r.df <- data.frame(extract(hard.surfs[[i]], suit.sim[, c("x", "y")]))
   c1 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   c2 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   return(c(spp.1 = c1$estimate,  spp.2 = c2$estimate))
}
pois.cors <- foreach(i = seq_along(hard.mod), .combine = rbind) %do% {
   r.df <- data.frame(extract(pois.surfs[[i]], suit.sim[, c("x", "y")]))
   c1 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   c2 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   return(c(spp.1 = c1$estimate,  spp.2 = c2$estimate))
}
ellip.cors <- foreach(i = seq_along(hard.mod), .combine = rbind) %do% {
   r.df <- data.frame(extract(ellip.surfs[[i]], suit.sim[, c("x", "y")]))
   names(r.df) <- c("spp.1", "spp.2")
   c1 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   c2 <- cor.test(r.df$spp.1, suit.sim$spp.1, method = "spearman", continuiuty = T)
   return(c(spp.1 = c1$estimate,  spp.2 = c2$estimate))
}

hard.cors <- data.frame(hard.cors)
pois.cors <- data.frame(pois.cors)
ellip.cors <- data.frame(ellip.cors)

cents.spp.1$corr.surf <- c(hard.cors$spp.1.rho, pois.cors$spp.1.rho, ellip.cors$spp.1.rho)
cents.spp.2$corr.surf <- c(hard.cors$spp.2.rho, pois.cors$spp.2.rho, ellip.cors$spp.2.rho)

perf.stats <- rbind(cents.spp.1, cents.spp.2)

write.csv(perf.stats, "../Resultados/Interacting-spp-suitability/Performance-stats.csv")




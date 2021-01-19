library(raster)
source("Functions-sim-true-spp/Sim-spp-funcitons.R")

bioclim <- stack(paste0("../Chelsa/CHELSA_bio10_", c("01", "02", "12", "14"), ".tif"))
cr <- crop(bioclim, y = extent(c(-105, -85, 7.5, 22.5)))

sum.bioc <- data.frame(summary(cr))

#Species 1
t.perf <- temp.optim(x = cr[[1]]/10, k = list(k1 = 1, k2 =0, k5 = 40,
                                              k3 = 2, k4 = 28.3, 
                                              k6 = 30, k7 = -0.02))
t.ran.perf <- 0.3 * cr[[2]]/10 * exp(-0.3 * cr[[2]]/10)
r.perf <- rain.resp(cr[[3]], k = list(k1 = 1, k2 = 0.001))
r.dry.perf <- 0.05 * cr[[4]] * exp(-0.05 * cr[[4]])

beta <- runif(4, 0, 2)

perf.layers <- stack(t.perf, t.ran.perf, r.perf, r.dry.perf) 

#Species 2
t.perf.2 <- temp.optim(x = cr[[1]]/10, k = list(k1 = 1, k2 =0, k5 = 40,
                                                k3 = 2, k4 = 28.3, 
                                                k6 = 30, k7 = -0.02))
t.ran.perf.2 <- 0.2 * cr[[2]] * exp(-0.2 * cr[[2]]/10)
r.perf.2 <- 0.0004 * cr[[3]] * exp(-0.0004* cr[[3]])
r.dry.perf.2 <- rain.resp(x = cr[[4]], k = list(k1 = 1, k2 = 0.02))

beta.2 <- runif(4, 0, 2)

perf.layers.2 <- stack(t.perf.2, t.ran.perf.2, r.perf.2, r.dry.perf.2) 

dir.create("../Resultados/Interacting-spp-suitability")

saveRDS(list(beta.1 = beta, beta.2 = beta.2,
        perf.layers = perf.layers,
        perf.layers.2 = perf.layers.2),
        "../Resultados/Interacting-spp-suitability/Spp-config.rds")

par(mfrow = c(2, 1))
plot(exp(sum(perf.layers * beta)))
plot(exp(sum(perf.layers.2 * beta.2)))

library(raster); library(spatstat); library(ntbox); library(foreach)

## Reading and formatting environmental data
env <- stack(paste0("../Chelsa/CHELSA_bio10_", c("01", "02", "12", "14"), ".tif"))
env <- crop(env, y = extent(c(-105, -85, 7.5, 22.5)))
env <- aggregate(env, 5)
names(env) <- paste0("bio", c(1, 2, 12, 14))

env.df <- data.frame(rasterToPoints(env))

# Spatstat formatting
ux = sort(unique(env.df$x)) #Extracting unique coordinates
uy = sort(unique(env.df$y))
nx = length(ux) #length of unique coordinates
ny = length(uy)
ref.cols = match(env.df$x, ux) #position of every data point
ref.lines = match(env.df$y, uy)
vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
vec[ref.vec] = 1
data.mask = matrix(vec, max(ref.lines), max(ref.cols), dimnames = list(uy, ux))
win = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window

X <- as.matrix(env.df[, 3:6])
names(X) <- names(env.df)

im.list <- foreach(i = 1:ncol(X)) %do% {
      vec.all = rep(NA, max(ref.lines)*max(ref.cols))
      vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
      vec.all[ref.vec] = X[,i]
      lay <- im(matrix(vec.all, max(ref.lines), max(ref.cols),
                       dimnames = list(uy, ux)), xcol = ux, yrow = uy)
      return(lay)
}
names(im.list) <- names(env)

## Reading occurrences and Fundamental niches

spp.dat <- readRDS("../Resultados/Interacting-spp-suitability/Simulated-interacting-spp.rds")

multi.hard <- spp.dat$multi.hard

pdf("../Graphs/Suitability-quantess-interacting-spp.pdf", width = 10, height = 10)
for(i in seq_along(multi.hard)){
   par(mfrow= c(2,2))
   plot(quantess(split(multi.hard[[i]])[[1]], spp.dat$suit[[1]], n = 5),
        main = paste0("Spp 1, r =", spp.dat$hradii$spp.1[i]))
   points(split(multi.hard[[i]])[[1]], col = "white", pch = "+", cex = 0.75)
   plot(rhohat(split(multi.hard[[i]])[[1]], spp.dat$suit[[1]]), main = "",
        xlab = "Suitability", ylab = "Intensity")

   plot(quantess(split(multi.hard[[i]])[[2]], spp.dat$suit[[2]], n = 5),
        main = paste0("Spp 2, r =", spp.dat$hradii$spp.2[i]))
   points(split(multi.hard[[i]])[[2]], col = "white", pch = "+", cex = 0.75)
   plot(rhohat(split(multi.hard[[i]])[[2]], spp.dat$suit[[2]]), main = "",
        xlab = "Suitability", ylab = "Intensity")
}
dev.off()

library(doParallel)
registerDoParallel(cores = 2)

hard.mod <- foreach(i = seq_along(multi.hard)) %dopar%{
   ppm(multi.hard[[i]] ~ 0 + marks + marks : (bio1 + I(bio1^2) + 
                                      bio2 + I(bio2^2) + 
                                      bio12 + I(bio12^2) + 
                                      bio14 + I(bio14^2)),
                 covariates = im.list, MultiHard())
   }

pois.mod <- foreach(i = seq_along(multi.hard)) %dopar%{
   ppm(multi.hard[[i]] ~ 0 + marks + marks : (bio1 + I(bio1^2) + 
                                      bio2 + I(bio2^2) + 
                                      bio12 + I(bio12^2) + 
                                      bio14 + I(bio14^2)),
                covariates = im.list)
   }

saveRDS(hard.mod, "../Resultados/Interacting-spp-suitability/Harcore-Models.rds")
saveRDS(pois.mod, "../Resultados/Interacting-spp-suitability/Poisson-Models.rds")
## Ellipsoids

split.spp <- lapply(multi.hard, split)

spp.1 <- lapply(split.spp, function(x)x$spp.1)
spp.2 <- lapply(split.spp, function(x)x$spp.2)

env.spp.1 <- lapply(spp.1, function(x){extract(env, as.data.frame(x))})
env.spp.2 <- lapply(spp.2, function(x){extract(env, as.data.frame(x))})

spp.cent.1 <- foreach(i = seq_along(multi.hard)) %dopar% {
   cov_center(data = env.spp.1[[i]], vars = names(env), mve = T, level = 0.99)}
spp.cent.2 <- foreach(i = seq_along(multi.hard)) %dopar% {
   cov_center(data = env.spp.2[[i]], vars = names(env), mve = T, level = 0.99)}

spp.ellip.1 <- foreach(i = seq_along(multi.hard)) %dopar% {
                ellipsoidfit(envlayers = env, 
                             centroid = spp.cent.1[[i]]$centroid,
                             covar = spp.cent.1[[i]]$covariance,
                             level = 0.99,
                             size = 1, plot = F)}

spp.ellip.2 <- foreach(i = seq_along(multi.hard)) %dopar%{
                ellipsoidfit(envlayers = env, 
                             centroid = spp.cent.2[[i]]$centroid,
                             covar = spp.cent.2[[i]]$covariance,
                             level = 0.99,
                             size = 1, plot = F)}

saveRDS(list(spp.1 = list(cent.cov = spp.cent.1,
                          ellip = spp.ellip.1),
             spp.2 = list(cent.cov = spp.cent.2,
                          ellip = spp.ellip.2)),
        "../Resultados/Interacting-spp-suitability/Ellisoid-models.rds")

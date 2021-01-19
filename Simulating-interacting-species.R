library(spatstat); library(raster); library(foreach)

spp.conf <- readRDS("../Resultados/Interacting-spp-suitability/Spp-config.rds")

FN.1 <- exp(sum(spp.conf$beta.1 * spp.conf$perf.layers))
FN.2 <- exp(sum(spp.conf$beta.2 * spp.conf$perf.layers.2))

FN.1 <- aggregate(FN.1, 5)
FN.2 <- aggregate(FN.2, 5)

FN.1 <- (FN.1 - cellStats(FN.1, min))/(cellStats(FN.1, max) - cellStats(FN.1, min))
FN.2 <- (FN.2 - cellStats(FN.2, min))/(cellStats(FN.2, max) - cellStats(FN.2, min))

FN.stack <- stack(FN.1, FN.2)
FN.df <- data.frame(rasterToPoints(FN.stack))

# Spatstat formatting
ux = sort(unique(FN.df$x)) #Extracting unique coordinates
uy = sort(unique(FN.df$y))
nx = length(ux) #length of unique coordinates
ny = length(uy)
ref.cols = match(FN.df$x, ux) #position of every data point
ref.lines = match(FN.df$y, uy)
vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
vec[ref.vec] = 1
data.mask = matrix(vec, max(ref.lines), max(ref.cols), dimnames = list(uy, ux))
win = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window

X <- cbind(FN.df[, 3], FN.df[, 4])
names(X) <- c("Spp.1", "Spp.2")

im.list <- foreach(i = 1:ncol(X)) %do% {
      vec.all = rep(NA, max(ref.lines)*max(ref.cols))
      vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
      vec.all[ref.vec] = X[,i]
      lay <- im(matrix(vec.all, max(ref.lines), max(ref.cols),
                       dimnames = list(uy, ux)), xcol = ux, yrow = uy)
      return(lay)
}
names(im.list) <- c("Spp.1", "Spp.2")

weighted.im.list <- lapply(im.list, function(x){rpois(1, 10) * x})

hradii <- expand.grid(spp.1 = seq(0, 0.5, len = 5),
                      spp.2 = seq(0, 0.5, len = 5))

rag.points <- lapply(1:nrow(hradii), function(x){
   ragsMultiHard(beta = weighted.im.list,
                            hradii = matrix(c(0, hradii$spp.1[x], hradii$spp.2[x], 0), nrow = 2), 
                            types = c("spp.1", "spp.2"))})

n.spp.1 <- sapply(rag.points, function(x){split(x)$spp.1$n})
n.spp.2 <- sapply(rag.points, function(x){split(x)$spp.2$n})

pdf("../Graphs/Interacting-species.pdf", width = 21, height = 7)
for(i in seq_along(rag.points)){
   par(mfrow = c(1, 2), mar = c(1.5, 1.5, 1.5, 1.5))
   plot(weighted.im.list[[1]], main = paste0("Spp 1, r = ", 
                                             hradii$spp.1[i], 
                                             ", n = ", n.spp.1[i]))
   points(split(rag.points[[i]])[[1]], col = "white")
   
   plot(weighted.im.list[[2]], main = paste0("Spp 2, r = ", hradii$spp.2[i],
                                             ", n = ", n.spp.2[i]))
   points(split(rag.points[[i]])[[2]], col = "white")
}
dev.off()

saveRDS(list(multi.hard = rag.points,
             hradii = hradii,
             suit = weighted.im.list),
        "../Resultados/Interacting-spp-suitability/Simulated-interacting-spp.rds")

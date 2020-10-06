#This script simulates 100 raster images with different statistical distributions: Normal, log-normal, beta and gamma
#To use for simulating species 
library(mvtnorm)   # Simulate multivariate distributions
library(raster)    # Handle raster objects
library(doParallel)   # Advanced loops

#The function to estimate the distance from one spatial unit to the rest
mat.distances <- function(len)
{
      coords.row <- rep(1:len, times=len)
      coords.col <- rep(1:len, each=len)
      row.col <- data.frame(coords.row, coords.col)
      D <- dist(row.col, method="euclidean", diag=TRUE, upper=TRUE)
      D <- as.matrix(D)
      return(D)
}

#Establecemos la función que va a simular las variables (Sólamente traducida al español de la original de Petr Keil)

corr.surf <- function(len, global.mean, lambda, sigma){
      
      require(mvtnorm)
      
      D <- mat.distances(len)
      # Aquí se escala la matríz por la función exponencial (puede ser otra función!)
      #Este es el efecto del espacio, las celdas vecinas se parecen más entre si
      SIGMA <- sigma^2 * exp(-lambda*D)
      mu <- rep(global.mean, times=len*len)
      # Aquí generamos el contenido de la matriz que después será un objeto ráster
      # El contenido será generado con una distribución normal multivariada
      M <- matrix(nrow=len, ncol=len)
      M[] <- rmvnorm(1, mu, SIGMA)
      return(M)
}

mat.to.rast <- function(mat, ext){
      require(raster)
      
      rast <- raster(mat)
      rast@extent@xmax <- ext
      rast@extent@ymax <- ext
      return(rast)
}


#######Simulating the surfaces

set.seed(34)

#Vamos a trabajar con 5 variables que representarán diversos aspectos del clima. El efecto de la
#autocorrelación espacial en cada variable será diferente y generado de manera aleatoria
#con una distribución uniforme

lambdas.l <- runif(10, 0.01, 0.25) # efectos de autocorrelación
lambdas.l.2 <- runif(10, 0.01, 0.25)
lambdas.s <- runif(10, 0.01, 0.1)
lambdas.s.2 <- runif(10, 0.01, 0.1)
lambdas.s.3 <- runif(10, 0.01, 0.1)
lambdas.s.4 <- runif(10, 0.01, 0.1)

sigma.l <- runif(10, 1, 3)
sigma.l.2 <- runif(10, 1, 3)
sigma.s <- runif(10, 0.2, 1)
sigma.s.2 <- runif(10, 0.2, 1)
sigma.s.3 <- runif(10, 0.2, 1)
sigma.s.4 <- runif(10, 0.2, 1)

large.means <- rnorm(10, 20, 10) #Medias de cada capa
large.means.2 <- rnorm(10, 100, 25)
small.means <- runif(10, -0.1, 0.5)
small.means.2 <- runif(10, -0.1, 0.5)
small.means.3 <- runif(10, -0.1, 0.5)
small.means.4 <- runif(10, -0.1, 0.5)

len <- 100 #Tamaño de las capas ráster
dists <- mat.distances(len = len)

registerDoParallel(cores = 3)

LM.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
      mat <- corr.surf(len = len, global.mean = large.means[i], lambda = lambdas.l[i], sigma = sigma.l[i])
      r <- mat.to.rast(mat, ext = len)
      return(r)
}
gc(reset = T)

LM2.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
   mat <- corr.surf(len = len, global.mean = large.means.2[i], lambda = lambdas.l[i], sigma = sigma.l.2[i])
   r <- mat.to.rast(mat, ext = len)
   return(r)
}
gc(reset = T)

SM.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
   mat <- corr.surf(len = len, global.mean = small.means[i], lambda = lambdas.s[i], sigma = sigma.s[i])
   r <- mat.to.rast(mat, ext = len)
   return(r)
}
gc(reset = T)

SM2.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
   mat <- corr.surf(len = len, global.mean = small.means.2[i], lambda = lambdas.s.2[i], sigma = sigma.s.2[i])
   r <- mat.to.rast(mat, ext = len)
   return(r)
}
gc(reset = T)

SM3.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
   mat <- corr.surf(len = len, global.mean = small.means.3[i], lambda = lambdas.s.3[i], sigma = sigma.s.3[i])
   r <- mat.to.rast(mat, ext = len)
   return(r)
}
gc(reset = T)

SM4.layers <- foreach(i = seq_along(lambdas.l)) %dopar% {
   mat <- corr.surf(len = len, global.mean = small.means.4[i], lambda = lambdas.s.4[i], sigma = sigma.s.4[i])
   r <- mat.to.rast(mat, ext = len)
   return(r)
}
gc(reset = T)

#Normally distributed layers
LM.r <- stack(LM.layers)
names(LM.r) <- paste0("Large-mean-", 1:10)

LM2.r <- stack(LM2.layers)
names(LM2.r) <- paste0("Very-large-mean-", 1:10)

SM.r <- stack(SM.layers)
names(SM.r) <- paste0("Small-mean.1-",1:10)

SM2.r <- stack(SM2.layers)
names(SM2.r) <- paste0("Small-mean.2-",1:10)

SM3.r <- stack(SM3.layers)
names(SM3.r) <- paste0("Small-mean.3-",1:10)

SM4.r <- stack(SM4.layers)
names(SM4.r) <- paste0("Small-mean.4-",1:10)

#Log-normally distributed layers
pos1.r <- exp(SM.r)
names(pos1.r) <- paste0("Log-normal.1-", 1:10)

pos2.r <- exp(SM2.r)
names(pos2.r) <- paste0("Log-normal.2-", 1:10)

pos3.r <- exp(SM3.r)
names(pos3.r) <- paste0("Log-normal.3-", 1:10)

pos4.r <- exp(SM4.r)
names(pos4.r) <- paste0("Log-normal.4-", 1:10)

#Beta-distributed layers
pos1.df <- data.frame(rasterToPoints(pos1.r))
pos2.df <- data.frame(rasterToPoints(pos2.r))
pos3.df <- data.frame(rasterToPoints(pos3.r))
pos4.df <- data.frame(rasterToPoints(pos4.r))


beta.1.df <- foreach(i = 3:12, .combine = cbind) %do% {
   sapply(1:nrow(pos1.df) , function(x){qbeta(p = 0.5, shape1 = pos1.df[x,i], shape2 = pos1.df[x,i] + pos4.df[x, i])})
}

beta.2.df <- foreach(i = 3:12, .combine = cbind) %do% {
   sapply(1:nrow(pos1.df) , function(x){qbeta(p = 0.5, shape1 = pos2.df[x,i], shape2 = pos3.df[x,i])})
}

beta.1.r <- rasterFromXYZ(data.frame(pos1.df[, c("x", "y")], beta.1.df))
names(beta.1.r) <- paste0("Beta.1-", 1:10)

beta.2.r <- rasterFromXYZ(data.frame(pos1.df[, c("x", "y")], beta.2.df))
names(beta.2.r) <- paste0("Beta.2-", 1:10)


#Gamma-distributed layers
LM2.df <- data.frame(rasterToPoints(LM2.r))

gamma.1 <- foreach(i = 3:12, .combine = cbind) %do% {
   sapply(1:nrow(pos1.df) , function(x){qgamma(p = 0.5, shape = beta.2.df[x,i-2], rate = pos2.df[x,i])})
}

gamma.2 <- foreach(i = 3:12, .combine = cbind) %do% {
   sapply(1:nrow(pos1.df) , function(x){qgamma(p = 0.5, shape = pos3.df[x, i-2], scale = LM2.df[x,i])})
}

gamma.1.r <- rasterFromXYZ(data.frame(pos1.df[, c("x", "y")], gamma.1))
names(gamma.1.r) <- paste0("Gamma.1-", 1:10)

gamma.2.r <- rasterFromXYZ(data.frame(pos1.df[, c("x", "y")], gamma.2))
names(gamma.2.r) <- paste0("Gamma.2-", 1:10)

dir.create("Simulated-layers")

all.Layers <- stack(LM.r, LM2.r,
                    SM3.r, SM4.r,
                    pos1.r, pos2.r,
                    beta.1.r, beta.2.r,
                    gamma.1.r, gamma.2.r)

means <- sapply(1:nlayers(all.Layers), function(x){
   d <- density(all.Layers[[x]])
   m <- d$x[which.max(d$y)]
   return(m)
})

variances <- sapply(1:nlayers(all.Layers), function(x){
   require(coda)
   int <- HPDinterval(as.mcmc(all.Layers[[x]][]), 0.69)
   var <- mean((int - means[x])^2)
   return(var)
})
 

saveRDS(all.Layers, "Simulated-layers/All-simulated-layers.rds")

## Layer summaries
Summaries <- data.frame(Layer.name = names(all.Layers),
                            means = means,
                            variance = variances)

write.csv(Summaries, "Simulated-layers/Layer-summaries.csv") 



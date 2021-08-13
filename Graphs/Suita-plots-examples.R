library(mvtnorm)

quants <- expand.grid(x1 = seq(-4, 4, len = 35),
                      x2 = seq(-4, 4, len = 35))

dist <- dmvnorm(x = quants, mean = c(0, 0),
                 sigma = diag(2))

mat <- matrix(dist, 35, 35)

distance <- with(quants, sqrt((0 - x1)^2 + (0 - x2)^2))

dist.mat <- matrix(distance, 35,35)

preds <- with(quants, exp(x1 * 0.001 - 0.3 * x1^2 +  x2 * 0.001 - 0.3 * x2^2))

mat2 <- matrix(preds, 35, 35)

library(plot3D)

pdf("../Graphs/Suitability-examples.pdf", width = 5, height = 5)

persp3D(x = unique(quants$x1),
        y = unique(quants$x2),
        z = mat,theta = 30,phi = 30, col = "lightblue",
      xlab = expression(x[1]), ylab = expression(x[2]), 
      zlab = "Density",
      main = "Density plot",
      border = "grey50", lwd = 0.75)

persp3D(x = unique(quants$x1),
        
        y = unique(quants$x2),
        z = dist.mat,theta = 30,phi = 30, col = "lightblue",
        xlab = expression(x[1]), ylab = expression(x[2]), 
        zlab = "Distance",
        main = "Distance to centroid",
        border = "grey50", lwd = 0.75)

persp3D(x = unique(quants$x1),
        y = unique(quants$x2),
        z = exp(-0.5 * dist.mat),theta = 30,phi = 30, col = "lightblue",
        xlab = expression(x[1]), ylab = expression(x[2]), 
        zlab = "Suitability",
        border = "grey50", lwd = 0.75)

persp3D(x = unique(quants$x1),
      y = unique(quants$x2),
      z = mat2,theta = 30, phi = 30, col = "lightblue",
      xlab = expression(x[1]), ylab = expression(x[2]), 
      zlab = expression(lambda(x[1], x[2])),
      main = "Point intensity",
      border = "grey50", lwd = 0.75)
dev.off()


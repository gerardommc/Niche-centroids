library(ggplot2); library(reshape2)

int.spp.dat <- read.csv("Interacting-spp/Performance-stats.csv")
int.spp.dat <- subset(int.spp.dat, select = c("spp", "method", "rad", "dist.t.cent", "corr.surf" ))

spp.1 <- subset(int.spp.dat, spp == "spp.1")
spp.2 <- subset(int.spp.dat, spp == "spp.2")

rad.sp1 <- spp.1$rad
rad.sp2 <- spp.2$rad

spp.1$rad.op <- rad.sp2
spp.2$rad.op <- rad.sp1

int.spp.dat <- rbind(spp.1, spp.2)
levels(int.spp.dat$method) <- c("Ellipses", "Hardcore", "Poisson")
levels(int.spp.dat$spp) <- c("Species 1", "Species 2")

library(viridis)
p.corr <- ggplot(int.spp.dat) + geom_tile(aes(x = rad, y = rad.op, fill = corr.surf)) + 
            facet_grid(spp ~ method) + 
            labs(fill = expression(rho), x = "Species 1 radius", y = "Species 2 radius") +
            scale_fill_gradientn(colours = magma(25)) +
            theme_dark()

pdf("../Graphs/Spp-radii-corr.pdf", width = 6, height = 4)
p.corr
dev.off()


hard <- subset(int.spp.dat, method == "Hardcore")
pois <- subset(int.spp.dat, method == "Poisson")
ellip <- subset(int.spp.dat, method == "Ellipses")

source("../Random functions/multiplot.R")

p1 <- ggplot(hard) + geom_tile(aes(x = rad, y = rad.op, fill = log10(dist.t.cent))) +
            facet_grid(spp ~ method) +
      labs(fill = expression(log[10]~ S), x = "", y = "Species 2 radius") + 
      scale_fill_gradientn(colours = magma(25)) +
      theme_dark() +
      theme(legend.position = "bottom",
            legend.text = element_text(angle = 45, vjust = 0.5),
            strip.text.y = element_blank())
p2 <- ggplot(pois) + geom_tile(aes(x = rad, y = rad.op, fill = log10(dist.t.cent))) +
      facet_grid(spp ~ method)+
      labs(fill = "", x = "Species 1 radius", y = "") + 
      scale_fill_gradientn(colours = magma(25)) +
      theme_dark() +
      theme(legend.position = "bottom",
            legend.text = element_text(angle = 45, vjust = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            strip.text.y = element_blank())
p3 <- ggplot(ellip) + geom_tile(aes(x = rad, y = rad.op, fill = log10(dist.t.cent))) +
      facet_grid(spp ~ method) +
      labs(fill = "", x = "    ", y = "") + 
      scale_fill_gradientn(colours = magma(25)) +
      theme_dark() +
      theme(legend.position = "bottom",
            legend.text = element_text(angle = 45, vjust = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

multiplot(p1, p2, p3, cols = 3)

pdf("../Graphs/Spp-radii-distance.pdf", width = 6, height = 4.5)
multiplot(p1, p2, p3, cols = 3)
dev.off()

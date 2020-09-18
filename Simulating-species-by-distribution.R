### Script to simulate species controlling the variable types used

## Reading up data
all.layers <- readRDS("Simulated-layers/All-simulated-layers.rds")
l.sum <- read.csv("Simulated-layers/Layer-summaries.csv")
l.sum$Distribution <- c(rep("normal", 40), 
                        rep("log-normal", 20),
                        rep("beta", 20),
                        rep("gamma", 20))

norm <- which(l.sum$Distribution == "normal")
log.norm <- which(l.sum$Distribution == "log-normal")
beta <- which(l.sum$Distribution == "beta")
gamma <- which(l.sum$Distribution == "gamma")

#Setting up single distribution sets
norm.combs <- combn(norm, 3)
l.norm.combs <- combn(log.norm, 3)
beta.combs <- combn(beta, 3)
gamma.combs <- combn(gamma, 3)

norm.samps <- sample(1:ncol(norm.combs), 500, replace = F)
l.norm.samps <- sample(1:ncol(l.norm.combs), 500, replace = F)
beta.samps <- sample(1:ncol(beta.combs), 500, replace = F)
gamma.samps <- sample(1:ncol(gamma.combs), 500, replace = F)

norm.com.samp <- norm.combs[, norm.samps]
l.norm.com.samp <- l.norm.combs[, l.norm.samps]
beta.com.samp <- beta.combs[, beta.samps]
gamma.com.samp <- gamma.combs[, gamma.samps]

##Setting up the variable combinations
var.types <- c(rep("normal", 3), rep("log-normal", 3), rep("gamma", 3), rep("beta", 3))

var.combs <- t(combn(var.types, 3, simplify = T))
var.combs.df <- data.frame(var.combs)
names(var.combs.df) <- c("var1", "var2", "var3")

for(i in 1:3){
      for(j in 1:n){
            
      }
}

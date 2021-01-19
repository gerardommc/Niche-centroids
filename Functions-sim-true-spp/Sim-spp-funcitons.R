temp.optim <- function(x, k = list(k1, k2, k3, k4, k5, k6, k7)){
      y <- with(k, (k1*(x - k2)^k3)/(k4^k3 + (x - k2)^k3) - exp(k7 - (k5 - (x - k2))/(k5 - k6)))
      return(y)
}

rain.resp <- function(x, k = list(k1, k2)){
      y <- with(k, k1 * (1 - exp( - k2 * x)))
      return(y)
}

 
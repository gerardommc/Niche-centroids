imFromStack <- function(r){
      require(spatstat); require(foreach)
      
      r.df <- data.frame(rasterToPoints(r))
      
      ux = sort(unique(r.df$x)) #Extracting unique coordinates
      uy = sort(unique(r.df$y))
      nx = length(ux) #length of unique coordinates
      ny = length(uy)
      ref.cols = match(r.df$x, ux) #position of every data point
      ref.lines = match(r.df$y, uy)
      vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
      ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
      vec[ref.vec] = 1
      data.mask = matrix(vec, max(ref.lines), max(ref.cols))
                         

      im.list <- foreach(i = 3:ncol(r.df)) %do% {
            vec.all = rep(NA, max(ref.lines)*max(ref.cols))
            vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
            vec.all[ref.vec] = r.df[,i]
            lay <- im(matrix(vec.all, max(ref.lines), max(ref.cols),
                             dimnames = list(uy, ux)), xcol = ux, yrow = uy)
            return(lay)
      }
      names(im.list) <- names(r)
      return(im.list)
}

winFromRaster <- function(r){      
      require(spatstat); require(foreach)
      
      r.df <- data.frame(rasterToPoints(r))
      
      ux = sort(unique(r.df$x)) #Extracting unique coordinates
      uy = sort(unique(r.df$y))
      nx = length(ux) #length of unique coordinates
      ny = length(uy)
      ref.cols = match(r.df$x, ux) #position of every data point
      ref.lines = match(r.df$y, uy)
      vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
      ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
      vec[ref.vec] = 1
      data.mask = matrix(vec, max(ref.lines), max(ref.cols))
      w = as.owin(im(data.mask, xcol = ux, yrow = uy)) #Data analysis window
      
      return(w)
}

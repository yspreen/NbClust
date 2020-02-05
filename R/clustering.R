clustering = function(method, game, nc, max_nc, hc, nn) {
  # init values
  cl0 = 0
  cl1 = 0
  cl2 = 0
  clall = 0
  clall1 = 0
  clmax = 0

  if (method == 8) 
  {
    set.seed(1)
    cl2 = kmeans(game,nc+1)$cluster
    set.seed(1)
    clmax = kmeans(game,max_nc)$cluster
    if (nc > 2)
    {
      set.seed(1)
      cl1 = kmeans(game,nc)$cluster
      clall = cbind(cl1, cl2)
      set.seed(1)
      cl0 = kmeans(game,nc-1)$cluster
      clall1 = cbind(cl0, cl1, cl2)
    }
    if (nc == 2)
    {
      set.seed(1)
      cl1 = kmeans(game,nc)$cluster
      clall = cbind(cl1, cl2)
      cl0 = rep(1,nn)
      clall1 = cbind(cl0, cl1, cl2)
    }
    if (nc == 1)
    {
      stop("Number of clusters must be higher than 2")
    }
  } else {
    cl1 = cutree(hc, k=nc)
    cl2 = cutree(hc, k=nc+1)
    clall = cbind(cl1, cl2)
    clmax = cutree(hc, k=max_nc)

    if (nc >= 2)
    {
      cl0 = cutree(hc, k=nc-1)
      clall1 = cbind(cl0, cl1, cl2)
    }
    if (nc == 1)
    {
      cl0 = rep(NA,nn)
      clall1 = cbind(cl0, cl1, cl2)
    }
  }

  return(list(cl0=cl0, cl1=cl1, cl2=cl2, clall=clall, clall1=clall1, clmax=clmax))
}

# mmiller: moved indices definitions to definitions.R
source("definitions.R")

NbClust <-function(data = NULL, diss=NULL, distance ="euclidean", min.nc=2, max.nc=15, method =NULL, index = "all", alphaBeale = 0.1)
{
    timestamp()
    time_start = Sys.time()
    time_last = time_start
    
    x<-0
    min_nc <- min.nc
    max_nc <- max.nc
    
    if(is.null(method))    
      stop("method is NULL")
    method <- pmatch(method, c("ward.D2", "single", "complete", "average", 
                               "mcquitty", "median", "centroid", "kmeans","ward.D"))
    
        
    indice <- pmatch(index, c("kl","ch","hartigan","ccc","scott","marriot","trcovw","tracew","friedman",
                              "rubin","cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball",
                              "ptbiserial","gap", "frey", "mcclain",  "gamma", "gplus", "tau", "dunn", 
                              "hubert", "sdindex", "dindex", "sdbw", "all","alllong"))
    if (is.na(indice))
      stop("invalid clustering index")
    
    if (indice == -1)
      stop("ambiguous index")
    
    if ((indice == 3)|| (indice == 5)|| (indice == 6)|| (indice == 7)|| (indice == 8)|| (indice == 9)|| (indice == 10)|| (indice == 11)|| (indice == 18)|| (indice == 27)|| (indice == 29)|| (indice == 31)|| (indice == 32))
    { 
      if((max.nc-min.nc)<2)
        stop("The difference between the minimum and the maximum number of clusters must be at least equal to 2")
    }
    
    
    
    if(is.null(data))
    {
       
      if(method==8)
      {
        stop("\n","method = kmeans, data matrix is needed")
      }
      else
      {  
        if ((indice == 1 )|| (indice == 2)|| (indice == 3 )|| (indice == 4 )|| (indice == 5)|| (indice == 6)|| (indice == 7)|| (indice == 8)|| (indice == 9)|| (indice == 10)|| (indice == 12)|| (indice == 14)|| (indice == 15)|| (indice == 16)|| (indice == 17)|| (indice == 18)
          || (indice == 19)|| (indice == 20)|| (indice == 23)|| (indice == 24)|| (indice == 25) || (indice == 27)|| (indice == 28)
          || (indice == 29) || (indice == 30)|| (indice == 31) || (indice == 32))
          stop("\n","Data matrix is needed. Only frey, mcclain, cindex, sihouette and dunn can be computed.", "\n")
       
        if(is.null(diss))
          stop("data matrix and dissimilarity matrix are both null")
        else  
          cat("\n","Only frey, mcclain, cindex, sihouette and dunn can be computed. To compute the other indices, data matrix is needed","\n") 
      }
    }
    
    else
    {
      jeu1 <- as.matrix(data)
      numberObsBefore <- dim(jeu1)[1]
      jeu <- na.omit(jeu1) # returns the object with incomplete cases removed 
      nn <- numberObsAfter <- dim(jeu)[1]
      pp <- dim(jeu)[2]    
      TT <- t(jeu)%*%jeu   
      sizeEigenTT <- length(eigen(TT)$value)
      eigenValues <- eigen(TT/(nn-1))$value
      
      # Only for indices using vv : CCC, Scott, marriot, tracecovw, tracew, friedman, rubin
      
      if (any(indice == 4) || (indice == 5) || (indice == 6) || (indice == 7) || (indice == 8) || (indice == 9) || (indice == 10) || (indice == 31) || (indice == 32))
      {
        for (i in 1:sizeEigenTT) 
        {
          if (eigenValues[i] < 0) {
            #cat(paste("There are only", numberObsAfter,"nonmissing observations out of a possible", numberObsBefore ,"observations."))
            stop("The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.")
          } 
        }
        s1 <- sqrt(eigenValues)
        ss <- rep(1,sizeEigenTT)
        for (i in 1:sizeEigenTT) 
        {
          if (s1[i]!=0) 
            ss[i]=s1[i]
        }
        vv <- prod(ss)  
      } 
            
    }
    
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "init", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()
    
    # mmiller: moved to distances.R
    source("distances.R")
    md = init_distances(distance, diss, jeu)
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "distances.R", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()

    # mmiller: moved to methods.R
    source("methods.R")
    init_methods_res = init_methods(max_nc, min_nc, method, md)
    hc = init_methods_res[[1]]
    res = init_methods_res[[2]]
    resCritical = init_methods_res[[3]]
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "methods.R", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()

    # mmiller: moved to main.R
    source("main.R")
    main_res = main(min_nc, max_nc, method, hc, nn, jeu, indice, res, md, TT, ss, vv, pi, pp, resCritical)
    res = main_res[[1]]
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "main.R", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()

    # mmiller: moved to estimate_n.R
    source("estimate_n.R")
    estimate_n_res = estimate_n(indice, min_nc, max_nc)
    ncs_and_indices = estimate_n_res[[1]]
    best.nc = estimate_n_res[[2]]
    attach(ncs_and_indices)
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "estimate_n.R", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()
    
    # mmiller: moved to display_results.R
    source("display_results.R")
    results.final = display_results(res, nc, indice, resCritical, min.nc, max.nc)
    print(sprintf("[ %-*s ] took %.2f seconds", 20, "display_results.R", difftime(Sys.time(), time_last, units = "secs")))
    time_last = Sys.time()
    
    time_end = Sys.time()
    print(sprintf("Total runtime: %.2f seconds", difftime(Sys.time(), time_last, units = "secs")))
    timestamp()
}

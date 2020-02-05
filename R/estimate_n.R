#########################################################################################################
#######################                      Best Number of Clusters                  ###################
#########################################################################################################



estimate_n <- function(indice, min_nc, max_nc) {
  nc.KL<-indice.KL<-0
  if (any(indice == 1) || (indice == 31) || (indice == 32)) 
  {  
    # KL - The value of u, which maximizes KL(u), is regarded as specifying the number of clusters [ClusterSim package].
    nc.KL <- (min_nc:max_nc)[which.max(res[,1])]
    indice.KL <- max(res[,1],na.rm = TRUE)
    best.nc<-nc.KL
  }
  
  nc.CH<-indice.CH<-0
  if (any(indice == 2) || (indice == 31) || (indice == 32)) 
  {
    # CH - The value of u, which maximizes CH(u), is regarded as specifying the number of clusters [ClusterSim package].
    nc.CH <- (min_nc:max_nc)[which.max(res[,2])]
    indice.CH <- max(res[,2],na.rm = TRUE)
    best.nc<-nc.CH
  }
  
  nc.CCC<-indice.CCC<-0
  if (any(indice == 4) || (indice == 31) || (indice == 32))
  {
    # CCC - The maximum value accross the hierarchy levels is used to indicate the optimal number of clusters in data [29].
    nc.CCC <- (min_nc:max_nc)[which.max(res[,4])]
    indice.CCC <- max(res[,4],na.rm = TRUE)
    best.nc<-nc.CCC
  }
  
  nc.DB<-indice.DB<-0 
  if (any(indice == 12) || (indice == 31) || (indice == 32)) 
  {
    # DB - The value of u, which minimizes DB(u), is regarded as specifying the number of clusters [clusterSim package].
    nc.DB <- (min_nc:max_nc)[which.min(res[,12])]
    indice.DB <- min(res[,12],na.rm = TRUE)
    best.nc<-nc.DB
  }
  
  nc.Silhouette<-indice.Silhouette<-0
  if (any(indice == 13) || (indice == 31) || (indice == 32)) 
  {
    # SILHOUETTE - The value of u, which maximizes S(u), is regarded as specifying the number of clusters [ClusterSim package].
    nc.Silhouette <- (min_nc:max_nc)[which.max(res[,13])]
    indice.Silhouette <- max(res[,13],na.rm = TRUE)
    best.nc<-nc.Silhouette
  }
  
  nc.Gap<-indice.Gap<-0
  # GAP - Choose the number of clusters via finding the smallest q such that: Gap(q)=Gap(q+1)-Sq+1 (q=1,\u{85},n-2) [ClusterSim package].
  if (any(indice == 20) || (indice == 32))
  {
    found <- FALSE
    for (ncG in min_nc:max_nc){
      if ((resCritical[ncG-min_nc+1,4] >=0) && (!found)){
        ncGap <- ncG
        indiceGap <- res[ncG-min_nc+1,20]
        found <- TRUE
      }
    }
    if (found){
      nc.Gap <- ncGap
      indice.Gap <- indiceGap
      best.nc<-nc.Gap
    }else{
      nc.Gap <- NA
      indice.Gap <- NA
    }
    
  }
  
  nc.Duda<-indice.Duda<-0
  # DUDA - Choose the number of clusters via finding the smallest q such that: duda >= critical_value [Duda and Hart (1973)].
  
  
  if (any(indice == 14) || (indice == 31) || (indice == 32))
  {
    
    foundDuda <- FALSE
    for (ncD in min_nc:max_nc)
    {
      
      if ((res[ncD-min_nc+1,14]>=resCritical[ncD-min_nc+1,1]) && (!foundDuda))
      {
        ncDuda <- ncD
        indiceDuda <- res[ncD-min_nc+1,14]
        foundDuda <- TRUE
      }
    }
    if (foundDuda)
    {
      nc.Duda <- ncDuda
      indice.Duda <- indiceDuda
      best.nc<-nc.Duda
    }
    else
    {
      nc.Duda <- NA
      indice.Duda <- NA
    }
    
    
  }  
  
  nc.Pseudo<-indice.Pseudo<-0  
  # PSEUDOT2 - Chooses the number of clusters via finding the smallest q such that: pseudot2 <= critical_value [SAS User's guide].
  if (any(indice == 15) || (indice == 31) || (indice == 32))
  {
    
    foundPseudo <- FALSE
    for (ncP in min_nc:max_nc)
    {
      
      if ((res[ncP-min_nc+1,15]<=resCritical[ncP-min_nc+1,2]) && (!foundPseudo))
      {
        ncPseudo <- ncP
        indicePseudo <- res[ncP-min_nc+1,15]
        foundPseudo <- TRUE
      }
    }
    if (foundPseudo)
    {
      nc.Pseudo <- ncPseudo
      indice.Pseudo <- indicePseudo
      best.nc<-nc.Pseudo
    }
    else
    {
      nc.Pseudo <- NA
      indice.Pseudo <- NA
    }
  }
  
  
  nc.Beale<-indice.Beale<-0
  if (any(indice == 16) || (indice == 31) || (indice == 32))
  {
    # BEALE - Chooses the number of clusters via finding the smallest q such that: Fvalue_beale >= 0.1 [Gordon (1999)].
    foundBeale <- FALSE
    for (ncB in min_nc:max_nc)
    {
      
      if ((resCritical[ncB-min_nc+1,3]>=alphaBeale) && (!foundBeale)){
        ncBeale <- ncB
        indiceBeale <- res[ncB-min_nc+1,16]
        foundBeale <- TRUE
      }
    }
    if (foundBeale){
      nc.Beale <- ncBeale
      indice.Beale <- indiceBeale
      best.nc<-nc.Beale
    }
    else
    {
      nc.Beale <- NA
      indice.Beale <- NA
    }
  }
  
  
  nc.ptbiserial<-indice.ptbiserial<-0
  if (any(indice == 19) || (indice == 31) || (indice == 32))
  {
    # POINT-BISERIAL - The maximum value was used to suggest the optimal number of clusters in the data [29].
    nc.ptbiserial <- (min_nc:max_nc)[which.max(res[,19])]
    indice.ptbiserial <- max(res[,19],na.rm = TRUE)
    best.nc<-nc.ptbiserial
  }
  
  foundNC<-foundIndice<-numeric(0)
  nc.Frey<-indice.Frey<-0
  if (any(indice == 21) || (indice == 31) || (indice == 32))
  {
    # FREY AND VAN GROENEWOUD - The best results occured when clustering was continued until the ratio fell below 1.00 for the last 
    #			      series of times. At this point, the cluster level before this series was taken as the optimal partition. 
    #			      If the ration never fell below 1.00, a one cluster solution was assumed [29].
    
    foundFrey <- FALSE
    i<-1
    for (ncF in min_nc:max_nc)
    {          
      
      if (res[ncF-min_nc+1,21] < 1) 
      {
        
        ncFrey <- ncF-1               
        indiceFrey <- res[ncF-1-min_nc+1,21]
        foundFrey <- TRUE
        foundNC[i]<-ncFrey
        foundIndice[i]<-indiceFrey
        i<-i+1
        
      }
      
    }
    if (foundFrey)
    {
      nc.Frey <- foundNC[1]
      indice.Frey <- foundIndice[1]
      best.nc<-nc.Frey
    }
    else 
    {
      nc.Frey <- NA
      indice.Frey <- NA
      print(paste("Frey index : No clustering structure in this data set"))
    }
    
    
  }  
  
  
  nc.McClain<-indice.McClain<-0
  if (any(indice == 22) || (indice == 31) || (indice == 32))
  {
    # MCCLAIN AND RAO - The minimum value of the index was found to give the best recovery information [29].
    nc.McClain <- (min_nc:max_nc)[which.min(res[,22])]
    indice.McClain <- min(res[,22],na.rm = TRUE)
    best.nc<-nc.McClain
    
  }
  
  nc.Gamma<-indice.Gamma<-0
  if (any(indice == 23) || (indice == 31) || (indice == 32))
  {
    # GAMMA - Maximum values were taken to represent the correct hierarchy level [29].
    nc.Gamma <- (min_nc:max_nc)[which.max(res[,23])]
    indice.Gamma <- max(res[,23],na.rm = TRUE)
    best.nc<-nc.Gamma
    
  }
  
  nc.Gplus<-indice.Gplus<-0
  if (any(indice == 24) || (indice == 31) || (indice == 32))
  {
    # GPLUS - Minimum values were used to determine the number of clusters in the data [29].
    nc.Gplus  <- (min_nc:max_nc)[which.min(res[,24])]
    indice.Gplus <- min(res[,24],na.rm = TRUE)
    best.nc<-nc.Gplus
  }
  
  nc.Tau<-indice.Tau<-0
  if (any(indice == 25) || (indice == 31) || (indice == 32))
  {
    # TAU - The maximum value in the hierarchy sequence was taken as indicating the correct number of clusters [29].
    nc.Tau <- (min_nc:max_nc)[which.max(res[,25])]
    indice.Tau <- max(res[,25],na.rm = TRUE)
    best.nc<-nc.Tau
  }
  
  
  #Some indices need to compute difference between hierarchy levels to identify relevant number of clusters
  
  
  if((indice==3)||(indice == 5)||(indice == 6)||(indice == 7)||(indice == 8)||(indice == 9)||(indice == 10)||(indice == 18)||(indice == 27)||(indice == 11)||(indice == 29)||(indice == 31)||(indice == 32))
  {
    
    DiffLev <- array(0, c(max_nc-min_nc+1,12))
    DiffLev[,1] <- min_nc:max_nc
    for (nc3 in min_nc:max_nc)      
    {
      if (nc3==min_nc)
      {    
        DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-NA)   # Hartigan
        DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-NA)   #Scott
        DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA)   # Marriot
        DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-NA)   #Trcovw
        DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA)   #Tracew
        DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-NA)   #Friedman
        DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA)  #Rubin
        DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-NA)  # Ball
        DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA) # Hubert   
        DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index
        
        
      }
      else
      {  
        if(nc3==max_nc)
        { 
          DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3])
          DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5])
          DiffLev[nc3-min_nc+1,4] <- abs(res[nc3-min_nc+1,6]-NA) # Marriot
          DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])  # trcovw
          DiffLev[nc3-min_nc+1,6] <- abs(res[nc3-min_nc+1,8]-NA) #traceW            
          DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
          DiffLev[nc3-min_nc+1,8] <- abs(res[nc3-min_nc+1,10]-NA) #Rubin
          DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])
          DiffLev[nc3-min_nc+1,10] <- abs(res[nc3-min_nc+1,27]-NA)
          DiffLev[nc3-min_nc+1,12] <- abs(res[nc3-min_nc+1,29]-NA) # D index  
          
          
        }
        else      
        {
          
          DiffLev[nc3-min_nc+1,2] <- abs(res[nc3-min_nc+1,3]-res[nc3-min_nc,3]) # Hartigan              
          DiffLev[nc3-min_nc+1,3] <- abs(res[nc3-min_nc+1,5]-res[nc3-min_nc,5]) 
          DiffLev[nc3-min_nc+1,4] <- ((res[nc3-min_nc+2,6]-res[nc3-min_nc+1,6])-(res[nc3-min_nc+1,6]-res[nc3-min_nc,6]))
          DiffLev[nc3-min_nc+1,5] <- abs(res[nc3-min_nc+1,7]-res[nc3-min_nc,7])
          DiffLev[nc3-min_nc+1,6] <- ((res[nc3-min_nc+2,8]-res[nc3-min_nc+1,8])-(res[nc3-min_nc+1,8]-res[nc3-min_nc,8]))
          DiffLev[nc3-min_nc+1,7] <- abs(res[nc3-min_nc+1,9]-res[nc3-min_nc,9])
          DiffLev[nc3-min_nc+1,8] <- ((res[nc3-min_nc+2,10]-res[nc3-min_nc+1,10])-(res[nc3-min_nc+1,10]-res[nc3-min_nc,10]))
          DiffLev[nc3-min_nc+1,9] <- abs(res[nc3-min_nc+1,18]-res[nc3-min_nc,18])  
          DiffLev[nc3-min_nc+1,10] <- abs((res[nc3-min_nc+1,27]-res[nc3-min_nc,27]))             
          DiffLev[nc3-min_nc+1,12] <-((res[nc3-min_nc+2,29]-res[nc3-min_nc+1,29])-(res[nc3-min_nc+1,29]-res[nc3-min_nc,29])) #Dindex     
          
        }         
      }         
    }
  }
  
  nc.Hartigan<-indice.Hartigan<-0
  if (any(indice == 3) || (indice == 31) || (indice == 32))
  {
    # HARTIGAN - The maximum differences between hierarchy levels were taken as indicating the correct number of clusters in the data [29].
    nc.Hartigan <- DiffLev[,1][which.max(DiffLev[,2])]
    indice.Hartigan <- max(DiffLev[,2],na.rm = TRUE)
    best.nc<-nc.Hartigan
  }
  
  nc.Ratkowsky<-indice.Ratkowsky<-0
  if (any(indice == 17) || (indice == 31) || (indice == 32))
  {
    # RATKOWSKY - The optimal number of groups is taken as the level where this criterion exhibits its maximum value [29].
    nc.Ratkowsky <- (min_nc:max_nc)[which.max(res[,17])]
    indice.Ratkowsky <- max(res[,17],na.rm = TRUE)
    best.nc<-nc.Ratkowsky
  }
  
  nc.cindex<-indice.cindex<-0
  if (any(indice == 11) || (indice == 31) || (indice == 32)) 
  {
    #CINDEX - The minimum value across the hierarchy levels was used to indicate the optimal number of clusters [29].
    nc.cindex <- (min_nc:max_nc)[which.min(res[,11])]
    indice.cindex <- min(res[,11],na.rm = TRUE)
    best.nc<-nc.cindex
  }  
  
  nc.Scott<-indice.Scott<-0
  if (any(indice == 5) || (indice == 31) || (indice == 32))
  {
    # SCOTT - The maximum difference between hierarchy levels was used to suggest the correct number of partitions [29].
    nc.Scott <- DiffLev[,1][which.max(DiffLev[,3])]
    indice.Scott <- max(DiffLev[,3],na.rm = TRUE)
    best.nc<-nc.Scott
  }
  
  nc.Marriot<-indice.Marriot<-0
  if (any(indice == 6) || (indice == 31) || (indice == 32))
  {
    # MARRIOT - The maximum difference between successive levels was used to determine the best partition level [29].
    nc.Marriot <- DiffLev[,1][which.max(DiffLev[,4])]
    round(nc.Marriot, digits=1)
    indice.Marriot <- max(DiffLev[,4],na.rm = TRUE)
    best.nc<-nc.Marriot
  }
  
  nc.TrCovW<-indice.TrCovW<-0
  if (any(indice == 7) || (indice == 31) || (indice == 32))
  {
    nc.TrCovW <- DiffLev[,1][which.max(DiffLev[,5])]
    indice.TrCovW <- max(DiffLev[,5],na.rm = TRUE)
    best.nc<-nc.TrCovW
  }
  
  
  nc.TraceW<-indice.TraceW<-0
  if (any(indice == 8) || (indice == 31) || (indice == 32))
  {
    # TRACE W - To determine the number of clusters in the data, maximum difference scores were used [29].
    nc.TraceW <- DiffLev[,1][which.max(DiffLev[,6])]
    indice.TraceW <- max(DiffLev[,6],na.rm = TRUE)
    best.nc<-nc.TraceW
  }
  
  nc.Friedman<-indice.Friedman<-0
  if (any(indice == 9) || (indice == 31) || (indice == 32))
  {
    # FRIEDMAN - The maximum difference in values of trace W-1B criterion was used to indicate the optimal number of clusters [29].
    nc.Friedman <- DiffLev[,1][which.max(DiffLev[,7])]
    indice.Friedman <- max(DiffLev[,7],na.rm = TRUE)
    best.nc<-nc.Friedman
  }
  
  nc.Rubin<-indice.Rubin<-0
  if (any(indice == 10) || (indice == 31) || (indice == 32))
  {
    # RUBIN - The difference between levels was used [29].
    nc.Rubin <- DiffLev[,1][which.min(DiffLev[,8])]
    indice.Rubin <- min(DiffLev[,8],na.rm = TRUE)
    best.nc<-nc.Rubin
  }
  
  nc.Ball<-indice.Ball<-0
  if (any(indice == 18) || (indice == 31) || (indice == 32))
  {
    # BALL - The largest difference between levels was used to indicate the optimal solution [29].
    nc.Ball <- DiffLev[,1][which.max(DiffLev[,9])]
    indice.Ball <- max(DiffLev[,9],na.rm = TRUE)
    best.nc<-nc.Ball
  }
  
  
  nc.Dunn<-indice.Dunn<-0
  if (any(indice == 26) || (indice == 31) || (indice == 32)) 
  {
    # Dunn - 
    nc.Dunn <- (min_nc:max_nc)[which.max(res[,26])]
    indice.Dunn <- max(res[,26],na.rm = TRUE)
    best.nc<-nc.Dunn
  }
  
  
  nc.Hubert<-indice.Hubert<-0
  if (any(indice == 27) || (indice == 31) || (indice == 32)) 
  {       
    # Hubert - 
    nc.Hubert  <- 0.00
    indice.Hubert  <- 0.00
    #x11()
    par(mfrow = c(1,2))
    plot(x_axis,res[,27], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert Statistic values")))
    plot(DiffLev[,1],DiffLev[,10], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Hubert statistic second differences")))
    cat(paste ("*** : The Hubert index is a graphical method of determining the number of clusters.
                  In the plot of Hubert index, we seek a significant knee that corresponds to a 
                  significant increase of the value of the measure i.e the significant peak in Hubert
                  index second differences plot.", "\n", "\n"))
  }
  
  nc.sdindex<-indice.sdindex<-0
  if (any(indice == 28) || (indice == 31) || (indice == 32)) 
  {
    # SD - 
    nc.sdindex <- (min_nc:max_nc)[which.min(res[,28])]
    indice.sdindex<- min(res[,28],na.rm = TRUE)
    best.nc<-nc.sdindex
  }
  
  
  nc.Dindex<-indice.Dindex<-0
  if (any(indice == 29) || (indice == 31) || (indice == 32)) 
  {
    
    nc.Dindex <- 0.00
    indice.Dindex<- 0.00
    #x11()
    par(mfrow = c(1,2))
    plot(x_axis,res[,29], tck=0, type="b", col="red", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Dindex Values")))
    plot(DiffLev[,1],DiffLev[,12], tck=0, type="b", col="blue", xlab= expression(paste("Number of clusters ")), ylab= expression(paste("Second differences Dindex Values")))
    cat(paste ("*** : The D index is a graphical method of determining the number of clusters. 
                  In the plot of D index, we seek a significant knee (the significant peak in Dindex
                  second differences plot) that corresponds to a significant increase of the value of
                  the measure.", "\n", "\n"))
  }
  
  nc.SDbw<-indice.SDbw<-0
  if (any(indice == 30) || (indice == 31) || (indice == 32)) 
  {
    # SDbw - 
    nc.SDbw <- (min_nc:max_nc)[which.min(res[,30])]
    indice.SDbw<- min(res[,30],na.rm = TRUE)  
    best.nc<-nc.SDbw
  }
  
  
  ncs = list(
    nc.KL = nc.KL,
    indice.KL = indice.KL,
    nc.CH = nc.CH,
    indice.CH = indice.CH,
    nc.Hartigan = nc.Hartigan,
    indice.Hartigan = indice.Hartigan,
    nc.CCC = nc.CCC,
    indice.CCC = indice.CCC,
    nc.Scott = nc.Scott,
    indice.Scott = indice.Scott,
    nc.Marriot = nc.Marriot,
    indice.Marriot = indice.Marriot,
    nc.TrCovW = nc.TrCovW,
    indice.TrCovW = indice.TrCovW,
    nc.TraceW = nc.TraceW,
    indice.TraceW = indice.TraceW,
    nc.Friedman = nc.Friedman,
    indice.Friedman = indice.Friedman,
    nc.Rubin = nc.Rubin,
    indice.Rubin = indice.Rubin,
    nc.cindex = nc.cindex,
    indice.cindex = indice.cindex,
    nc.DB = nc.DB,
    indice.DB = indice.DB,
    nc.Silhouette = nc.Silhouette,
    indice.Silhouette = indice.Silhouette,
    nc.Duda = nc.Duda,
    indice.Duda = indice.Duda,
    nc.Pseudo = nc.Pseudo,
    indice.Pseudo = indice.Pseudo,
    nc.Beale = nc.Beale,
    indice.Beale = indice.Beale,
    nc.Ratkowsky = nc.Ratkowsky,
    indice.Ratkowsky = indice.Ratkowsky,
    nc.Ball = nc.Ball,
    indice.Ball = indice.Ball,
    nc.ptbiserial = nc.ptbiserial,
    indice.ptbiserial = indice.ptbiserial,
    nc.Gap = nc.Gap,
    indice.Gap = indice.Gap,
    nc.Frey = nc.Frey,
    indice.Frey = indice.Frey,
    nc.McClain = nc.McClain,
    indice.McClain = indice.McClain,
    nc.Gamma = nc.Gamma,
    indice.Gamma = indice.Gamma,
    nc.Gplus = nc.Gplus,
    indice.Gplus = indice.Gplus,
    nc.Tau = nc.Tau,
    indice.Tau = indice.Tau,
    nc.Dunn = nc.Dunn,
    indice.Dunn = indice.Dunn,
    nc.Hubert = nc.Hubert,
    indice.Hubert = indice.Hubert,
    nc.sdindex = nc.sdindex,
    indice.sdindex = indice.sdindex,
    nc.Dindex = nc.Dindex,
    indice.Dindex = indice.Dindex,
    nc.SDbw = nc.SDbw,
    indice.SDbw = indice.SDbw
  )
  
  return(list(ncs, best.nc))
}
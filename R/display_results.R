######################################################################################################################
########################                Displaying results             #########################################
######################################################################################################################

display_results <- function(res, nc, indice, resCritical, min.nc, max.nc) {
  if (indice < 31)
  {
    res <- res[,c(indice)]
    
    if (indice == 14) { resCritical <- resCritical[,1]  }
    if (indice == 15) { resCritical <- resCritical[,2] }
    if (indice == 16) { resCritical <- resCritical[,3] }
    if (indice == 20) { resCritical <- resCritical[,4] }        
    
  }
  
  if (indice == 31) 
  { 
    res <- res[,c(1:19,21:22,26:30)]
    resCritical <- resCritical[,c(1:3)]        
    
  }
  
  
  
  if (any(indice == 20) || (indice == 23) || (indice == 24) || (indice == 25) || (indice == 32))
  {
    
    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
                 nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman, 
                 indice.Friedman, nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette,
                 indice.Silhouette, nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky,
                 indice.Ratkowsky, nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Gap, indice.Gap, 
                 nc.Frey, indice.Frey, nc.McClain, indice.McClain, nc.Gamma, indice.Gamma, nc.Gplus, indice.Gplus,
                 nc.Tau, indice.Tau, nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw)
    results1<-matrix(c(results),nrow=2,ncol=30)
    resultats <- matrix(c(results),nrow=2,ncol=30,dimnames = list(c("Number_clusters","Value_Index"),
                                                                  c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
                                                                    "TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
                                                                    "Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
                                                                    "Gap", "Frey", "McClain", "Gamma", "Gplus", "Tau", "Dunn", 
                                                                    "Hubert", "SDindex", "Dindex", "SDbw")))
  }
  else
  {
    
    results <- c(nc.KL, indice.KL, nc.CH, indice.CH, nc.Hartigan, indice.Hartigan, nc.CCC, indice.CCC, nc.Scott, indice.Scott,
                 nc.Marriot, indice.Marriot, nc.TrCovW, indice.TrCovW, nc.TraceW, indice.TraceW, nc.Friedman, indice.Friedman, 
                 nc.Rubin, indice.Rubin, nc.cindex, indice.cindex, nc.DB, indice.DB, nc.Silhouette, indice.Silhouette,
                 nc.Duda, indice.Duda, nc.Pseudo, indice.Pseudo, nc.Beale, indice.Beale, nc.Ratkowsky, indice.Ratkowsky,
                 nc.Ball, indice.Ball, nc.ptbiserial, indice.ptbiserial, nc.Frey, indice.Frey, nc.McClain, indice.McClain, 
                 nc.Dunn, indice.Dunn, nc.Hubert, indice.Hubert, nc.sdindex, indice.sdindex, nc.Dindex, indice.Dindex, nc.SDbw, indice.SDbw 
    )
    results1<-matrix(c(results),nrow=2,ncol=26)
    resultats <- matrix(c(results),nrow=2,ncol=26,dimnames = list(c("Number_clusters","Value_Index"),
                                                                  c("KL","CH","Hartigan","CCC", "Scott", "Marriot", "TrCovW",
                                                                    "TraceW", "Friedman", "Rubin", "Cindex", "DB", "Silhouette",
                                                                    "Duda","PseudoT2", "Beale", "Ratkowsky", "Ball", "PtBiserial",
                                                                    "Frey", "McClain", "Dunn", 		"Hubert", "SDindex", "Dindex", "SDbw")))
    
  }
  
  
  if (any(indice <= 20)||(indice == 23)||(indice == 24)||(indice == 25)) 
  {   
    resultats <- resultats[,c(indice)]   
  }
  
  if (any(indice == 21)|| (indice == 22)) 
  {  
    indice3 <-indice-1
    resultats <- resultats[,c(indice3)]   
  }
  
  if (any(indice == 26) || (indice == 27) || (indice == 28) || ( indice == 29)|| ( indice == 30)) 
  { 
    indice4 <- indice-4     
    resultats <- resultats[,c(indice4)] 
  }
  
  
  resultats<-round(resultats, digits=4)
  res<-round(res, digits=4)
  resCritical<-round(resCritical, digits=4)
  
  #  if (numberObsAfter != numberObsBefore) 
  #  {
  #     cat(paste(numberObsAfter,"observations were used out of", numberObsBefore ,"possible observations due to missing values."))
  #  }
  
  #  if (numberObsAfter == numberObsBefore) 
  #  {
  #     cat(paste("All", numberObsAfter,"observations were used.", "\n", "\n"))
  #  }
  
  
  
  ######################## Summary results #####################################
  
  
  if(any(indice == 31) || (indice == 32))
  {
    cat("*******************************************************************", "\n")
    cat("* Among all indices:                                               ", "\n")
    BestCluster<-results1[1,]
    c=0
    for(i in min.nc:max.nc)
    {
      vect<-which(BestCluster==i)
      if(length(vect)>0)
        cat("*",length(vect), "proposed", i,"as the best number of clusters", "\n")
      
      if(c<length(vect))
      { 
        j=i 
        c<-length(vect)
      }
    }
    
    cat("\n","                  ***** Conclusion *****                           ", "\n", "\n")
    cat("* According to the majority rule, the best number of clusters is ",j , "\n", "\n", "\n")
    cat("*******************************************************************", "\n")
    
    
    ########################## The Best partition    ###################
    
    if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
        (method == 5) || (method == 6) || (method == 7)||(method == 9))         
      partition<- cutree(hc, k=j)
    
    else
    {
      set.seed(1)
      partition<-kmeans(jeu,j)$cluster
    }
    
  }
  
  
  if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
      ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)||(indice==14)
      ||(indice==15)||(indice==16)||(indice==17)||(indice==18)||(indice==19)||(indice==20)
      ||(indice==21)||(indice==22)||(indice==23)||(indice==24)||(indice==25)||(indice==26)
      ||(indice==28)||(indice==30))
  {
    if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
        (method == 5) || (method == 6) || (method == 7) || (method == 9)) 
      
      partition<- cutree(hc, k=best.nc)
    
    else
    {
      set.seed(1)
      partition<-kmeans(jeu,best.nc)$cluster
    }
    
  }
  
  
  #########################  Summary results   ############################
  
  
  
  if ((indice == 14)|| (indice == 15)|| (indice == 16)|| (indice == 20)|| (indice == 31)|| (indice == 32))
  { 
    results.final <- list(All.index=res,All.CriticalValues=resCritical,Best.nc=resultats, Best.partition=partition)
  }
  
  if ((indice == 27)|| (indice == 29))
    results.final <- list(All.index=res)
  
  if (any(indice==1)||(indice==2)||(indice==3)||(indice==4)||(indice==5)||(indice==6)||(indice==7)
      ||(indice==8)||(indice==9)||(indice==10)||(indice==11)||(indice==12)||(indice==13)
      ||(indice==17)||(indice==18)||(indice==19)||(indice==21)||(indice==22)||(indice==23)||(indice==24)
      ||(indice==25)||(indice==26)||(indice==28)||(indice==30))  
    
    results.final <- list(All.index=res,Best.nc=resultats, Best.partition=partition)
  
  
  
  return(results.final)
}
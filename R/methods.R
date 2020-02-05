#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                                                                      #
#                                              Methods                                                                 #
#                                                                                                                      #
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


init_methods <- function(max_nc, min_nc, method, md) {
  res <- array(0, c(max_nc-min_nc+1,30))
  x_axis <- min_nc:max_nc
  resCritical <- array(0, c(max_nc-min_nc+1,4))
  rownames(resCritical) <- min_nc:max_nc
  colnames(resCritical) <- c("CritValue_Duda", "CritValue_PseudoT2", "Fvalue_Beale", "CritValue_Gap")
  rownames(res) <- min_nc:max_nc
  colnames(res) <- c("KL","CH","Hartigan","CCC","Scott","Marriot", "TrCovW", "TraceW","Friedman","Rubin","Cindex","DB",
                     "Silhouette", "Duda", "Pseudot2", "Beale", "Ratkowsky", "Ball", "Ptbiserial", "Gap", "Frey", "McClain","Gamma", "Gplus", "Tau", "Dunn", 
                     "Hubert", "SDindex", "Dindex", "SDbw")   
  
  if (is.na(method))
    stop("invalid clustering method")
  if (method == -1)
    stop("ambiguous method")
  if (method == 1) 
  {
    hc<-hclust(md,method = "ward.D2")      
  }
  if (method == 2) 
  {
    hc<-hclust(md,method = "single")		
  }
  if (method == 3)
  {
    hc<-hclust(md,method = "complete")		
  }
  
  if (method == 4) 
  {
    hc<-hclust(md,method = "average")
  }
  
  if (method == 5) 
  {
    hc<-hclust(md,method = "mcquitty")
    
  }
  if (method == 6) 
  {
    hc<-hclust(md,method = "median")
    
  }
  if (method == 7) 
  {
    hc<-hclust(md,method = "centroid")
    
  }
  if (method == 9) 
  {
    hc<-hclust(md,method = "ward.D")
    
  }
  
  return(list(hc, res, resCritical))
}
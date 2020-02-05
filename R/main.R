main = function(min_nc, max_nc, method, hc, nn, jeu, indice, res, md, TT, ss, vv, pi, pp, resCritical) {
  for (nc in min_nc:max_nc)
  {  
    
    if (any(method == 1) || (method == 2) || (method == 3) || (method == 4) || 
        (method == 5) || (method == 6) || (method == 7)||(method == 9)) 
    {
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
    
    if (method == 8) 
    {
      set.seed(1)
      cl2 = kmeans(jeu,nc+1)$cluster
      set.seed(1)
      clmax = kmeans(jeu,max_nc)$cluster
      if (nc > 2)
      {
        set.seed(1)
        cl1 = kmeans(jeu,nc)$cluster
        clall = cbind(cl1, cl2)
        set.seed(1)
        cl0 = kmeans(jeu,nc-1)$cluster
        clall1 = cbind(cl0, cl1, cl2)
      }
      if (nc == 2)
      {
        set.seed(1)
        cl1 = kmeans(jeu,nc)$cluster
        clall = cbind(cl1, cl2)
        cl0 = rep(1,nn)
        clall1 = cbind(cl0, cl1, cl2)
      }
      if (nc == 1)
      {
        stop("Number of clusters must be higher than 2")
      }
      
    }
    
    j = table(cl1)  # table uses the cross-classifying factors to build a contingency table of the counts at each combination of factor levels.
    s = sum(j==1)    
    j2 = table(cl2)
    s2 = sum(j2==1)
    
    
    ########### Indices.Traces-hartigan - 3e colonne de res ############ 
    if (any(indice == 3) || (indice == 31) || (indice == 32))
    {    
      res[nc-min_nc+1,3] = Indices.Traces(jeu, md, clall1, index = "hart")	
    } 
    
    ########### Cubic Clustering Criterion-CCC  - 4e colonne de res ############
    if (any(indice == 4) || (indice == 31) || (indice == 32))
    {    	  
      res[nc-min_nc+1,4] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$ccc
    }
    
    ########### Scott and Symons - 5e colonne de res ############
    if (any(indice == 5) || (indice == 31) || (indice == 32))
    {     	  
      res[nc-min_nc+1,5] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$scott
    }
    
    ########### Marriot - 6e colonne de res ############
    if (any(indice == 6) || (indice == 31) || (indice == 32))
    {   	  
      res[nc-min_nc+1,6] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$marriot
    }	
    
    ########### Trace Cov W - 7e colonne de res ############
    if (any(indice == 7) || (indice == 31) || (indice == 32))
    {   	 
      res[nc-min_nc+1,7] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$trcovw	  
    }
    
    ########### Trace W - 8e colonne de res ############
    if (any(indice == 8) || (indice == 31) || (indice == 32))
    {  	  
      res[nc-min_nc+1,8] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$tracew
    }
    
    ########### Friedman - 9e colonne de res ############
    if (any(indice == 9) || (indice == 31) || (indice == 32))
    {     	  
      res[nc-min_nc+1,9] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$friedman
    }
    
    ########### Rubin - 10e colonne de res ############
    if (any(indice == 10) || (indice == 31) || (indice == 32))
    {     	  
      res[nc-min_nc+1,10] = Indices.WBT(x=jeu, cl=cl1, P=TT,s=ss,vv=vv)$rubin
    }
    
    
    ########### Indices.WKWL-duda - 14e colonne de res ############
    if (any(indice == 14) || (indice == 31) || (indice == 32))
    {  
      res[nc-min_nc+1,14] = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$duda	
    }
    
    
    ########### Indices.WKWL-pseudot2 - 15e colonne de res ############
    if (any(indice == 15) || (indice == 31) || (indice == 32))
    {   	  
      res[nc-min_nc+1,15] = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$pseudot2	
    }
    
    ########### Indices.WKWL-beale - 16e colonne de res ############
    if (any(indice == 16) || (indice == 31) || (indice == 32))
    {             	  
      res[nc-min_nc+1,16] = beale = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$beale
    }
    
    ########### Indices.WKWL- duda or pseudot2 or beale ############   
    if (any(indice == 14) || (indice == 15) || (indice == 16) || (indice == 31) || (indice == 32))
    {      	  
      NM = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NM
      NK = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NK
      NL = Indices.WKWL(x=jeu,cl1=cl1,cl2=cl2)$NL
      zz = 3.20 # Best standard score in Milligan and Cooper 1985
      zzz = zz*sqrt(2*(1-8/((pi^2)*pp))/(NM*pp))
      
      
      
      if (any(indice == 14) || (indice == 31) || (indice == 32))
      {
        resCritical[nc-min_nc+1,1] = critValue = 1-(2/(pi*pp))-zzz
      }
      
      if ((indice == 15)|| (indice == 31) || (indice == 32))
      {
        critValue = 1-(2/(pi*pp))-zzz
        resCritical[nc-min_nc+1,2] = ((1-critValue)/critValue)*(NK+NL-2)
        
      }
      
      
      if (any(indice == 16) || (indice == 31) || (indice == 32))
      {
        df2 = (NM-2)*pp
        resCritical[nc-min_nc+1,3] = 1-pf(beale,pp,df2)
      }
    }
    
    ########### Indices.TracesL-ball - 18e colonne de res ############
    if (any(indice == 18) || (indice == 31) || (indice == 32))
    {        	  
      res[nc-min_nc+1,18] = Indices.Traces(jeu, md, clall1, index = "ball")
    }
    
    ########### Indice.Point-Biserial - 19e colonne de res ############ 
    if (any(indice == 19) || (indice == 31) || (indice == 32))
    {       
      res[nc-min_nc+1,19] = Indice.ptbiserial(x=jeu, md=md, cl1=cl1)     
    }
    
    ########### Gap - 20e colonne de res ############       	
    if (any(indice == 20) || (indice == 32))
    {
      
      if (method == 1) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "ward.D2", d = NULL, centrotypes = "centroids")
      }
      if (method == 2) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "single", d = NULL, centrotypes = "centroids")
      }
      if (method == 3) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "complete", d = NULL, centrotypes = "centroids")
      }
      if (method == 4) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "average", d = NULL, centrotypes = "centroids")
      }
      if (method == 5) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "mcquitty", d = NULL, centrotypes = "centroids")
      }
      if (method == 6) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "median", d = NULL, centrotypes = "centroids")
      }
      if (method == 7) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "centroid", d = NULL, centrotypes = "centroids")
      }
      if (method == 9) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "ward.D", d = NULL, centrotypes = "centroids")
      }
      if (method == 8) {
        resultSGAP = Indice.Gap(x=jeu, clall=clall, reference.distribution = "unif", B = 10, method = "k-means", d = NULL, centrotypes = "centroids")
      }
      res[nc-min_nc+1,20] = resultSGAP$gap
      resCritical[nc-min_nc+1,4] = resultSGAP$diffu
    }
    
    if (nc >=2)
    {
      ########### Indices.Traces-kl - 1e colonne de res ############ 	  
      if (any(indice == 1) || (indice == 31) || (indice == 32)) 
      {	 
        res[nc-min_nc+1,1] = Indices.Traces(jeu, md, clall1, index = "kl")
      }
      
      ########### Indices.Traces-ch - 2e colonne de res ############
      if (any(indice == 2) || (indice == 31) || (indice == 32)) 
      {		   
        res[nc-min_nc+1,2] = Indices.Traces(jeu, md, clall1, index = "ch")
      }
      
      ########### Indice.cindex - 11e colonne de res ############
      if (any(indice == 11) || (indice == 31) || (indice == 32)) 
      {	  	   
        res[nc-min_nc+1,11] = Indice.cindex(d=md, cl=cl1)
      }
      
      ########### Indice.DB  - 12e colonne de res ############
      if (any(indice == 12) || (indice == 31) || (indice == 32)) 
      {		   
        res[nc-min_nc+1,12] = Indice.DB(x=jeu, cl=cl1, d = NULL, centrotypes = "centroids", p = 2, q = 2)$DB
      }                         
      
      ########### Silhouette - 13e colonne de res ############
      if (any(indice == 13) || (indice == 31) || (indice == 32)) 
      {		   
        res[nc-min_nc+1,13] = Indice.S(d=md, cl=cl1)
      }
      
      ########### Indices.Traces-ratkowsky- 17e colonne de res ############
      if (any(indice == 17) || (indice == 31) || (indice == 32))
      {  	  
        res[nc-min_nc+1,17] = Indices.Traces(jeu, md, clall1, index = "ratkowsky")
      }
      
      ########### Indice.Frey - 21e colonne de res ############
      if (any(indice == 21) || (indice == 31) || (indice == 32))
      {      
        res[nc-min_nc+1,21] = Index.15and28(cl1=cl1,cl2=cl2,md=md)$frey
      }
      
      ########### Indice.McClain - 22e colonne de res ############
      if (any(indice == 22) || (indice == 31) || (indice == 32))
      {  	     
        res[nc-min_nc+1,22] = Index.15and28(cl1=cl1,cl2=cl2,md=md)$mcclain
      }
      
      ########### Indice.Gamma - 23e colonne de res ############ 
      if (any(indice == 23) || (indice == 32))
      {            
        
        res[nc-min_nc+1,23] = Index.sPlussMoins(cl1=cl1,md=md)$gamma
      }
      
      ########### Indice.Gplus- 24e colonne de res ############
      if (any(indice == 24) || (indice == 32))
      {    	     
        res[nc-min_nc+1,24] = Index.sPlussMoins(cl1=cl1,md=md)$gplus
      }
      
      ########### Indice.Tau  - 25e colonne de res ############
      if (any(indice == 25) || (indice == 32))
      {   	     
        res[nc-min_nc+1,25] = Index.sPlussMoins(cl1=cl1,md=md)$tau
      } 
      
      ########### Indices.Dunn  - 26e colonne de res ############
      if (any(indice == 26 ) || (indice == 31) || (indice == 32))
      {    	    
        res[nc-min_nc+1,26] = Index.dunn(md, cl1, Data=jeu, method=NULL)	
      } 
      
      ########### Indices.Hubert - 27e colonne de res ############
      if (any(indice == 27 ) || (indice == 31) || (indice == 32))
      {         	     
        res[nc-min_nc+1,27] = Index.Hubert(jeu, cl1)	
      }	 
      
      ########### Indices.SD - 28e colonne de res ############
      if (any(indice == 28 ) || (indice == 31) || (indice == 32))
      {	    
        res[nc-min_nc+1,28] = Index.sdindex(jeu, clmax, cl1)
      }	
      
      ########### Indices.Dindex - 29e colonne de res ############ 
      if (any(indice == 29 ) || (indice == 31) || (indice == 32))
      {        	        
        res[nc-min_nc+1,29] = Index.Dindex(cl1, jeu)	
      }  
      
      ########### Indices.SDbw - 30e colonne de res ############
      if (any(indice == 30 ) || (indice == 31) || (indice == 32))
      { 	      
        res[nc-min_nc+1,30] = Index.SDbw(jeu, cl1)	
      }      	
      
    }
    
    else
    {
      res[nc-min_nc+1,1] = NA
      res[nc-min_nc+1,2] = NA
      res[nc-min_nc+1,11] = NA
      res[nc-min_nc+1,12] = NA
      res[nc-min_nc+1,13] = NA
      res[nc-min_nc+1,17] = NA
      res[nc-min_nc+1,21] = NA
      res[nc-min_nc+1,22] = NA
      res[nc-min_nc+1,23] = NA
      res[nc-min_nc+1,24] = NA
      res[nc-min_nc+1,25] = NA
      res[nc-min_nc+1,26] = NA
      res[nc-min_nc+1,27] = NA
      res[nc-min_nc+1,28] = NA
      res[nc-min_nc+1,29] = NA
      res[nc-min_nc+1,30] = NA
    }
  }
  return(list(res=res))
}
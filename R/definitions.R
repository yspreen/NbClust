#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                                                                      #
#                                              Indices                                                                 #
#                                                                                                                      #
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#




##############################
#                            #
#        SD and SDbw         #
#                            #
##############################    


centers = function(cl,x)
{
  x = as.matrix(x)
  n = length(cl)
  k = max(cl)
  centers = matrix(nrow = k, ncol = ncol(x))
  {
    for (i in 1:k) 
    {
      for (j in 1:ncol(x)) 
      {
        centers[i, j] = mean(x[cl == i, j])
      }
    }
  }
  return(centers)
}    

Average.scattering = function (cl, x)
{
  x = as.matrix(x)
  n = length(cl)
  k = max(cl)
  centers.matrix = centers(cl,x)
  
  cluster.size = numeric(0)  
  variance.clusters = matrix(0, ncol = ncol(x), nrow = k)
  var = matrix(0, ncol = ncol(x), nrow = k)
  
  for (u in 1:k) 
    cluster.size[u] = sum(cl == u)
  
  for (u in 1:k) 
  {  
    for (j in 1:ncol(x)) 
    { 
      for(i in 1:n) 
      {     				   
        if(cl[i]==u)                   
          variance.clusters[u,j] = variance.clusters[u,j]+(x[i, j]-centers.matrix[u,j])^2 
      }
    }            
  }
  
  for (u in 1:k) 
  {    
    for (j in 1:ncol(x)) 
      variance.clusters[u,j]= variance.clusters[u,j]/ cluster.size[u]   
  }
  
  
  variance.matrix = numeric(0)
  for(j in 1:ncol(x)) 
    variance.matrix[j]=var(x[,j])*(n-1)/n
  
  
  Somme.variance.clusters = 0
  for (u in 1:k) 
    Somme.variance.clusters = Somme.variance.clusters+sqrt((variance.clusters[u,]%*%(variance.clusters[u,])))
  
  
  # Standard deviation
  stdev = (1/k)*sqrt(Somme.variance.clusters)
  
  #Average scattering for clusters  
  scat = (1/k)* (Somme.variance.clusters /sqrt(variance.matrix %*% variance.matrix))
  
  scat = list(stdev=stdev, centers=centers.matrix, variance.intraclusters= variance.clusters, scatt=scat)
  return(scat)
}

density.clusters = function(cl, x)
{
  x = as.matrix(x)
  k = max(cl)
  n = length(cl)
  
  distance = matrix(0, ncol = 1, nrow = n)
  density =  matrix(0, ncol = 1, nrow = k)
  centers.matrix = centers(cl,x)
  stdev = Average.scattering(cl,x)$stdev 
  for(i in 1:n) 
  {        
    u=1
    while(cl[i] != u )
      u = u+1
    for (j in 1:ncol(x))   
    {               
      distance[i] = distance[i]+(x[i,j]-centers.matrix[u,j])^2 
    }     
    distance[i] = sqrt(distance[i])            
    if (distance[i] <= stdev)
      density[u]= density[u]+1                      
  }  
  dens = list(distance=distance, density=density)    
  return(dens)          
  
}


density.bw = function(cl, x)
{
  x = as.matrix(x)
  k = max(cl)
  n = length(cl)   
  centers.matrix = centers(cl,x)
  stdev = Average.scattering(cl,x)$stdev 
  density.bw = matrix(0, ncol = k, nrow = k)
  u = 1
  
  for(u in 1:k)
  {
    for(v in 1:k)
    {
      if(v!=u)
      {  
        distance = matrix(0, ncol = 1, nrow = n)
        moy = (centers.matrix[u,]+centers.matrix[v,])/2
        for(i in 1:n)
        {
          if((cl[i]==u)||(cl[i]==v))
          {
            for (j in 1:ncol(x))   
            {               
              distance[i] = distance[i]+(x[i,j]-moy[j])^2 
            }   
            distance[i] = sqrt(distance[i])
            if(distance[i]<= stdev)
            {
              density.bw[u,v] = density.bw[u,v]+1                  
            }  
          }           
        }
      }       
    }
  }
  density.clust = density.clusters(cl,x)$density 
  S = 0
  for(u in 1:k)
    for(v in 1:k)
    {  
      if(max(density.clust[u], density.clust[v])!=0)
        S=S+ (density.bw[u,v]/max(density.clust[u], density.clust[v]))
    }   
  density.bw = S/(k*(k-1))
  return(density.bw) 
  
}      



Dis = function (cl, x)
{   # Dis : Total separation between clusters
  
  x = as.matrix(x)
  k = max(cl)
  centers.matrix = centers(cl,x)
  Distance.centers = dist(centers.matrix)
  Dmin = min(Distance.centers)
  Dmax = max(Distance.centers)
  Distance.centers = as.matrix(Distance.centers)
  s2 = 0
  for (u in 1:k)
  {
    s1=0
    for(j in 1:ncol(Distance.centers))
    {s1 = s1 + Distance.centers[u,j]       
    }
    s2 = s2+1/s1      
  }  
  Dis = (Dmax/Dmin)*s2  
  return(Dis)
}  



##################################
#                                #  
#         Hubert index           #
#                                #
##################################



Index.Hubert = function(x, cl)
{
  
  k = max(cl)
  n = dim(x)[1]
  y = matrix(0, ncol = dim(x)[2], nrow = n)
  P = as.matrix(md)
  meanP = mean(P)
  variance.matrix = numeric(0)
  md = dist(x, method="euclidean")
  for(j in 1:n) 
    variance.matrix[j]=var(P[,j])*(n-1)/n
  varP = sqrt(variance.matrix %*% variance.matrix)
  
  centers.clusters = centers(cl,x)
  for(i in 1:n)
  {
    for(u in 1:k)
    {
      if(cl[i]==u)
        y[i,] = centers.clusters[u,]
    }   
  }  
  
  Q = as.matrix(dist(y, method="euclidean"))
  meanQ = mean(Q)
  for(j in 1:n) 
    variance.matrix[j]=var(Q[,j])*(n-1)/n
  varQ = sqrt(variance.matrix %*% variance.matrix)
  
  M = n*(n-1)/2
  S = 0
  n1 = n-1
  
  for(i in 1:n1)
  { 
    j = i+1
    while(j<=n)
    {
      S = S+(P[i,j]-meanP)*(Q[i,j]-meanQ)
      j = j+1
    }
    
  } 
  gamma = S/(M*varP*varQ)
  
  return(gamma)
}   



##################################
#                                #  
#   Gamma, Gplus and Tau         #
#                                #
##################################    


Index.sPlussMoins = function (cl1,md)
{
  cn1 = max(cl1)
  n1 = length(cl1)
  dmat = as.matrix(md)
  average.distance = median.distance = separation = cluster.size = within.dist1 = between.dist1 = numeric(0)
  separation.matrix = matrix(0, ncol = cn1, nrow = cn1)
  di = list()
  for (u in 1:cn1) {
    cluster.size[u] = sum(cl1 == u)
    du = as.dist(dmat[cl1 == u, cl1 == u])
    within.dist1 = c(within.dist1, du)
    average.distance[u] = mean(du)
    median.distance[u] = median(du)
    bv = numeric(0)
    for (v in 1:cn1) {
      if (v != u) {
        suv = dmat[cl1 == u, cl1 == v]
        bv = c(bv, suv)
        if (u < v) {
          separation.matrix[u, v] = separation.matrix[v,u] = min(suv)
          between.dist1 = c(between.dist1, suv)
        }
      }
    }
  }
  
  nwithin1 = length(within.dist1)
  nbetween1 = length(between.dist1)
  meanwithin1 = mean(within.dist1)
  meanbetween1 = mean(between.dist1)
  
  s.plus = s.moins = 0 
  #s.moins = sum(rank(c(within.dist1,between.dist1),ties="first")[1:nwithin1]-rank(within.dist1,ties="first"))
  #s.plus  = sum(rank(c(-within.dist1,-between.dist1),ties="first")[1:nwithin1]-rank(-within.dist1,ties="first"))
  for (k in 1: nwithin1)
  {
    s.plus = s.plus+(colSums(outer(between.dist1,within.dist1[k], ">")))
    s.moins = s.moins+(colSums(outer(between.dist1,within.dist1[k], "<")))
  }    
  
  Index.Gamma = (s.plus-s.moins)/(s.plus+s.moins)
  Index.Gplus = (2*s.moins)/(n1*(n1-1))
  t.tau  = (nwithin1*nbetween1)-(s.plus+s.moins)
  Index.Tau = (s.plus-s.moins)/(((n1*(n1-1)/2-t.tau)*(n1*(n1-1)/2))^(1/2))
  
  results = list(gamma=Index.Gamma, gplus=Index.Gplus, tau=Index.Tau)
  return(results)
}




##################################
#                                #  
#      Frey and McClain          #
#                                #
################################## 




Index.15and28  = function (cl1,cl2,md)
{
  cn1 = max(cl1)
  n1 = length(cl1)
  dmat = as.matrix(md)
  average.distance = median.distance = separation = cluster.size = within.dist1 = between.dist1 = numeric(0)
  separation.matrix = matrix(0, ncol = cn1, nrow = cn1)
  di = list()
  for (u in 1:cn1) 
  {
    cluster.size[u] = sum(cl1 == u)
    du = as.dist(dmat[cl1 == u, cl1 == u])
    within.dist1 = c(within.dist1, du)
    #average.distance[u] = mean(du)
    #median.distance[u] = median(du)
    #bv = numeric(0)
    for (v in 1:cn1) {
      if (v != u) {
        suv = dmat[cl1 == u, cl1 == v]
        #bv = c(bv, suv)
        if (u < v) {
          separation.matrix[u, v] = separation.matrix[v,u] = min(suv)
          between.dist1 = c(between.dist1, suv)
        }
      }
    }
  }
  cn2 = max(cl2)
  n2 = length(cl2)
  dmat = as.matrix(md)
  average.distance = median.distance = separation = cluster.size = within.dist2 = between.dist2 = numeric(0)
  separation.matrix = matrix(0, ncol = cn2, nrow = cn2)
  di = list()
  for (w in 1:cn2) {
    cluster.size[w] = sum(cl2 == w)
    dw = as.dist(dmat[cl2 == w, cl2 == w])
    within.dist2 = c(within.dist2, dw)
    #average.distance[w] = mean(dw)
    #median.distance[w] = median(dw)
    bx = numeric(0)
    for (x in 1:cn2) {
      if (x != w) {
        swx = dmat[cl2 == w, cl2 == x]
        bx = c(bx, swx)
        if (w < x) {
          separation.matrix[w, x] = separation.matrix[x,w] = min(swx)
          between.dist2 = c(between.dist2, swx)
        }
      }
    }
  }
  nwithin1 = length(within.dist1)
  nbetween1 = length(between.dist1)
  meanwithin1 = mean(within.dist1)
  meanbetween1 = mean(between.dist1)
  meanwithin2 = mean(within.dist2)
  meanbetween2 = mean(between.dist2)
  Index.15 = (meanbetween2-meanbetween1)/(meanwithin2-meanwithin1)
  Index.28 = (meanwithin1/nwithin1)/(meanbetween1/nbetween1)
  
  results = list(frey=Index.15,mcclain=Index.28)
  return(results)
}


##################################
#                                #  
#      Point-biserial            #
#                                #
##################################  



Indice.ptbiserial = function (x,md,cl1)
{
  nn = dim(x)[1]
  pp = dim(x)[2]
  
  md2 = as.matrix(md)
  m01 = array(NA, c(nn,nn))
  nbr = (nn*(nn-1))/2
  pb = array(0,c(nbr,2))
  
  m3 = 1
  for (m1 in 2:nn)
  {
    m12 = m1-1
    for (m2 in 1:m12)
    {
      if (cl1[m1]==cl1[m2]) m01[m1,m2] = 0
      if (cl1[m1]!=cl1[m2]) m01[m1,m2] = 1
      pb[m3,1] = m01[m1,m2]
      pb[m3,2] = md2[m1,m2]
      m3 = m3+1
    }
  }
  
  y = pb[,1]
  x = pb[,2] 
  
  biserial.cor = function (x, y, use = c("all.obs", "complete.obs"), level = 1) 
  {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric variable.\n")
    y = as.factor(y)
    if (length(levs = levels(y)) > 2) 
      stop("'y' must be a dichotomous variable.\n")
    if (length(x) != length(y)) 
      stop("'x' and 'y' do not have the same length")
    use = match.arg(use)
    if (use == "complete.obs") {
      cc.ind = complete.cases(x, y)
      x = x[cc.ind]
      y = y[cc.ind]
    }
    ind = y == levs[level]
    diff.mu = mean(x[ind]) - mean(x[!ind])
    prob = mean(ind)
    diff.mu * sqrt(prob * (1 - prob))/sd(x)
  }
  
  ptbiserial = biserial.cor(x=pb[,2], y=pb[,1], level = 2)
  return(ptbiserial)
}


##########################################
#                                        #
#       Duda, pseudot2 and beale         #
#                                        #
##########################################


Indices.WKWL = function (x,cl1=cl1,cl2=cl2)
{
  dim2 = dim(x)[2]
  wss = function(x) 
  {
    x = as.matrix(x)
    n = length(x)
    centers = matrix(nrow = 1, ncol = ncol(x))
    
    if (ncol(x) == 1) 
    {	centers[1, ] = mean(x) 	}
    if (is.null(dim(x))) 
    {
      bb = matrix(x,byrow=FALSE,nrow=1,ncol=ncol(x))
      centers[1, ] = apply(bb, 2, mean)
    }
    else 
    {
      centers[1, ] = apply(x, 2, mean)
    }
    
    x.2 = sweep(x,2,centers[1,],"-")
    withins = sum(x.2^2)
    wss = sum(withins)
    return(wss)
  }
  
  ncg1 = 1
  ncg1max = max(cl1)
  while((sum(cl1==ncg1)==sum(cl2==ncg1)) && ncg1 <=ncg1max) 
  { ncg1 = ncg1+1 }
  g1 = ncg1
  
  ncg2 = max(cl2)
  nc2g2 = ncg2-1
  while((sum(cl1==nc2g2)==sum(cl2==ncg2)) && nc2g2 >=1) 
  { 
    ncg2 = ncg2-1 
    nc2g2 = nc2g2-1
  }
  g2 = ncg2
  
  NK = sum(cl2==g1)
  WK.x = x[cl2==g1,]
  WK = wss(x=WK.x)
  
  NL = sum(cl2==g2)
  WL.x = x[cl2==g2,]
  WL = wss(x=WL.x)
  
  NM = sum(cl1==g1)
  WM.x = x[cl1==g1,]
  WM = wss(x=WM.x)
  
  duda = (WK+WL)/WM
  
  BKL = WM-WK-WL
  pseudot2 = BKL/((WK+WL)/(NK+NL-2))
  
  beale = (BKL/(WK+WL))/(((NM-1)/(NM-2))*(2^(2/dim2)-1))
  
  results = list(duda=duda,pseudot2=pseudot2,NM=NM,NK=NK,NL=NL,beale=beale)
  return(results)
}


########################################################################
#                                                                      #
#       ccc, scott, marriot, trcovw, tracew, friedman and rubin        #
#                                                                      #
########################################################################    



Indices.WBT = function(x,cl,P,s,vv) 
{
  n = dim(x)[1]
  pp = dim(x)[2]
  qq = max(cl)
  z = matrix(0,ncol=qq,nrow=n)
  clX = as.matrix(cl)
  
  for (i in 1:n)
    for (j in 1:qq)
    {
      z[i,j]==0
      if (clX[i,1]==j) 
      {z[i,j]=1}
    }
  
  xbar = solve(t(z)%*%z)%*%t(z)%*%x
  B = t(xbar)%*%t(z)%*%z%*%xbar
  W = P-B
  marriot = (qq^2)*det(W)
  trcovw = sum(diag(cov(W)))
  tracew = sum(diag(W))
  if(det(W)!=0)
    scott = n*log(det(P)/det(W))
  else {cat("Error: division by zero!")}
  friedman = sum(diag(solve(W)*B))
  rubin = sum(diag(P))/sum(diag(W))
  
  
  R2 = 1-sum(diag(W))/sum(diag(P))
  v1 = 1
  u = rep(0,pp)
  c = (vv/(qq))^(1/pp)
  u = s/c
  k1 = sum((u>=1)==TRUE)
  p1 = min(k1,qq-1)
  if (all(p1>0,p1<pp))
  {
    for (i in 1:p1)
      v1 = v1*s[i]
    c = (v1/(qq))^(1/p1)
    u = s/c
    b1 = sum(1/(n+u[1:p1]))
    b2 = sum(u[p1+1:pp]^2/(n+u[p1+1:pp]),na.rm=TRUE)
    E_R2 = 1-((b1+b2)/sum(u^2))*((n-qq)^2/n)*(1+4/n)
    ccc = log((1-E_R2)/(1-R2))*(sqrt(n*p1/2)/((0.001+E_R2)^1.2))
  }else 
  {
    b1 = sum(1/(n+u))
    E_R2 = 1-(b1/sum(u^2))*((n-qq)^2/n)*(1+4/n)
    ccc = log((1-E_R2)/(1-R2))*(sqrt(n*pp/2)/((0.001+E_R2)^1.2))
  }
  results = list(ccc=ccc,scott=scott,marriot=marriot,trcovw=trcovw,tracew=tracew,friedman=friedman,rubin=rubin)
  return(results)
}



########################################################################
#                                                                      #
#                   Kl, Ch, Hartigan, Ratkowsky and Ball               #
#                                                                      #
########################################################################     




Indices.Traces = function(data, d, clall, index = "all") 
{
  x = data
  cl0 = clall[,1]
  cl1 = clall[,2]
  cl2 = clall[,3]
  clall = clall
  nb.cl0 = table(cl0)
  nb.cl1 = table(cl1)
  nb.cl2 = table(cl2)
  nb1.cl0 = sum(nb.cl0==1)
  nb1.cl1 = sum(nb.cl1==1)
  nb1.cl2 = sum(nb.cl2==1)
  
  gss = function(x, cl, d) 
  {
    n = length(cl)
    k = max(cl)
    centers = matrix(nrow = k, ncol = ncol(x))
    for (i in 1:k) 
    {
      
      if (ncol(x) == 1)
      {
        centers[i, ] = mean(x[cl == i, ])
      }
      if (is.null(dim(x[cl == i, ])))
      {
        bb = matrix(x[cl == i, ],byrow=FALSE,nrow=1,ncol=ncol(x))
        centers[i, ] = apply(bb, 2, mean)
      }
      else 
      {
        centers[i, ] = apply(x[cl == i, ], 2, mean)
      }
      
    }
    allmean = apply(x, 2, mean)
    dmean = sweep(x, 2, allmean, "-")
    allmeandist = sum(dmean^2)
    withins = rep(0, k)
    x.2 = (x - centers[cl, ])^2
    for (i in 1:k) {
      withins[i] = sum(x.2[cl == i, ])
    }
    wgss = sum(withins)
    bgss = allmeandist - wgss
    
    results = list(wgss=wgss,bgss=bgss, centers=centers)
    return(results)
  }
  
  vargss = function(x, clsize, varwithins) 
  {
    nvar = dim(x)[2]
    n = sum(clsize)
    k = length(clsize)
    varallmean = rep(0, nvar)
    varallmeandist = rep(0, nvar)
    varwgss = rep(0, nvar)
    for (l in 1:nvar) varallmean[l] = mean(x[, l])
    vardmean = sweep(x, 2, varallmean, "-")
    for (l in 1:nvar) {
      varallmeandist[l] = sum((vardmean[, l])^2)
      varwgss[l] = sum(varwithins[, l])
    }
    varbgss = varallmeandist - varwgss
    vartss = varbgss + varwgss
    zvargss = list(vartss = vartss, varbgss = varbgss)
    return(zvargss)
  }
  varwithinss = function(x, centers, cluster) {
    nrow = dim(centers)[1]
    nvar = dim(x)[2]
    varwithins = matrix(0, nrow, nvar)
    x = (x - centers[cluster, ])^2
    for (l in 1:nvar) {
      for (k in 1:nrow) {
        varwithins[k, l] = sum(x[cluster == k, l])
      }
    }
    return(varwithins)
  }
  
  
  
  indice.kl = function (x, clall, d = NULL, centrotypes = "centroids"){
    if (nb1.cl1 > 0){
      KL = NA
    }
    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
      stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d)) 
      stop("For argument centrotypes = 'medoids' d cannot be null")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d = as.matrix(d)
      }
      row.names(d) = row.names(x)
    }
    #if (is.null(dim(x))) {
    #	    dim(x) = c(length(x), 1)
    #}
    m = ncol(x)
    g = k = max(clall[, 2])
    KL = abs((g - 1)^(2/m) * gss(x, clall[, 1], d)$wgss - 
                g^(2/m) * gss(x, clall[, 2], d)$wgss)/abs((g)^(2/m) * 
                                                            gss(x, clall[, 2], d)$wgss - (g + 1)^(2/m) * 
                                                            gss(x, clall[, 3], d)$wgss)
    return(KL)
  }
  
  
  
  indice.ch = function (x, cl, d = NULL, centrotypes = "centroids"){
    if (nb1.cl1 > 0){
      CH = NA
    }
    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
      stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d)) 
      stop("For argument centrotypes = 'medoids' d cannot be null")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d = as.matrix(d)
      }
      row.names(d) = row.names(x)
    }
    #if (is.null(dim(x))) {
    #	    dim(x) = c(length(x), 1)
    #}
    n = length(cl)
    k = max(cl)
    CH = (gss(x, cl, d)$bgss/(k-1))/
      (gss(x, cl, d)$wgss/(n-k))
    return(CH)
  }
  
  
  # hartigan
  
  indice.hart = function(x, clall, d = NULL, centrotypes = "centroids"){
    if (sum(c("centroids", "medoids") == centrotypes) == 0) 
      stop("Wrong centrotypes argument")
    if ("medoids" == centrotypes && is.null(d)) 
      stop("For argument centrotypes = 'medoids' d cannot be null")
    if (!is.null(d)) {
      if (!is.matrix(d)) {
        d = as.matrix(d)
      }
      row.names(d) = row.names(x)
    }
    #if (is.null(dim(x))) {
    #    dim(x) = c(length(x), 1)
    #}
    n = nrow(x)
    g = max(clall[, 1])
    HART = (gss(x, clall[, 2], d)$wgss/gss(x, clall[, 3],d)$wgss - 1) * (n - g - 1)
    return(HART)
  }
  
  
  
  
  indice.ball = function(x, cl, d = NULL, centrotypes = "centroids"){
    wgssB = gss(x, cl, d)$wgss
    qq = max(cl)
    ball = wgssB/qq
    return(ball)
  }
  
  
  
  
  indice.ratkowsky = function(x, cl, d, centrotypes = "centroids"){
    qq = max(cl)
    clsize = table(cl)
    centers = gss(x, cl, d)$centers
    varwithins = varwithinss(x, centers, cl)
    zvargss = vargss(x, clsize, varwithins)
    ratio = mean(sqrt(zvargss$varbgss/zvargss$vartss))
    ratkowsky = ratio/sqrt(qq)
    return(ratkowsky)
  }
  
  indice = pmatch(index, c("kl", "ch", "hart", "ratkowsky", "ball", "all"))
  if (is.na(indice)) 
    stop("invalid clustering index")
  if (indice == -1) 
    stop("ambiguous index")
  vecallindex = numeric(5)
  if (any(indice == 1) || (indice == 6)) 
    vecallindex[1] = indice.kl(x,clall,d)
  if (any(indice == 2) || (indice == 6)) 
    vecallindex[2] = indice.ch(x,cl=clall[,2],d)
  if (any(indice == 3) || (indice == 6)) 
    vecallindex[3] = indice.hart(x,clall,d)
  if (any(indice == 4) || (indice == 6)) 
    vecallindex[4] = indice.ratkowsky(x,cl=cl1, d)
  if (any(indice == 5) || (indice == 6)) 
    vecallindex[5] = indice.ball(x,cl=cl1,d)
  names(vecallindex) = c("kl", "ch", "hart", "ratkowsky", "ball")
  if (indice < 6) 
    vecallindex = vecallindex[indice]
  return(vecallindex)
  
}



########################################################################
#                                                                      #
#                              C-index                                 #
#                                                                      #
######################################################################## 



Indice.cindex = function (d, cl) 
{
  d = data.matrix(d)
  DU = 0
  r = 0
  v_max = array(1, max(cl))
  v_min = array(1, max(cl))
  for (i in 1:max(cl)) {
    n = sum(cl == i)
    if (n > 1) {
      t = d[cl == i, cl == i]
      DU = DU + sum(t)/2
      v_max[i] = max(t)
      if (sum(t == 0) == n) 
        v_min[i] = min(t[t != 0])
      else v_min[i] = 0
      r = r + n * (n - 1)/2
    }
  }
  Dmin = min(v_min)
  Dmax = max(v_max)
  if (Dmin == Dmax) 
    result = NA
  else result = (DU - r * Dmin)/(Dmax * r - Dmin * r)
  result
}



########################################################################
#                                                                      #
#                                 DB                                   #
#                                                                      #
########################################################################     


Indice.DB = function (x, cl, d = NULL, centrotypes = "centroids", p = 2, q = 2) 
{
  if (sum(c("centroids") == centrotypes) == 0) 
    stop("Wrong centrotypes argument")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d = as.matrix(d)
    }
    row.names(d) = row.names(x)
  }
  if (is.null(dim(x))) {
    dim(x) = c(length(x), 1)
  }
  x = as.matrix(x)
  n = length(cl)
  k = max(cl)
  dAm = d
  centers = matrix(nrow = k, ncol = ncol(x))
  if (centrotypes == "centroids") {
    for (i in 1:k) {
      for (j in 1:ncol(x)) {
        centers[i, j] = mean(x[cl == i, j])
      }
    }
  }
  else {
    stop("wrong centrotypes argument")
  }
  S = rep(0, k)
  for (i in 1:k) {
    ind = (cl == i)
    if (sum(ind) > 1) {
      centerI = centers[i, ]
      centerI = rep(centerI, sum(ind))
      centerI = matrix(centerI, nrow = sum(ind), ncol = ncol(x), 
                        byrow = TRUE)
      S[i] = mean(sqrt(apply((x[ind, ] - centerI)^2, 1, 
                              sum))^q)^(1/q)
    }
    else S[i] = 0
  }
  M = as.matrix(dist(centers, p = p))
  R = array(Inf, c(k, k))
  r = rep(0, k)
  for (i in 1:k) {
    for (j in 1:k) {
      R[i, j] = (S[i] + S[j])/M[i, j]
    }
    r[i] = max(R[i, ][is.finite(R[i, ])])
  }
  DB = mean(r[is.finite(r)])
  resul = list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
  resul
}




########################################################################
#                                                                      #
#                             Silhouette                               #
#                                                                      #
########################################################################     



Indice.S = function (d, cl) 
{
  d = as.matrix(d)
  Si = 0
  for (k in 1:max(cl)) {
    if ((sum(cl == k)) <= 1) 
      Sil = 1
    else {
      Sil = 0
      for (i in 1:length(cl)) {
        if (cl[i] == k) {
          ai = sum(d[i, cl == k])/(sum(cl == k) - 1)
          dips = NULL
          for (j in 1:max(cl)) if (cl[i] != j) 
            if (sum(cl == j) != 1) 
              dips = cbind(dips, c((sum(d[i, cl == j]))/(sum(cl == 
                                                                j))))
          else dips = cbind(dips, c((sum(d[i, cl == 
                                              j]))))
          bi = min(dips)
          Sil = Sil + (bi - ai)/max(c(ai, bi))
        }
      }
    }
    Si = Si + Sil
  }
  Si/length(cl)
}



########################################################################
#                                                                      #
#                                  Gap                                 #
#                                                                      #
######################################################################## 

Indice.Gap = function (x, clall, reference.distribution = "unif", B = 10, 
                        method = "ward.D2", d = NULL, centrotypes = "centroids") 
{
  GAP = function(X, cl, referenceDistribution, B, method, d, centrotypes) 
  {
    set.seed(1)
    simgap = function(Xvec) 
    {
      ma = max(Xvec)
      mi = min(Xvec)
      set.seed(1)
      Xout = runif(length(Xvec), min = mi, max = ma)
      return(Xout)
    }
    pcsim = function(X, d, centrotypes) 
    {
      if (centrotypes == "centroids") 
      {
        Xmm = apply(X, 2, mean)
      }
      
      for (k in (1:dim(X)[2])) 
      {
        X[, k] = X[, k] - Xmm[k]
      }
      ss = svd(X)
      Xs = X %*% ss$v
      Xnew = apply(Xs, 2, simgap)
      Xt = Xnew %*% t(ss$v)
      for (k in (1:dim(X)[2])) {
        Xt[, k] = Xt[, k] + Xmm[k]
      }
      return(Xt)
    }
    if (is.null(dim(x))) 
    {
      dim(x) = c(length(x), 1)
    }
    ClassNr = max(cl)
    Wk0 = 0
    WkB = matrix(0, 1, B)
    for (bb in (1:B)) {
      if (reference.distribution == "unif") 
        Xnew = apply(X, 2, simgap)
      else if (reference.distribution == "pc") 
        Xnew = pcsim(X, d, centrotypes)
      else stop("Wrong reference distribution type")
      if (bb == 1) {
        pp = cl
        if (ClassNr == length(cl)) 
          pp2 = 1:ClassNr
        else if (method == "k-means") 
        { set.seed(1)
          pp2 = kmeans(Xnew, ClassNr, 100)$cluster
        }
        else if (method == "single" || method == "complete" || 
                 method == "average" || method == "ward.D2" || 
                 method == "mcquitty" || method == "median" || 
                 method == "centroid"|| method=="ward.D") 
          pp2 = cutree(hclust(dist(Xnew), method = method), 
                        ClassNr)
        else stop("Wrong clustering method")
        if (ClassNr > 1) {
          for (zz in (1:ClassNr)) {
            Xuse = X[pp == zz, ]
            Wk0 = Wk0 + sum(diag(var(Xuse))) * (length(pp[pp == 
                                                             zz]) - 1)/(dim(X)[1] - ClassNr)
            Xuse2 = Xnew[pp2 == zz, ]
            WkB[1, bb] = WkB[1, bb] + sum(diag(var(Xuse2))) * 
              (length(pp2[pp2 == zz]) - 1)/(dim(X)[1] - 
                                              ClassNr)
          }
        }
        if (ClassNr == 1) 
        {
          Wk0 = sum(diag(var(X)))
          WkB[1, bb] = sum(diag(var(Xnew)))
        }
      }
      if (bb > 1) {
        if (ClassNr == length(cl)) 
          pp2 = 1:ClassNr
        else if (method == "k-means")
        {
          set.seed(1)
          pp2 = kmeans(Xnew, ClassNr, 100)$cluster
        }
        else if (method == "single" || method == "complete" || 
                 method == "average" || method == "ward.D2" || 
                 method == "mcquitty" || method == "median" || 
                 method == "centroid"||method == "ward.D") 
          pp2 = cutree(hclust(dist(Xnew), method = method), 
                        ClassNr)
        else stop("Wrong clustering method")
        if (ClassNr > 1) {
          for (zz in (1:ClassNr)) {
            Xuse2 = Xnew[pp2 == zz, ]
            WkB[1, bb] = WkB[1, bb] + sum(diag(var(Xuse2))) * 
              length(pp2[pp2 == zz])/(dim(X)[1] - ClassNr)
          }
        }
        if (ClassNr == 1) {
          WkB[1, bb] = sum(diag(var(Xnew)))
        }
      }
    }
    Sgap = mean(log(WkB[1, ])) - log(Wk0)
    Sdgap = sqrt(1 + 1/B) * sqrt(var(log(WkB[1, ]))) * sqrt((B - 
                                                                1)/B)
    resul = list(Sgap = Sgap, Sdgap = Sdgap)
    resul
  }
  if (sum(c("centroids", "medoids") == centrotypes) == 0) 
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d)) 
    stop("For argument centrotypes = 'medoids' d can not be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d = as.matrix(d)
    }
    row.names(d) = row.names(x)
  }
  X = as.matrix(x)
  gap1 = GAP(X, clall[, 1], reference.distribution, B, method, 
              d, centrotypes)
  gap = gap1$Sgap
  gap2 = GAP(X, clall[, 2], reference.distribution, B, method, 
              d, centrotypes)
  diffu = gap - (gap2$Sgap - gap2$Sdgap)
  resul = list(gap = gap, diffu = diffu)
  resul
  
}




########################################################################
#                                                                      #
#                              SD, sdbw, dunn                          #
#                                                                      #
########################################################################   




Index.sdindex = function(x, clmax, cl)
{  
  x = as.matrix(x)
  Alpha = Dis(clmax,x)
  Scatt = Average.scattering(cl,x)$scatt
  Dis0 = Dis(cl,x)
  SD.indice = Alpha*Scatt + Dis0
  return(SD.indice)
}

Index.SDbw = function(x, cl)
{
  x = as.matrix(x)
  Scatt = Average.scattering(cl,x)$scatt
  Dens.bw = density.bw(cl,x)
  SDbw = Scatt+Dens.bw
  return(SDbw)
}    



########################################################################
#                                                                      #
#                              D index                                 #
#                                                                      #
######################################################################## 




Index.Dindex = function(cl, x)
{
  x = as.matrix(x)
  distance = density.clusters(cl, x)$distance
  n = length(distance)
  S = 0
  for(i in 1:n)
    S = S+distance[i]
  inertieIntra = S/n
  return(inertieIntra)
}    


#####################################################################
#                                                                   #
#                            Dunn index                             #
#                                                                   #
#####################################################################



Index.dunn = function(md, clusters, Data=NULL, method="euclidean")
{
  
  distance = as.matrix(md)
  nc = max(clusters)
  interClust = matrix(NA, nc, nc)
  intraClust = rep(NA, nc)
  
  for (i in 1:nc) 
  {
    c1 = which(clusters==i)
    for (j in i:nc) {
      if (j==i) intraClust[i] = max(distance[c1,c1])
      if (j>i) {
        c2 = which(clusters==j)
        interClust[i,j] = min(distance[c1,c2])
      }
    }
  }
  dunn = min(interClust,na.rm=TRUE)/max(intraClust)
  return(dunn)
}



################
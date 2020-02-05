#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                                                                      #
#                                              Distances                                                               #
#                                                                                                                      #
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

init_distances  =  function(distance, diss, jeu) {
  if(is.null(distance))
    distanceM = 7
  if(!is.null(distance))
    distanceM  =  pmatch(distance, c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  
  if (is.na(distanceM)) 
  {
    stop("invalid distance")
  } 
  
  if(is.null(diss))
  {  
    
    if (distanceM == 1) 
    {
      md  =  dist(jeu, method="euclidean")	# "dist" function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
    }
    if (distanceM == 2) 
    {
      md  =  dist(jeu, method="maximum")	
    }
    if (distanceM == 3) 
    {
      md  =  dist(jeu, method="manhattan")	
    }
    if (distanceM == 4) 
    {
      md  =  dist(jeu, method="canberra")	
    }
    if (distanceM == 5) 
    {
      md  =  dist(jeu, method="binary")	
    }
    if (distanceM == 6) 
    {
      md  =  dist(jeu, method="minkowski")	
    }
    
    if (distanceM == 7) 
    {		  
      stop("dissimilarity matrix and distance are both NULL")		
    } 
  }
  
  if(!is.null(diss))
  {
    if((distanceM==1)||(distanceM==2)|| (distanceM==3)|| (distanceM==4)|| (distanceM==5)|| (distanceM==6))
      stop("dissimilarity matrix and distance are both not null")
    else
      md  =  diss
  }
  
  return(md)
}

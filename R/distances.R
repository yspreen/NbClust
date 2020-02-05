#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                                                                      #
#                                              Distances                                                               #
#                                                                                                                      #
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

init_distances = function(distance, diss, game) {
  if(is.null(distance))
    distanceM = 7
  if(!is.null(distance))
    distanceM = pmatch(distance, c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  
  if (is.na(distanceM)) 
  {
    stop("invalid distance")
  } 
  
  if(is.null(diss))
  {  
    
    if (distanceM == 1) 
    {
      md = dist(game, method="euclidean")	# "dist" function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
    }
    if (distanceM == 2) 
    {
      md = dist(game, method="maximum")	
    }
    if (distanceM == 3) 
    {
      md = dist(game, method="manhattan")	
    }
    if (distanceM == 4) 
    {
      md = dist(game, method="canberra")	
    }
    if (distanceM == 5) 
    {
      md = dist(game, method="binary")	
    }
    if (distanceM == 6) 
    {
      md = dist(game, method="minkowski")	
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
      md = diss
  }
  
  return(md)
}

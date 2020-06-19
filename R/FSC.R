
# This file includes all the FSC methods for clustering functional data
# modified 29 May 2020
# --------------------

require(fda)
require(fda.usc)

#===================
# Spectral Clustering function FSC

# This function performs functional spectral clustering 
# and returns the cluster memebers of the curves 
# The function reads the actual data, timeline, basis, number of clusters,
# and the value d, where d can be 0, or 1 or 2.

FSC.fun <- function(matrix.data , timeline , basis ,  nclusters , d){ 
  
  if(class(matrix.data)=="matrix"){
    
    # find the number of curves
    n <- dim(matrix.data)[2] 
    # smooth the data according to the supplied basis
    smoothedfd <- smooth.basis( timeline , matrix.data, basis)
    deriv1 <- deriv.fd(smoothedfd$fd, 1) # for plotting purposes
    deriv2 <- deriv.fd(smoothedfd$fd, 2) # for plotting purposes
    
    # evaluate the smoothed data
    evaluatedata <- eval.fd(timeline, as.fd(smoothedfd) )
    # calculate the distance using fda.usc metric.lp() measure
    fobject <- fdata(t(evaluatedata))  
    # assign the type of clustering according to d
    if(d==0){  
      M <- metric.lp(fobject) 
    }else if(d==1){  
      M <- semimetric.deriv(fobject, nderiv = 1)
    }else if(d==2){  
      M <- semimetric.deriv(fobject, nderiv = 2)
    }else{
      print("Warning: d must take the value: 0 or 1 or 2.")
    } 
    # calculate sd of the elements of M
    sdM <- sd(M)
    # form the afinity matrix
    A <- exp(-M/(2*sdM))
    # form the Diagonal Matrix D, where D = sum of W across rows
    D <- diag(apply(A,1,sum), n, n)  
    # form the L matrix, such that L=D*AD*, where *=-1/2
    L <- diag(1/sqrt(diag(D))) %*% A %*% diag(1/sqrt(diag(D)))
    # form the matrix X from the eigenvectors
    X <- eigen(L, symmetric = TRUE)$vectors[,1:nclusters] 
    # form Y from X
    Y <- X/sqrt(rowSums(X^2))
    # Cluster Y rows using (k-means)
    kmY <- kmeans(Y, centers = nclusters, iter.max = 100L, nstart = 1000L)
    # get the clusters members
    cls.kmY <- kmY$cluster
    
    # to plot the curves in their clusters
    par(mfrow=c(1,3))
    plot(smoothedfd, col=cls.kmY ,lwd=1,lty=1, main="Smoothed curves")
    plot(deriv1, col=cls.kmY ,lwd=1,lty=1, main= "first derivatives")
    plot(deriv2, col=cls.kmY ,lwd=1,lty=1, main="second derivatives")
    
    return(cls.kmY)
    
  }else{
    print("the data must be in a matrix format")
  }
  
  
} 


#===================



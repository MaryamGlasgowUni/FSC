# This file includes all the FSC methods for clustering functional data
# modified 29 May 2020
# --------------------

require(fda)
require(fda.usc)

#===================
# Spectral Clustering function FSC-S (Do)

# This function does the spectral clustering using original curves.
# The function reads the actual data, basis, and the number of clusters
# The function returns a matrix with the clusters members. 

fsc.clst0.fun <- function(matrix.data, basis,  nclusters){ 
  
  n <- dim(matrix.data[[1]])[2] 
  v <- length(matrix.data) 
  
  smoothedfd <- list()
  evaluatedata <- list()
  fobject <- list()
  M <- array(rep(0,n^2), dim=c(n,n,v))
  sdM <- numeric(v)
  A <- array(rep(0,n^2), dim=c(n,n,v)) 
  D <- list()
  L <- list()
  X <- list()
  Y <- list()
  kmY <- list()
  cls.kmY <- matrix(0,n,v) 
  
  
  for(i in 1:v){
    
    smoothedfd[[i]] <- smooth.basis( t , matrix.data[[i]], basis)
    # evaluate the smoothed data
    evaluatedata[[i]] <- eval.fd(t, as.fd(smoothedfd[[i]]) )
    # calculate the distance using fda.usc metric.lp() measure
    fobject[[i]] <- fdata(t(evaluatedata[[i]]))  
    
    M[,,i] <- metric.lp(fobject[[i]]) 
    
    sdM[i] <- sd(M[,,i])
    
    # form the afinity matrix
    A[,,i] <- exp(-M[,,i]/(2*sdM[i]))
    
    # form the Diagonal Matrix D, where D = sum of W across rows
    D[[i]] = diag(apply(A[,,i],1,sum), n, n)  
    
    # form the L matrix, such that L=D*AD*, where *=-1/2
    L[[i]] = diag(1/sqrt(diag(D[[i]]))) %*% A[,,i] %*% diag(1/sqrt(diag(D[[i]])))
    
    # form the matrix X from the eigenvectors
    X[[i]] = eigen(L[[i]], symmetric = TRUE)$vectors[,1:nclusters] 
    
    # form Y from X
    Y[[i]] = X[[i]]/sqrt(rowSums(X[[i]]^2))
    
    # Cluster Y rows using (k-means)
    kmY[[i]] = kmeans(Y[[i]], centers = nclusters, iter.max = 100L, nstart = 1000L)
    # get the clusters members
    cls.kmY[,i] <- kmY[[i]]$cluster
    
    
    # to plot the curves in their clusters
    #par(mfrow=c(1,2))
    #matplot(matrix.data[[i]] ,  type = "l", col = cls.kmY[,i],main="Original curves")
    #plot(smoothedfd[[i]], col=cls.kmY[,i] ,lwd=1,lty=1, main="Smoothed curves")
    cat("Cluster data set", i , "of", v, "data sets is done!\n")
    
  }
  
  return(cls.kmY)
  
} 


#===================
# Spectral Clustering function FSC-S (D1)

# This function does the spectral clustering using first derivatives.
# The function reads the actual data, basis, and the number of clusters
# The function returns a matrix with the clusters members. 

fsc.clst1.fun <- function(matrix.data, basis, nclusters){ 
  
  n <- dim(matrix.data[[1]])[2] 
  v <- length(matrix.data) 
  
  smoothedfd <- list()
  evaluatedata <- list()
  fobject <- list()
  M <- array(rep(0,n^2), dim=c(n,n,v))
  sdM <- numeric(v)
  A <- array(rep(0,n^2), dim=c(n,n,v)) 
  D <- list()
  L <- list()
  X <- list()
  Y <- list()
  kmY <- list()
  cls.kmY <- matrix(0,n,v) 
  
  
  
  for(i in 1:v){
    
    smoothedfd[[i]] <- smooth.basis( t , matrix.data[[i]], penfd)
    # evaluate the smoothed data
    evaluatedata[[i]] <- eval.fd(t, as.fd(smoothedfd[[i]]) )
    # calculate the distance using fda.usc metric.lp() measure
    fobject[[i]] <- fdata(t(evaluatedata[[i]]))  
    
    M[,,i] <- semimetric.deriv(fobject[[i]], nderiv = 1)
    sdM[i] <- sd(M[,,i])
    
    # form the afinity matrix
    A[,,i] <- exp(-M[,,i]/(2*sdM[i]))
    
    # form the Diagonal Matrix D, where D = sum of W across rows
    D[[i]] = diag(apply(A[,,i],1,sum), n, n)  
    
    # form the L matrix, such that L=D*AD*, where *=-1/2
    L[[i]] = diag(1/sqrt(diag(D[[i]]))) %*% A[,,i] %*% diag(1/sqrt(diag(D[[i]])))
    
    # form the matrix X from the eigenvectors
    X[[i]] = eigen(L[[i]], symmetric = TRUE)$vectors[,1:3] 
    
    # form Y from X
    Y[[i]] = X[[i]]/sqrt(rowSums(X[[i]]^2))
    
    # Cluster Y rows using (k-means)
    kmY[[i]] = kmeans(Y[[i]], centers = nclusters, iter.max = 100L, nstart = 1000L)
    # get the clusters members
    cls.kmY[,i] <- kmY[[i]]$cluster
    
    # to plot the curves in their clusters
    #par(mfrow=c(1,2))
    #matplot(matrix.data[[i]] ,  type = "l", col = cls.kmY[,i],main="Original curves")
    #plot(smoothedfd[[i]], col=cls.kmY[,i] ,lwd=1,lty=1, main="Smoothed curves")
    cat("Cluster data set", i , "of", v, "data sets is done!\n")
    
  }
  
  return(cls.kmY)
  
} 

#==================
# Spectral Clustering function FSC-S (D2)

# This function does the spectral clustering using second derivatives.
# The function reads the actual data, basis, and the number of clusters
# The function returns a matrix with the clusters members. 

fsc.clst2.fun <- function(matrix.data, basis, nclusters){ 
  
  n <- dim(matrix.data[[1]])[2] 
  v <- length(matrix.data) 
  
  smoothedfd <- list()
  evaluatedata <- list()
  fobject <- list()
  M0 <- array(rep(0,n^2), dim=c(n,n,v))
  M1 <- array(rep(0,n^2), dim=c(n,n,v))
  M <- array(rep(0,n^2), dim=c(n,n,v))
  sdM <- numeric(v)
  A <- array(rep(0,n^2), dim=c(n,n,v)) 
  D <- list()
  L <- list()
  X <- list()
  Y <- list()
  kmY <- list()
  cls.kmY <- matrix(0,n,v) 
  
  
  
  for(i in 1:v){
    
    smoothedfd[[i]] <- smooth.basis( t , matrix.data[[i]], penfd)
    # evaluate the smoothed data
    evaluatedata[[i]] <- eval.fd(t, as.fd(smoothedfd[[i]]) )
    # calculate the distance using fda.usc metric.lp() measure
    fobject[[i]] <- fdata(t(evaluatedata[[i]]))  
    
    M[,,i] <- semimetric.deriv(fobject[[i]], nderiv = 2)
    sdM[i] <- sd(M[,,i])
    
    # form the afinity matrix
    A[,,i] <- exp(-M[,,i]/(2*sdM[i]))
    
    # form the Diagonal Matrix D, where D = sum of W across rows
    D[[i]] = diag(apply(A[,,i],1,sum), n, n)  
    
    # form the L matrix, such that L=D*AD*, where *=-1/2
    L[[i]] = diag(1/sqrt(diag(D[[i]]))) %*% A[,,i] %*% diag(1/sqrt(diag(D[[i]])))
    
    # form the matrix X from the eigenvectors
    X[[i]] = eigen(L[[i]], symmetric = TRUE)$vectors[,1:3] 
    
    # form Y from X
    Y[[i]] = X[[i]]/sqrt(rowSums(X[[i]]^2))
    
    # Cluster Y rows using (k-means)
    kmY[[i]] = kmeans(Y[[i]], centers = nclusters, iter.max = 100L, nstart = 1000L)
    
    # get the clusters members
    cls.kmY[,i] <- kmY[[i]]$cluster
    
    # to plot the curves in their clusters
    #par(mfrow=c(1,2))
    #matplot(matrix.data[[i]] ,  type = "l", col = cls.kmY[,i],main="Original curves")
    #plot(smoothedfd[[i]], col=cls.kmY[,i] ,lwd=1,lty=1, main="Smoothed curves")
    cat("Cluster data set", i , "of", v, "data sets is done!\n")
    
  }
  
  return(cls.kmY)
  
} 


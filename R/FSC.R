FSC.fun <- function(matrix.data , timeline , basis ,  nclusters , d, ...){

   if(class(matrix.data)!="matrix")
     stop("the data must be in a matrix format")

   if(!missing(timeline)){
     timeline <- timeline
   }else{
     timeline <- 1:nrow(matrix.data)
   }

  ifelse(!missing(nclusters) , nclusters <- nclusters, print("Please choose number of clusters"))

      # find the number of curves
      n <- dim(matrix.data)[2]
      # smooth the data according to the supplied basis
      smoothedfd <- smooth.basis( timeline , matrix.data, basis)
      deriv1 <- deriv.fd(smoothedfd$fd, 1)
      deriv2 <- deriv.fd(smoothedfd$fd, 2)

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
      cls <- kmY$cluster
      # more information about the clusters
      cls.size <- kmY$size

         # to plot the curves in their clusters
         par(mfrow=c(1,3))
         plot(smoothedfd, col=cls ,lwd=1,lty=1, main="Smoothed curves")
         plot(deriv1, col=cls ,lwd=1,lty=1, main= "first derivatives")
         plot(deriv2, col=cls ,lwd=1,lty=1, main="second derivatives")

    returnlist <- list(cls , cls.size, smoothedfd$fd$coefs, deriv1$coefs, deriv2$coefs, M)
    names(returnlist) <- c("clusters", "clusters.size","coefs.fd", "coefs.deriv1", "coefs.deriv2", "dist.matrix")

    return(returnlist)


}



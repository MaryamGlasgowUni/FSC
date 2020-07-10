#' Functional Spectral Clustering
#'
#' This function introduces an approach for clustering functional data based on spectral clustering analysis.
#' @param matrix.data is the data to be clustered in a matrix format.
#' @param timeline is the time interval of the data.
#' @param basis is some generated smooth bases appropriate for the data, see example below.
#' @param nclusters the number of clusters in the data.
#' @param d d can take the value 0 or 1 or 2, depends on the chosen set of curves for clustering the data.
#' @param ... extra arguments can be passed to the function.
#' @return FSC.fun returns a list of objects
#' \item{clusters}{the curves labels according to their clusters.}
#' \item{clusters.size}{the number of curves in each cluster.}
#' \item{coefs.fd}{the coefficients of the smoothed curves.}
#' \item{coefs.deriv1}{the coefficients of the first derivatives.}
#' \item{coefs.deriv2}{the coefficients of the second derivatives.}
#' \item{dist.matrix}{the distance matrix of the curves based on the choice of d.}
#' @details The clustering procedure is straightforward once the smoothing technique is chosen with care. The smoothing stage plays an important role in this technique and can determine the success of the clustering results, see \href{https://cran.r-project.org/web/packages/fda/fda.pdf}{fda pacakge} for creating basis functions. While the number of clusters in the data can be determined by different ways, for instance \href{https://rdrr.io/cran/NbClust/man/NbClust.html}{NbClust} can be used to estimate the number of clusters of the original data. Also, in some cases the number of clusters is known priori from the data. Finally, regarding the choice of d; usually clustering functional data is based on the original trajectories. However, it has been shown that the first derivatives (rate of change) can sometimes hold more information about the data and accordingly it can detect the similarities/dissimilarities better. Similary the same concept applied for the second derivatives (accelaration).
#' @examples # Apply the FSC technique of the Canadian weather temperature data
#' y <- CanadianWeather$dailyAv[,,1]
#' t <- 1:365
#' # plot the original observations over the timeline
#' matplot(y, type="l", main=" temperature observations for 35 Canadian cities", xlab="days", ylab="temp C")
#'
#' # choose smoothing technique for the data y,
#' bbasis.y <- create.bspline.basis(rangeval= c(0, 365), nbasis = 367, norder= 4)
#' penfd.y <- fdPar(bbasis.y,Lfdobj=int2Lfd(2),lambda= 10^{4})
#' smooth.y <- smooth.basis(1:365, y , penfd.y)
#' # plot the smoothed curves after applying the basis functions
#' plot(smooth.y ,main="smoothed curves of the daily tempreture curves", xlab="over a year", ylab="temp C")
#'
#' # find the first and the second derivatives of the smoothed curves:
#' deriv1.y <- deriv.fd(smooth.y$fd, 1)
#' deriv2.y <- deriv.fd(smooth.y$fd, 2)
#' # plot the first and the second derivatives to examine them
#' par(mfrow=c(1,2))
#' plot(deriv1.y, main="first derivatives")
#' plot(deriv2.y, main="second derivatives")
#'
#' # Apply the clustering technique. According to (Ramsay, and Silverman, 2005), we will assume the number of clusters = 4.
#'
#' # using original curves for clustering
#' clust.d0 <- FSC.fun(y, t, penfd.y, 4, 0)
#' # using first derivatives for clustering
#' clust.d1 <- FSC.fun(y, t, penfd.y, 4, 1)
#' # using second derivatives for clustering
#' clust.d2 <- FSC.fun(y, t, penfd.y, 4, 2)
#'
#'
#' @import fda fda.usc
#'
#' @export

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



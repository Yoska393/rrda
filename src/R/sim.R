#' @title Generate simulated data for Ridge Redundancy Analysis (RDA).
#' @description The function generates simulated data for Ridge Redundancy Analysis (RDA). It creates two data matrices, X and Y, based on a set of shared latent variables H. The function adds noise to the data and returns a list containing the matrices X, Y, the latent variables H, and the regression coefficients theta.y used for generating Y.
#' @param n The number of samples.
#' @param p The number of variables of X.
#' @param q The number of variables of Y.
#' @param k The number of latent variables.
#' @param s2n The numeric parameters of signal to noise ratio for X and Y, default value is c(1,1).
#' @return A list containing matrices X, Y, H, and theta.y.
#' @importFrom stats rnorm
#' @export
#' @examples
#' # Example usage of rdasim1
#' set.seed(10)
#' sim_data <- rdasim1(n = 10, p = 5, q = 3, k = 2)
#' str(sim_data)
rdasim1<- function(n,p,q,k,s2n = c(5, 5)){
  H  <- matrix(stats::rnorm(n = n*k, mean = 0, sd = 1),n,k)
  
  theta.x  <- matrix(stats::rnorm(n = k*p, mean = 0, sd = 1),k,p)
  noise.x <-  matrix(stats::rnorm(n = n*p, mean = 0, sd = 1),n,p)
  
  signal.x <- H %*% theta.x
  sigma.x <- sqrt(sum(signal.x^2) / (sum(noise.x^2) * s2n[1]))
  noise.x <- sigma.x * noise.x

  X  <- signal.x + noise.x
  
  theta.y  <- matrix(stats::rnorm(n = k*q, mean = 0, sd = 1),k,q)
  noise.y  <- matrix(stats::rnorm(n = n*q, mean = 0, sd = 1),n,q)
  
  signal.y <- H %*% theta.y
  sigma.y <- sqrt(sum(signal.y^2) / (sum(noise.y^2) * s2n[2]))
  noise.y <- sigma.y * noise.y


  Y  <- signal.y + noise.y

  return(list(X=X,Y=Y,H=H,theta.y=theta.y))
}

#' @title Generate simulated data for Ridge Redundancy Analysis (RDA).
#' @description The rdasim2 function generates simulated data for Ridge Redundancy Analysis (RDA) with adjustable signal-to-noise ratio and covariance structure for X. The data matrix Y is created by a low-rank model, where the rank is set by the product of two matrices A and C corresponding to the number of latent variables (k). The function allows control over the signal-to-noise ratio (s2n) and off-diagonal elements of the covariance matrix for X (xofd). It returns a list containing the matrices X, Y, the regression coefficient matrix B (obtained as the product of A and C), and the error matrix E.
#' @param n The number of samples.
#' @param p The number of variables of X.
#' @param q The number of variables of Y.
#' @param k The number of latent variables.
#' @param s2n The numeric parameter of signal to noise ratio, default value is 5.
#' @param xofd The numeric parameter of the off-diagnal elements of covariance matrix of X, default is 0.
#' @return A list containing matrices X, Y, B, E.
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @export
#' @examples
#' # Example usage of rdasim2
#' set.seed(10)
#' sim_data2 <- rdasim2(n = 10, p = 5, q = 3, k = 2)
#' str(sim_data2)
rdasim2<- function(n,p,q,k,s2n=5,xofd=0){
  cov <- matrix(xofd, nrow = p, ncol = p)
  diag(cov) <- 1
  X  <- MASS::mvrnorm(n, rep(0, p), cov)

  A <- matrix(stats::rnorm(n = p * k, mean = 0, sd = 1), p, k)
  A <- apply(A, 2, function(col) col / sqrt(sum(col^2)))
  C <- matrix(stats::rnorm(n = k * q, mean = 0, sd = 1), k, q)
  C <- apply(C, 2, function(col) col / sqrt(sum(col^2)))

  B  <- crossprod(t(A),C)
  E  <- matrix(stats::rnorm(n = n*q, mean = 0, sd = 1),n,q)
  XB <- crossprod(t(X),B)
  sigma <- sqrt(sum(XB^2) / (sum(E^2) * s2n))
  E  <- E * sigma
  Y <- XB + E
  return(list(X=X,Y=Y,B=B,E=E))
}



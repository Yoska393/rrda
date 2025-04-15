#' @title Compute the square root of the inverse of (d^2 + lambda).
#' @description Compute the square root of the inverse of (d^2 + lambda).
#' @param d A scalar of singular value of X.
#' @param lambda A scalar of the ridge penalty.
#' @return A scalar of the square root of the inverse of (d^2 + lambda).
sqrt_inv_d2_lambda <- function(d,lambda) {
  return(1 / sqrt(d**2 + lambda))
}


# rrda Ridge Redundancy Analysis
#' @title Calculate the coefficient Bhat by Ridge Redundancy Analysis.
#' @description This function performs Ridge Redundancy Analysis (RRDA) to obtain the coefficient Bhat, which models the relationship between a matrix of response variables (Y; \eqn{n \times q} matrix) and a matrix of explanatory variables (X; \eqn{n \times p} matrix). Especially, the function is designed to facilitate a high-dimensional computation and storing (for the details, refer to the article Yoshioka et al. 2025).
#'
#' The Ridge Redundancy Analysis model is represented as:
#' \deqn{Y = XB + E}
#' where:
#' - \eqn{Y} is the response matrix (\eqn{n \times q}),
#' - \eqn{X} is the predictor matrix (\eqn{n \times p}),
#' - \eqn{B} is the regression coefficient matrix (\eqn{p \times q}), and
#' - \eqn{E} is the error matrix (\eqn{n \times q}).
#'
#' The regularized estimate of \eqn{B} is described as:
#' \deqn{\hat{B}(\lambda) = \left(X'X + \lambda P_{X'}\right)^{-} X'Y}
#'
#' Additionally, the regularized-rank-restricted estimation of \eqn{B} is represented as:
#' \deqn{\hat{B}(\lambda, r) = U_{\hat{B}(\lambda)}^{[r]} D_{\hat{B}(\lambda)}^{[r]} V_{\hat{B}(\lambda)}^{[r]'}}
#' Here:
#' - \eqn{U_{\hat{B}(\lambda)}^{[r]}} is a \eqn{p \times r} matrix,
#' - \eqn{D_{\hat{B}(\lambda)}^{[r]}} is a \eqn{r \times r} diagonal matrix, and
#' - \eqn{V_{\hat{B}(\lambda)}^{[r]}} is a \eqn{q \times r} matrix.
#'
#' The user can specify ranks (nrank), ridge penalty values (lambda), and whether to center and scale the X and Y matrices.
#'
#' The Bhat can be returned as either component vectors or matrices. To store a large size of matrix, the coefficient Bhat is by default stored as LeftBhatlambda_k (F; \eqn{p \times r} matrix) and RightBhatlambda_k (G; \eqn{q \times r}). Here, r is the specified rank (nrank) in the Ridge Redundancy Analysis formula.
#'
#' For \eqn{i = 1, \ldots, r}, the matrices \eqn{F} and \eqn{G} are defined as:
#' \deqn{F_{.i} = U_{\hat{B}(\lambda)}^{[i]}D_{\hat{B}(\lambda)}^{[i]}, \quad G_{.i} = V_{\hat{B}(\lambda)}^{[i]}}
#'
#' These definitions allow the decomposition of \eqn{\hat{B}(\lambda)} into rank-specific components, facilitating the storing of the high-dimensional regression coefficients. To reconstruct the matrix form of Bhat, you can use the `rrda.coef` function.
#'
#'
#' @param Y A numeric matrix of response variables.
#' @param X A numeric matrix of explanatory variables.
#' @param nrank A numeric vector specifying the ranks of Bhat. Default is `NULL`, which sets it to `(1:min(15, min(dim(X), dim(Y))))`.
#' @param lambda A numeric vector of ridge penalty values. Default value is 1.
#' @param component Logical indicating if Bhat is returned as vectors or matrices. If `TRUE`, returns Bhat as component vectors. If `FALSE`, returns Bhat as matrices.
#' @param center.X Logical indicating if `X` should be centered. If `TRUE`, scales `X`. Default is `TRUE`.
#' @param center.Y Logical indicating if `Y` should be centered. If `TRUE`, scales `Y`. Default is `TRUE`.
#' @param scale.X Logical indicating if `X` should be scaled. If `TRUE`, scales `X`. Default is `FALSE`.
#' @param scale.Y Logical indicating if `Y` should be scaled. If `TRUE`, scales `Y`. Default is `FALSE`.
#' @return A list containing Bhat components or Bhat matrices (the coefficient of Ridge Redundancy Analysis for each parameter lambda and nrank), ranks, and lambda values.
#' @export
#' @importFrom stats sd
#' @importFrom furrr future_map
#' @importFrom RSpectra svds
#' @examples
#' set.seed(10)
#' simdata<-rdasim1(n = 100,p = 200,q = 200,k = 5)
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' # Sequential
#' Bhat <- rrda.fit(Y = Y, X = X, nrank = c(1:10))
#' names(Bhat)
rrda.fit<- function(Y, X, nrank = NULL, lambda= 1, component = TRUE, center.X = TRUE, center.Y = TRUE, scale.X = FALSE, scale.Y = FALSE){

  Y <- as.matrix(Y)
  X <- as.matrix(X)

  if (is.null(nrank)){
    nrank <- seq_len(min(15,min(dim(X),dim(Y))))
  }
  if (!is.numeric(nrank)){
    stop("nrank must be numeric \n")
  }else if (max(nrank) > min(dim(X),dim(Y)) ){
    stop("rank(B) must be less than or equal to rank(X) and rank(Y) \n")
  }

  if (any(lambda < 0)){
    stop("lambdas should be non-negative")
  }

  nrank<-sort(nrank)
  lambda<- sort(lambda)

  meansXtr <- NULL
  stdXtr <- NULL
  meansYtr <- NULL
  stdYtr <- NULL

  if(center.X){
    meansXtr <- colMeans(X)
    X <- scale(X, center = meansXtr, scale = FALSE)
    if (scale.X){
      stdXtr <- apply(X, 2, stats::sd)
      X <- scale(X, center = meansXtr, scale = stdXtr)
      X[is.nan(X)] <- 0
    }
  }else{
    if (scale.X){
      warning("X Scaling is not performed when center.X = FALSE")
    }
  }

  if(center.Y){
    meansYtr <- colMeans(Y)
    Y <- scale(Y, center = meansYtr, scale = FALSE)
    if (scale.Y){
      stdYtr <- apply(Y, 2, stats::sd)
      Y <- scale(Y, center = meansYtr, scale = stdYtr)
      Y[is.nan(Y)] <- 0
    }
  }else{
    if (scale.Y){
      warning("Y Scaling is not performed when center.Y = FALSE")
    }
  }

  scalelist<-list(meansXtr=meansXtr,
                  stdXtr=stdXtr,
                  meansYtr=meansYtr,
                  stdYtr=stdYtr)

  svd_resultx<- svd(X)
  U <- svd_resultx$u
  d <- svd_resultx$d
  V <- svd_resultx$v
  tV<-t(V)
  DUtY<- crossprod(U,Y) * d

  res_sqrt_inv_d2_lambda <- as.data.frame(outer(d, lambda,
                                                FUN = sqrt_inv_d2_lambda))
  names(res_sqrt_inv_d2_lambda)<- paste0("lambda",lambda)
  Bhat_comp <- furrr::future_map(res_sqrt_inv_d2_lambda,
                                 ~get_Bhat_comp(.x,
                                                DUtY = DUtY,
                                                nrank = nrank,
                                                tV = tV))
  Res <- list(Bhat_comp = Bhat_comp,
              rank = nrank,
              lambda=lambda,
              Scale=scalelist)

  if (!component){
    Bhat_mat <- furrr::future_map(Bhat_comp,
                                  ~Bhat_mat_rlist(.x, nrank = nrank))
    Res <- list(Bhat_mat = Bhat_mat,
                rank = nrank,
                lambda=lambda,
                Scale=scalelist)
  }
  return(Res)
}



#' @title Calculate the Bhat matrix from the return of the `rrda.fit` function.
#' @description This function calculates the coefficient Bhat (the coefficient of Ridge Redundancy Analysis for each parameter lambda and nrank) as a matrix form by using the Bhat components calculated by the rrda.fit function.

#' This function obtain the matrix form of Bhat as follows
#' \deqn{\hat{B}(\lambda, r) = FG^{\prime}}
#'
#' Here, the Bhat components F and G are obtained from the `rrda.fit` function as follows
#'
#' For \eqn{i = 1, \ldots, r}, the matrices \eqn{F} and \eqn{G} are defined as:
#' \deqn{F_{.i} = U_{\hat{B}(\lambda)}^{[i]}D_{\hat{B}(\lambda)}^{[i]}, \quad G_{.i} = V_{\hat{B}(\lambda)}^{[i]}}
#'
#' If the input already contains Bhat as matrix form (Bhat_mat), the function selects the preferred matrix from the list of Bhat matrices.
#'
#'  The function can handle different ranks (nrank) and ridge penalty values (lambda) based on the input. If nrank or lambda is NULL, the function will use the values from the Bhat components. Note that if lambda = NULL and B matrix is large (nrow(B)*ncol(B) > 100000), the function is performed for the minimum lambda value only.
#'
#' @param Bhat A list of vectors of Bhat components, obtained by the `rrda.fit` function. If Bhat_mat is detected in Input, it selects the preferred matrix from the list of Bhat.
#' @param nrank A numeric vector specifying the ranks of Bhat. Default is `NULL`, which sets it to the ranks defined in the Bhat components.
#' @param lambda A numeric vector of ridge penalty values. Default is `NULL`, which sets it to the lambda values defined in the Bhat components.
#' @return A list containing the Bhat matrix
#' @export
#' @importFrom furrr future_map
#' @examples
#' set.seed(10)
#' simdata<-rdasim1(n = 100,p = 200,q = 200,k = 5)
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' Bhat <- rrda.fit(Y = Y, X = X, nrank = c(1:10))
#' Bhat_mat <- rrda.coef(Bhat = Bhat, nrank = 10)
rrda.coef<- function(Bhat, nrank = NULL, lambda= NULL) {

  if(is.null(Bhat$rank)){
    stop("nrank is not found in Bhat \n")
  }

  if(is.null(Bhat$lambda)){
    stop("lambda is not found in Bhat \n")
  }

  if(is.null(nrank)){
    nrank <- Bhat$rank
  }

  lambdaNULL<- FALSE

  if(is.null(lambda)){
    lambdaNULL <- TRUE
    lambda <- Bhat$lambda
  }

  nrank_cd <- nrank %in% Bhat$rank

  if (all((nrank_cd))){
    # do nothing
  } else if (all(!(nrank_cd))){
    stop("nrank must be included in the Bhat \n")
  } else if (any((nrank_cd))){
    nrank_outvec<- paste(nrank[!(nrank_cd)],collapse=",")
    nrank <- nrank[(nrank_cd)]
    warning(paste0("nrank (= ",nrank_outvec,") not found in Bhat, this function is performed for the other values  \n"))
  }

  lambda_cd <- lambda %in% Bhat$lambda

  if (all((lambda_cd))){
    # do nothing
  } else if (all(!(lambda_cd))){
    stop("lambda must be included in the Bhat \n")
  } else if (any((lambda_cd))){
    lambda_outvec<- paste(lambda[!(lambda_cd)],collapse=",")
    lambda<-lambda[(lambda_cd)]
    warning(paste0("lambda (= ",lambda_outvec,") not found in Bhat, this function is performed for the other values \n"))
  }

  if (names(Bhat)[1] == "Bhat_mat") {
    message("Note: Bhat_mat is detected as the input of Bhat \n")

    if (all(Bhat$rank %in% nrank) && all(Bhat$lambda %in% lambda)) {
      return(Bhat$Bhat_mat)
    } else {
      mat <- Bhat$Bhat_mat
      lambda_indices <- which(names(mat) %in% paste0("lambda", lambda))

      matLR <- lapply(mat[lambda_indices], function(sub_mat) {
        sub_mat[names(sub_mat) %in% paste0("rank", nrank)]
      })

      Res <- matLR

      cat("rank:\n")
      print(nrank)
      cat("lambda:\n")
      print(lambda)

      return(Res)

    }

  } else if (names(Bhat)[1] == "Bhat_comp"){

  nrank<-sort(nrank)
  lambda<- sort(lambda)

  if(lambdaNULL){
    if (nrow(Bhat[[1]][[1]][[1]]) * nrow(Bhat[[1]][[1]][[2]]) > 100000){
      if(length(lambda) > 1){
        lambda<- lambda[1]
        message("Large B matrix; the smallest lambda is used, the other(s) are suppressed \n")
      }
    }
  }

  Bhat_lambdas<-Bhat[["Bhat_comp"]][Bhat$lambda %in% lambda]
  Bhat_mat<- furrr::future_map(Bhat_lambdas, ~Bhat_mat_rlist(.x, nrank = nrank))
  cat("rank:\n")
  print(nrank)
  cat("lambda:\n")
  print(lambda)

  return(Bhat_mat)


} else {
  stop("Bhat is not detected")
}

}


#' @title Compute the components of the coefficient Bhat using SVD.
#' @description Compute the components of Bhat using SVD. In our formula, Bhat is stored as LeftBhatlambda_k (p times r matrix) and RightBhatlambda_k (q times r). Here, n is the number of samples, p is the number of variables of X, q is the number of variables of Y, and r is the specified rank in the Ridge Redundancy Analysis.
#' @param rsid2_1L A numeric vector used for each lambda value.
#' @param DUtY A numeric matrix (n times q).
#' @param nrank A numeric vector indicating the rank(s) of Bhat.
#' @param tV A numeric matrix.
#' @return A list containing the left and right components of Bhat of (`LeftBhatlambda_k` and `RightBhatlambda_k`).
get_Bhat_comp<- function(rsid2_1L,DUtY,nrank,tV) {
  maxrank <- max(nrank)
  M   <- DUtY * rsid2_1L
  if (nrow(M) < 3 || ncol(M) < 3) {
    svd_M <- svd(x = M)
  } else {
    svd_M <- suppressWarnings(RSpectra::svds(A = M, k = maxrank))
  }

  BU <- svd_M$u[, 1:maxrank, drop = FALSE]
  Bd <- svd_M$d[1:maxrank]
  BV <- svd_M$v[, 1:maxrank, drop = FALSE]

  LeftBhatlambda_k <-  t(Bd * t(crossprod(tV, (BU * rsid2_1L))))
  RightBhatlambda_k <- BV
  # if the gsvd result D is included in Right..
  # LeftBhatlambda_k <- crossprod(tV, BU * rsid2_1L)
  # RightBhatlambda_k <- t(t(BV) * Bd)
  return(list(LeftBhatlambda_k=LeftBhatlambda_k,
              RightBhatlambda_k=RightBhatlambda_k))
}


#' @title Generate a list of rank-specific Bhat matrices (the coefficient of Ridge Redundancy Analysis for each parameter lambda and nrank).
#' @description Generate a list of rank-specific Bhat matrices (the coefficient of Ridge Redundancy Analysis for each parameter lambda and nrank). In our formula, Bhat is stored as LeftBhatlambda_k (p times r matrix) and RightBhatlambda_k (q times r). Here, n is the number of samples, p is the number of variables of X, q is the number of variables of Y, and r is the specified rank in the Ridge Redundancy Analysis.
#' @param Bhat_comp_1L A list containing components of Bhat for each lambda value.
#' @param nrank A numeric vector indicating the rank(s) of Bhat.
#' @return A list of matrices, each representing a rank-specific Bhat matrix.
Bhat_mat_rlist <- function(Bhat_comp_1L, nrank){
  LV <- Bhat_comp_1L$LeftBhatlambda_k
  RV<- Bhat_comp_1L$RightBhatlambda_k
  rlist<-get_rlist(LV=LV, RV=RV, nrank = nrank)
  rlist_nrank<-rlist[names(rlist) %in% nrank]
  names(rlist_nrank)<-paste0("rank",nrank)
  return(rlist_nrank)
}


#' @title Generate rank-specific matrices by combining the left and right components.
#' @description Generate rank-specific matrices by combining the left and right components of the coefficeint Bhat. In our formula, Bhat is stored as LeftBhatlambda_k (the left component vectors, p times r matrix) and RightBhatlambda_k (the right component vectors, q times r) for each lambda value. Here, n is the number of samples, p is the number of variables of X, q is the number of variables of Y, and r is the specified rank in the Ridge Redundancy Analysis.
#' @param LV A numeric matrix of the left component vectors.
#' @param RV A numeric matrix of the right component vectors.
#' @param nrank A numeric vector indicating the rank(s) of Bhat.
#' @return A list of matrices of rank-specific combinations of LV and RV.
get_rlist <- function(LV, RV, nrank) {
  maxrank <- max(nrank)
  rlist_each_rank <- furrr::future_map(seq_len(maxrank),
                                       ~outer(LV[, .x], RV[, .x], FUN = "*"))
  rlist <- Reduce(function(a, b) {a + b}, rlist_each_rank, accumulate = TRUE)
  names(rlist) <- c(seq_len(maxrank))
  return(rlist)
}



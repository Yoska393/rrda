
#' @title Calculate the predicted matrix Yhat using the coefficient Bhat obtained from the `rrda.fit` function.
#' @description This function calculates the predicted response matrix (Yhat) using Bhat (the coefficient of Ridge Redundancy Analysis) obtained from the rrda.fit function. The user can specify the ranks and lambda values to be used for the prediction. If not provided, the ranks and lambda values are set based on the Bhat components. The function returns the predicted Yhat matrix, as well as the ranks and lambda values used.
#'
#' The predicted response matrix is calculated as:
#' \deqn{\hat{Y}_r = X \left( \sum_{i=1}^{r} (F_{.i} G_{.i}^{\prime}) \right)
#' = \sum_{i=1}^{r} (X F_{.i} G_{.i}^{\prime})
#' = \sum_{i=1}^{r} (\tilde{X}_{.i} G_{.i}^{\prime})}
#'
#' Here, \eqn{\tilde{X} = X F}, and \eqn{r} is the rank of \eqn{\hat{B}(\lambda, r)}.
#' The regularized-rank-restricted estimation of \eqn{B} is obtained from the \code{rrda.fit} function:
#' \deqn{\hat{B}(\lambda, r) = FG^{\prime}}
#'
#' @param Bhat A list of vectors of Bhat components, obtained by the `rrda.fit` function.
#' @param X A numeric matrix of explanatory variables.
#' @param nrank A numeric vector specifying the ranks of Bhat. Default is `NULL`, which sets it to the ranks defined in the Bhat components.
#' @param lambda A numeric vector of ridge penalty values. Default is `NULL`, which sets it to the lambda values defined in the Bhat components.
#' @return A list containing the predicted Yhat matrix, ranks, and lambda values.
#' @export
#' @importFrom furrr future_map
#' @examples
#' set.seed(10)
#' simdata<-rdasim1(n = 100,p = 200,q = 200,k = 5)
#' X <- simdata$X
#' Y <- simdata$Y
#' Bhat <- rrda.fit(Y = Y, X = X, nrank = c(1:10))
#' Yhat_mat <- rrda.predict(Bhat = Bhat, X = X, nrank = 10)
rrda.predict <- function(Bhat, X, nrank = NULL, lambda = NULL) {

  if(is.null(nrank)){
    nrank <- Bhat$rank
  }
  if(is.null(lambda)){
    lambda <- Bhat$lambda
  }
  if(names(Bhat)[1] == "Bhat_mat"){
    stop("Bhat_mat is detected as the input of Bhat \n")
  }
  if (!is.numeric(nrank)){
    stop("nrank must be numeric")
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

  if (!is.numeric(lambda)){
    stop("lambda must be numeric") # put to other fuctions
  }
  lambda_cd <- lambda %in% Bhat$lambda
  if (all((lambda_cd))){
    # do nothing
  } else if (all(!(lambda_cd))){
    stop("lambda must be included in the Bhat \n")
  } else if (any((lambda_cd))){
    lambda_outvec<- paste(lambda[!(lambda_cd)],collapse=",")
    lambda<-lambda[(lambda_cd)]
    warning(paste0("lambda (= ",lambda_outvec,") not found in Bhat, this function is performed for the other values  \n"))
  }

  nrank<-sort(nrank)
  lambda<- sort(lambda)

  Bhat_lambdas<-Bhat[["Bhat_comp"]][Bhat$lambda %in% lambda]

  meansXtr <- NULL
  stdXtr <- NULL
  meansYtr <- NULL
  stdYtr <- NULL

  if (!is.null(Bhat$Scale$meansXtr)){
    meansXtr <- Bhat$Scale$meansXtr
    if (is.null(Bhat$Scale$stdXtr)){
      X <- scale(X, center = meansXtr, scale = FALSE)
    } else {
      stdXtr   <- Bhat$Scale$stdXtr
      X <- scale(X, center = meansXtr, scale = stdXtr)
    }
  }

  Yhat_mat<- furrr::future_map(Bhat_lambdas,
                               ~Yhat_mat_rlist(.x, nrank = nrank, X = X))

  if (!is.null(Bhat$Scale$meansYtr)){
    meansYtr <- Bhat$Scale$meansYtr
    if (is.null(Bhat$Scale$stdYtr)){
      Yhat_mat <- unscale_nested_matrices_map(mat=Yhat_mat,
                                              means=meansYtr)
    } else {
      stdYtr <- Bhat$Scale$stdYtr
      Yhat_mat <- unscale_nested_matrices_map(mat=Yhat_mat,
                                              means=meansYtr,
                                              std=stdYtr)
    }
  }

  Res <- list(Yhat_mat = Yhat_mat, rank = nrank, lambda=lambda)

  return(Res)
}


#' @title Generate a list of rank-specific Yhat matrices.
#' @description Generate a list of rank-specific Yhat matrices.
#' @param Bhat_comp_1L A list containing components of Bhat.
#' @param X A numeric matrix of explanatory variables.
#' @param nrank A numeric vector indicating the rank(s) of Bhat.
#' @return A list of matrices, each representing rank-specific predicted values Yhat.
Yhat_mat_rlist <- function(Bhat_comp_1L, X, nrank) {
  LV<- crossprod(t(X),Bhat_comp_1L$LeftBhatlambda_k)
  RV<- Bhat_comp_1L$RightBhatlambda_k
  rlist <- get_rlist(LV=LV, RV=RV, nrank = nrank)
  rlist_nrank<-rlist[names(rlist) %in% nrank]
  names(rlist_nrank)<-paste0("rank",nrank)
  return(rlist_nrank)
}



#' @title Apply unscaling to a nested list of matrices using specified mean and standard deviation values.
#' @description Apply unscaling to a nested list of matrices using specified mean and standard deviation values.
#' @param mat A nested list of matrices.
#' @param means A numeric vector of means for each matrix.
#' @param std A numeric vector of standard deviations for each matrix (optional).
#' @return A nested list of unscaled matrices.
unscale_nested_matrices_map <- function(mat, means, std=NULL) {
  furrr::future_map(mat, ~ furrr::future_map(.x,
                                             ~ unscale_matrices(.x,
                                                                means=means,
                                                                std=std)))
}




#' @title Unscale a matrix based on provided mean and standard deviation values.
#' @description Unscale a matrix based on provided mean and standard deviation values.
#' @param mat A numeric matrix to be unscaled.
#' @param means A numeric vector of means for each column.
#' @param std A numeric vector of standard deviations for each column (optional).
#' @return An unscaled numeric matrix.
unscale_matrices <- function(mat, means, std=NULL) {
  if(!(is.null(std))){
    mat <- sweep(mat, 2, std, FUN = "*")
  }
  mat<- sweep(mat, 2, means, FUN = "+")
  return(mat)
}


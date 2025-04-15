globalVariables(c("rank", "MSE","SEM","label"))

#' @title Cross-validation for Ridge Redundancy Analysis
#' @description This function performs cross-validation to evaluate the performance of Ridge Redundancy Analysis (RDA) models. It calculates the mean squared error (MSE) for different ranks and ridge penalty values through cross-validation folds. The function also supports centering and scaling of the input matrices.
#'
#' The range of lambda for the cross-validation is automatically calculated following the method of "glmnet" (Friedman et al., 2010). When we have a matrix of response variables (Y; n times q matrix) and a matrix of explanatory variables (X; n times p matrix), the largest lambda for the validation is obtained as follows
#'
#' \deqn{ \lambda_{\text{max}} = \frac{\max_{j \in \{1, 2, \dots, p\}} \sqrt{\sum_{k=1}^{q} \left( \sum_{i=1}^{n}  (x_{ij}\cdot y_{ik})  \right)^2}}{N \times 10^{-3}}}
#'
#' Then, we define \eqn{\lambda_{min}=10^{-4}\lambda_{max}}, and the sequence \eqn{\lambda} is generated based on the range.
#'
#' Also, to reduce the computation, the variable sampling is performed for the large matrix of X and Y (by default, when the number of the variables is over 1000). Alternatively, the range of lambda can be specified manually.
#' @param Y A numeric matrix of response variables.
#' @param X A numeric matrix of explanatory variables.
#' @param maxrank A numeric vector specifying the maximum rank of the coefficient Bhat. Default is `NULL`, which sets it to `(min(15, min(dim(X), dim(Y))))`.
#' @param lambda A numeric vector of ridge penalty values. Default is `NULL`, where the lambda values are automatically chosen.
#' @param nfold The number of folds for cross-validation. Default is 5.
#' @param folds A vector specifying the folds. Default is `NULL`, which randomly assigns folds.
#' @param sample.X A number of variables sampled from X for the lamdba range estimate. Default is 1000.
#' @param sample.Y A number of variables sampled from Y for the lamdba range estimate. Default is 1000.
#' @param center.X Logical indicating if `X` should be centered. If `TRUE`, scales `X`. Default is `TRUE`.
#' @param center.Y Logical indicating if `Y` should be centered. If `TRUE`, scales `Y`. Default is `TRUE`.
#' @param scale.X Logical indicating if `X` should be scaled. If `TRUE`, scales `X`. Default is `FALSE`.
#' @param scale.Y Logical indicating if `Y` should be scaled. If `TRUE`, scales `Y`. Default is `FALSE`.
#' @param verbose Logical indicating. If `TRUE`, the function displays information about the function call. Default is `TRUE`.
#' @return A list containing the cross-validated MSE matrix, lambda values, rank values, and the optimal lambda and rank.
#' @importFrom furrr future_map
#' @importFrom dplyr bind_cols
#' @importFrom RSpectra svds
#' @importFrom stats sd
#' @export
#' @examples
#'
#' set.seed(10)
#' simdata<-rdasim1(n = 10,p = 30,q = 30,k = 3)
#' X <- simdata$X
#' Y <- simdata$Y
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5)
#' rrda.summary(cv_result = cv_result)
#'
#' ##Complete Example##
#' set.seed(10)
#' simdata<-rdasim1(n = 10,p = 30,q = 30,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' rrda.summary(cv_result = cv_result) # cv result
#'
#' p <- rrda.plot(cv_result) # cv result plot
#' print(p)
#' h <- rrda.heatmap(cv_result) # cv result heatmao
#' print(h)
#'
#' estimated_lambda<-cv_result$opt_min$lambda  # selected parameter
#' estimated_rank<-cv_result$opt_min$rank # selected parameter
#'
#' Bhat <- rrda.fit(Y = Y, X = X, nrank = estimated_rank,lambda = estimated_lambda) # fitting
#' Bhat_mat<-rrda.coef(Bhat)
#' Yhat_mat <- rrda.predict(Bhat = Bhat, X = X) # prediction
#' Yhat<-Yhat_mat[[1]][[1]][[1]] # predicted values
#'
#' cor_Y_Yhat<-diag(cor(Y,Yhat)) # correlation
#' summary(cor_Y_Yhat)

rrda.cv <- function(Y, X, maxrank=NULL, lambda=NULL, nfold=5, folds = NULL, sample.X = 1000, sample.Y = 1000, scale.X = FALSE, scale.Y = FALSE, center.X = TRUE, center.Y = TRUE, verbose = TRUE){

  Y_cont <- deparse(substitute(Y))
  X_cont <- deparse(substitute(X))

  Y <- as.matrix(Y)
  X <- as.matrix(X)

  if (is.null(maxrank)){
    maxrank <- min(15,min(dim(X),dim(Y)))
  }
  if (maxrank > min(dim(X),dim(Y)) ){
    maxrank <- min(15,min(dim(X),dim(Y)))
    warning("rank(B) must be less than or equal to rank(X) and rank(Y) \n")
  }
  if (is.null(lambda)){
    lambda <- get_lambda(Y=Y,X=X,scale.X=scale.X, sample.X = sample.X, sample.Y = sample.Y)
  }
  if (any(lambda < 0)){
    stop("lambdas should be non-negative")
  }
  if(is.null(folds)){
    folds   <- sample(rep(seq_len(nfold), length.out = nrow(X)))
  }

  # Determine the smallest training set size across all folds
  min_n_all_folds <- min(sapply(seq_len(nfold), function(i) {
  	sum(folds != i)
  }))

  # Adjust maxrank if it's too large for the training set
  if (maxrank > min_n_all_folds) {
  	message(paste0("maxrank (", maxrank, ") is greater than the smallest training fold size (", min_n_all_folds, "). Adjusting maxrank to ", min_n_all_folds, "."))
  	maxrank <- min(maxrank, min_n_all_folds)
  }

  if (nfold > nrow(X)) {
  	message(paste0("nfold (", nfold, ") is greater than n (", nrow(X), "). Using nfold = ", nrow(X), " instead."))
  	nfold <- nrow(X)
  }

  nrank <- seq_len(maxrank)
  lambda <- sort(lambda)
  if (nfold  < 2){
    stop("nfold must be at least 2 for cross-validation.")
  } else {



    cv_MSE <- furrr::future_map(seq_len(nfold), function(i) {

      X_train <- X[folds != i,, drop=FALSE]
      X_test  <- X[folds == i,, drop=FALSE]
      Y_train <- Y[folds != i,, drop=FALSE]
      Y_test  <- Y[folds == i,, drop=FALSE]

      if (center.X){
        meansXtr <- colMeans(X_train)
        X_train <- scale(X_train, center = meansXtr, scale = FALSE)
        X_test  <- scale(X_test, center = meansXtr, scale = FALSE)

        if (scale.X){
        	constant_cols <- apply(X_train, 2, stats::sd) == 0
        	if (any(constant_cols)) {
        		message("Removing constant columns in X_train")
        		X_train <- X_train[, !constant_cols, drop=FALSE]
        		X_test <- X_test[, !constant_cols, drop=FALSE]
        	}

        	meansXtr <- colMeans(X_train)
          stdXtr   <- apply(X_train, 2, stats::sd)
          X_train <- scale(X_train, center = meansXtr, scale = stdXtr)
          X_test  <- scale(X_test, center = meansXtr, scale = stdXtr)

        }
      }else{
        if (scale.X){
          warning("X Scaling is not performed when center.X = FALSE")
        }
      }

      if (center.Y){
        meansYtr <- colMeans(Y_train)
        Y_train <- scale(Y_train, center = meansYtr, scale = FALSE)
        Y_test  <- scale(Y_test, center = meansYtr, scale = FALSE)
        if (scale.Y){
          constant_cols <- apply(Y_train, 2, stats::sd) == 0
          if (any(constant_cols)) {
            message("Removing constant columns in Y_train")
            Y_train <- Y_train[, !constant_cols, drop=FALSE]
            Y_test <- Y_test[, !constant_cols, drop=FALSE]
          }

          meansYtr <- colMeans(Y_train)
          stdYtr   <- apply(Y_train, 2, stats::sd)
          Y_train <- scale(Y_train, center = meansYtr, scale = stdYtr)
          Y_test  <- scale(Y_test, center = meansYtr, scale = stdYtr)

        }
      }else{
        if (scale.Y){
          warning("Y Scaling is not performed when center.Y = FALSE")
        }
      }

      Bhat <- rrda.fit(X = X_train, Y = Y_train,
                       nrank = nrank, lambda = lambda,
                       component = TRUE,
                       center.X = FALSE, center.Y = FALSE,
                       scale.X = FALSE, scale.Y = FALSE)

      res_MSE <- furrr::future_map(Bhat[[1]],
                                    ~MSE_lambda_rank(.x,
                                                      X = X_test,
                                                      Y = Y_test,
                                                      nrank =nrank))
      res_MSE_matrix <- dplyr::bind_cols(res_MSE)

      return(res_MSE_matrix)
    })
    cv_MSE_mean <- Reduce(`+`, cv_MSE) / nfold
    cv_MSE_sd <- sqrt(Reduce(`+`, lapply(cv_MSE, function(x) (x - cv_MSE_mean)^2)) / (nfold - 1))
    cv_MSE_se <- cv_MSE_sd / sqrt(nfold)

  }

  rownames(cv_MSE_mean) <- nrank
  colnames(cv_MSE_mean) <- lambda
  rownames(cv_MSE_se) <- nrank
  colnames(cv_MSE_se) <- lambda

  min_value <- min(cv_MSE_mean)
  min_position <- matrix(which(cv_MSE_mean == min_value, arr.ind = TRUE)[1,],1,2)
  row_name <- rownames(cv_MSE_mean)[min_position[1]]
  col_name <- colnames(cv_MSE_mean)[min_position[2]]

  bestrank <- nrank[min_position[1]]
  bestlambda <-lambda[min_position[2]]

  threshold <-  min_value + cv_MSE_se[min_position]

  lse_num <- max(which(cv_MSE_mean[min_position[1],] <= threshold))
  rse_num <- min(which(cv_MSE_mean[,min_position[2]] <= threshold))

  lambda.1se <- lambda[lse_num]
  rank.1se   <- nrank[rse_num]

  cv_lambda.1se<- cv_MSE_mean[min_position[1],lse_num]
  cv_rank.1se<- cv_MSE_mean[rse_num,min_position[2]]

  Res<- list(MSE = cv_MSE_mean,
             SEM = cv_MSE_se,
             rank = nrank,
             lambda = lambda,
             opt_min =  list(MSE = min_value,
                        lambda = bestlambda,
                        rank = bestrank),
             opt_lambda.1se = list(MSE = cv_lambda.1se,
                            lambda = lambda.1se,
                            rank = bestrank),
             opt_rank.1se = list(MSE = cv_rank.1se,
                            lambda = bestlambda,
                            rank = rank.1se)
             )

  if (verbose){

  cat("Call:\n")
  cat(paste0("rrda.cv(Y = ", Y_cont,
             ", X = ", X_cont,
             ", lambda = ",lambda[1]," - ",lambda[length(lambda)],
             ", maxrank = ", maxrank,", nfold = ", nfold,")\n"))

  }
  return(Res)
}


#' @title Estimate an appropriate value for the ridge penalty (lambda).
#' @description Estimate an appropriate value for the ridge penalty (lambda).
#' @param Y A numeric matrix of response variables.
#' @param X A numeric matrix of explanatory variables.
#' @param sample.X A number of variables sampled from X for the lamdba range estimate. Default is 1000.
#' @param sample.Y A number of variables sampled from Y for the lamdba range estimate. Default is 1000.
#' @param scale.X Logical indicating if `X` is scaled. Default is `FALSE`.
#' @return A numeric vector of the range of lambda values.
get_lambda <- function(Y, X, scale.X = FALSE, sample.X = sample.X, sample.Y = sample.Y) {
  Y<-as.matrix(Y)
  X<-as.matrix(X)
  y_rate<- 1
  if(ncol(X)>sample.X){
    X<-X[,sample(ncol(X), sample.X, replace = FALSE)]
  }
  if(ncol(Y)>sample.Y){
    y_rate<- ncol(Y)/sample.Y
    Y<-Y[,sample(ncol(Y), sample.Y, replace = FALSE)]
  }

  if (scale.X){
    X<-unbiased_scale(X)
    #X[is.nan(X)] <- 0
  }

  if(ncol(Y)==1){
    result <- max(abs(colSums(X*Y[,1])))/nrow(Y)
  } else{
    tmp_sum<- colSums((crossprod(Y,X))^2) * y_rate
    result <- max(sqrt(tmp_sum)) / (nrow(Y) * 10^(-3))
  }

  e<- 10^seq(-4,0, length.out = 50)
  lambda_def<- result * e
  lambda_def<- signif(lambda_def,4)
  return(lambda_def)
  }


#' @title Scale a matrix using unbiased estimators for the mean and standard deviation.
#' @description Scale a matrix using unbiased estimators for the mean and standard deviation.
#' @param x A numeric matrix to be scaled.
#' @importFrom stats sd
#' @return A scaled numeric matrix.
unbiased_scale <- function(x) {
  l <- nrow(x)
  xx<-scale(x) / sqrt((l-1)/l)
  xx[is.na(xx)] <- 0
  return(xx)
}

#' @title Compute MSE for different ranks of the coefficient Bhat and lambda.
#' @description For cross-validation, compute MSE for different ranks of the coefficient Bhat and lambda.
#' @param Bhat_comp_1L A list containing components of Bhat.
#' @param X A numeric matrix of explanatory variables.
#' @param Y A numeric matrix of response variables.
#' @param nrank A numeric vector indicating the rank(s) of Bhat.
#' @return A numeric vector of MSE values for each rank.
MSE_lambda_rank <- function(Bhat_comp_1L, X, Y, nrank) {
  rlist <- Yhat_mat_rlist(Bhat_comp_1L, X, nrank)
  MSE_ranks <- vapply(rlist, function(s) {
    rmse <- mean((Y - s)^2)
    return(rmse)
  }, numeric(1))
  return(MSE_ranks)
}

#' @title Summarize the results of cross-validation for the coefficient Bhat obtained from the `rrda.cv` function.
#' @description The function provides a summary of the results from cross-validation for Bhat obtained using the rrda.cv function. It displays the Mean Squared Error (MSE), rank, and lambda values for different cross-validation options.
#' @param cv_result A result list from the function `rrda.cv`, containing a matrix of MSE values for each rank and lambda, and a vector of lambda values.
#' @return Prints the estimated Coefficients, rank, lambda, and MSE from cross-validation.
#' @export
rrda.summary <- function(cv_result = cv_result) {
  cat("=== opt_min ===\n")
  cat("MSE: \n")
  print(cv_result$opt_min$MSE)
  cat("rank: \n")
  print(cv_result$opt_min$rank)
  cat("lambda: \n")
  print(cv_result$opt_min$lambda)

  cat("\n=== opt_lambda.1se ===\n")
  cat("MSE: \n")
  print(cv_result$opt_lambda.1se$MSE)
  cat("rank: \n")
  print(cv_result$opt_lambda.1se$rank)
  cat("lambda: \n")
  print(cv_result$opt_lambda.1se$lambda)

  cat("\n=== opt_rank.1se ===\n")
  cat("MSE: \n")
  print(cv_result$opt_rank.1se$MSE)
  cat("rank: \n")
  print(cv_result$opt_rank.1se$rank)
  cat("lambda: \n")
  print(cv_result$opt_rank.1se$lambda)
}

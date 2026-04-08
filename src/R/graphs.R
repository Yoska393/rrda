
#' @title Plot the results of cross-validation for Bhat obtained from the `rrda.cv` function.
#' @description This function visualizes the results of cross-validation for the estimated Bhat matrix obtained from the rrda.cv function. It creates a plot of the Mean Squared Error (MSE) for each combination of rank and lambda regularization parameter, allowing for the selection of specific ranks and lambda ranges to be plotted. Error bars representing the standard error of the MSE can be displayed for the best rank.
#' @param cv_result A result list from the function `rrda.cv`, containing a matrix of MSE values for each rank and lambda, and a vector of lambda values.
#' @param nrank A numeric vector specifying the ranks of Bhat to be plotted. Default is `NULL`, which plots all ranks.
#' @param min_l Minimum lambda value to be plotted. Default is `NULL`, which uses the minimum lambda value in `cv_result`.
#' @param max_l Maximum lambda value to be plotted. Default is `NULL`, which uses the maximum lambda value in `cv_result`.
#' @param show_error_bar Logical value indicating if the error bar is shown on the line that gives the best MSE value.
#' @param title Title of the figure
#' @return A plot of MSE cross-validation results.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom scales hue_pal
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#' set.seed(10)
#' simdata<-rdasim1(n = 100,p = 200,q = 200,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' rrda.summary(cv_result = cv_result)
#' rrda.plot(cv_result = cv_result)
#' }
rrda.plot <- function(cv_result, nrank = NULL, min_l = NULL, max_l = NULL, show_error_bar = FALSE, title = NULL) {
  res_MSE_matrix <- cv_result[["MSE"]]
  res_SEM_matrix <- cv_result[["SEM"]]  # Standard error matrix
  lambda <- cv_result[["lambda"]]
  bestrank <- cv_result$opt_min$rank

  if (is.null(min_l)) {
    min_l <- lambda[1]
  }
  if (is.null(max_l)) {
    max_l <- lambda[length(lambda)]
  }
  range <- min_l <= lambda & lambda <= max_l
  lambda <- lambda[range]
  res_MSE_matrix <- res_MSE_matrix[, range]
  res_SEM_matrix <- res_SEM_matrix[, range]  # Apply the same range for SEM

  if (is.null(nrank)) {
    nrank <- cv_result[["rank"]]
  }

  if (is.null(title)) {
   title = "MSE - Lambda and Rank"
  }

  res_MSE_matrix <- res_MSE_matrix[rownames(res_MSE_matrix) %in% nrank, ]
  res_SEM_matrix <- res_SEM_matrix[rownames(res_SEM_matrix) %in% nrank, ]  # Filter SEM

  # Convert the data to long format for ggplot
  res_MSE_df <- as.data.frame(res_MSE_matrix)
  res_MSE_df$rank <- rownames(res_MSE_matrix)

  res_SEM_df <- as.data.frame(res_SEM_matrix)
  res_SEM_df$rank <- rownames(res_SEM_matrix)

  data_long <- reshape2::melt(res_MSE_df,
                              id.vars = "rank",
                              variable.name = "lambda",
                              value.name = "MSE")
  data_long$lambda <- as.numeric(gsub("V", "", data_long$lambda))
  data_long$rank <- factor(data_long$rank,
                           levels = sort(as.numeric(unique(data_long$rank))))

  # SEM long format
  data_long_SEM <- reshape2::melt(res_SEM_df,
                                  id.vars = "rank",
                                  variable.name = "lambda",
                                  value.name = "SEM")
  data_long_SEM$lambda <- as.numeric(gsub("V", "", data_long_SEM$lambda))

  # Merge MSE and SEM
  data_long <- merge(data_long, data_long_SEM, by = c("rank", "lambda"))

  # Define a color palette with sufficient colors
  num_colors <- 10
  base_colors <- scales::hue_pal()(num_colors)
  rank_levels <- sort(as.numeric(unique(cv_result$rank)))
  color_palette <- stats::setNames(rep(base_colors,
                                       length.out = length(rank_levels)),
                                   rank_levels)

  # Define line types
  line_types <- stats::setNames(
    rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
        length.out = length(rank_levels),each=10),
    rank_levels)

  plot <- ggplot2::ggplot(data_long, ggplot2::aes(
    x = log(lambda),
    y = MSE,
    color = rank,
    linetype = rank)
  ) +

    ggplot2::geom_line() +
    ggplot2::labs(x = expression(log(lambda)),
                  y = "MSE",
                  title = title,
                  color = "Rank",     # Set legend title for color
                  linetype = "Rank") +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::scale_linetype_manual(values = line_types) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(size = 18),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if(min(lambda) <= cv_result$opt_min$lambda & cv_result$opt_min$lambda <= max(lambda)) {
    plot <- plot + ggplot2::geom_vline(xintercept = log(cv_result$opt_min$lambda),
                        linetype = "dotted",
                        color = "black", linewidth = 0.4)
  }

  if(min(lambda) <= cv_result$opt_lambda.1se$lambda & cv_result$opt_lambda.1se$lambda <= max(lambda)) {
    plot <- plot + ggplot2::geom_vline(xintercept = log(cv_result$opt_lambda.1se$lambda),
                        linetype = "dotted",
                        color = "black", linewidth  = 0.4)
  }

  # Add error bars for rank == 5
  if(show_error_bar){

  plot <- plot + ggplot2::geom_errorbar(
    data = subset(data_long, rank == bestrank),  # Only for rank 5
    ggplot2::aes(
      ymin = MSE - SEM,
      ymax = MSE + SEM
    ),
    width = 0.1,
    #color = color_palette[bestrank],
    color = "gray",
    linetype = "solid"
  )

  }

  return(plot)
}


#' @title Heatmap of the results of cross-validation for Bhat obtained from the `rrda.cv` function.
#' @description This function creates a heatmap to visualize the Mean Squared Error (MSE) results from the cross-validation of the Bhat matrix obtained from the rrda.cv function. The heatmap displays the MSE for different ranks of Bhat and values of the regularization parameter lambda, allowing users to visually assess the best combination of rank and lambda. The function also allows the user to highlight the points corresponding to the minimum MSE and the 1-standard error rule, helping to identify optimal model parameters.
#' @param cv_result A result list from the function `rrda.cv`, containing a matrix of MSE values for each rank and lambda, and a vector of lambda values.
#' @param nrank A numeric vector specifying the ranks of Bhat to be plotted. Default is `NULL`, which plots all ranks.
#' @param min_l Minimum lambda value to be plotted. Default is `NULL`, which uses the minimum lambda value in `cv_result`.
#' @param max_l Maximum lambda value to be plotted. Default is `NULL`, which uses the maximum lambda value in `cv_result`.
#' @param highlight_min Logical indicating if the marks should be plotted on the best prediction point, and 1se point. Default is `TRUE`.
#' @param title Title of the figure
#' @return A heatmap of MSE cross-validation results.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' set.seed(10)
#' simdata<-rdasim1(n = 100,p = 200,q = 200,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' rrda.summary(cv_result = cv_result)
#' rrda.heatmap(cv_result=cv_result)
#' }

rrda.heatmap <- function(cv_result, nrank = NULL, min_l = NULL, max_l = NULL, highlight_min = TRUE, title = NULL) {
  res_MSE_matrix <- cv_result[["MSE"]]
  lambda <- cv_result[["lambda"]]

  if (is.null(min_l)) {
    min_l <- lambda[1]
  }
  if (is.null(max_l)) {
    max_l <- lambda[length(lambda)]
  }
  range <- min_l <= lambda & lambda <= max_l
  lambda <- lambda[range]
  res_MSE_matrix <- res_MSE_matrix[, range]

  if (is.null(nrank)) {
    nrank <- cv_result[["rank"]]
  }

  if (is.null(title)) {
    title = "Heatmap of MSE - Lambda and Rank"
  }

  res_MSE_matrix <- res_MSE_matrix[rownames(res_MSE_matrix) %in% nrank, ]

  res_MSE_df <- as.data.frame(res_MSE_matrix)
  res_MSE_df$rank <- rownames(res_MSE_matrix)

  data_long <- reshape2::melt(
    res_MSE_df,
    id.vars = "rank",
    variable.name = "lambda",
    value.name = "MSE")
  data_long$lambda <- as.numeric(gsub("V", "", data_long$lambda))
  data_long$rank <- factor(data_long$rank,
                           levels = sort(as.numeric(unique(data_long$rank))))


  min_MSE <- min(data_long$MSE, na.rm = TRUE)
  min_MSE_point <- data_long[which.min(data_long$MSE), ]
  min_MSE_lse <- data.frame()
  min_MSE_rse <- data.frame()

  if (!is.null(cv_result$opt_lambda.1se$MSE)) {
    min_MSE_lse <- data_long[which(data_long$MSE == cv_result$opt_lambda.1se$MSE), ]
  }
  if (!is.null(cv_result$opt_rank.1se$MSE)) {
    min_MSE_rse <- data_long[which(data_long$MSE == cv_result$opt_rank.1se$MSE), ]
  }

  p <- ggplot2::ggplot(data_long, ggplot2::aes(x = log(lambda),
                                               y = rank,
                                               fill = MSE)) +
    ggplot2::geom_tile(color = "grey") +
    ggplot2::scale_fill_gradient(
      low = "darkblue",
      high = "white",
      limits = c(min(data_long$MSE), max(data_long$MSE) * 1.0)
    ) +
    ggplot2::labs(x = expression(log(lambda)), y = "Rank",
                  title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(size = 18),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 1
      ),
      panel.grid.major = ggplot2::element_line(linewidth = 0.5, color = "grey"),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.25, color = "grey")
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0))

  if (highlight_min) {

    min_MSE_point$label <- "Min"

    if (nrow(min_MSE_lse) > 0) {
      min_MSE_lse$label <- "Lambda 1se"
    }
    if (nrow(min_MSE_rse) > 0) {
      min_MSE_rse$label <- "Rank 1se"
    }

    points_df <- rbind(min_MSE_point, min_MSE_lse, min_MSE_rse)

    p <- p + ggplot2::geom_point(
      data = points_df,
      ggplot2::aes(x = log(lambda), y = rank, color = label, shape = label),
      size = 3)

    # Add legends for color and shape
    p <- p + ggplot2::scale_color_manual(
      name = "Highlights",
      values = c("Min" = "coral", "Lambda 1se" = "coral", "Rank 1se" = "coral"),
      breaks = c("Min", "Lambda 1se", "Rank 1se")) +
      ggplot2::scale_shape_manual(
        name = "Highlights",
        values = c("Min" = 18, "Lambda 1se" = 4, "Rank 1se" = 1),
        breaks = c("Min", "Lambda 1se", "Rank 1se"))
  }

  return(p)
}





#' @title Top feature interactions visualization with rank and lambda penalty
#' @description Visualizes the most influential feature interactions (based on the L2 norm) from Ridge Redundancy Analysis (RRDA) as a heatmap.
#'
#' Let the (rank-\eqn{r} truncated) decomposition of \eqn{\hat{B}(\lambda)} be
#' \deqn{\hat{B}(\lambda, r) = U_{\hat{B}(\lambda)} \, D_{\hat{B}(\lambda)} \, V_{\hat{B}(\lambda)}^{\prime}.}
#'
#' The following three biplot scalings are defined:
#'
#' **Symmetric scaling (default)**:
#' \deqn{\tilde{F} = U_{\hat{B}(\lambda)} \, D_{\hat{B}(\lambda)}^{1/2}, \qquad
#'       \tilde{G} = V_{\hat{B}(\lambda)} \, D_{\hat{B}(\lambda)}^{1/2}.}
#'
#' **X scaling**:
#' \deqn{\tilde{F} = U_{\hat{B}(\lambda)} \, D_{\hat{B}(\lambda)}, \qquad
#'       \tilde{G} = V_{\hat{B}(\lambda)}.}
#'
#' **Y scaling**:
#' \deqn{\tilde{F} = U_{\hat{B}(\lambda)}, \qquad
#'       \tilde{G} = V_{\hat{B}(\lambda)} \, D_{\hat{B}(\lambda)}.}
#'
#' In all three cases, \eqn{\hat{B}(\lambda, r) = \tilde{F} \, \tilde{G}^{\prime}.}
#'
#' Variable importance is scored by the row-wise \eqn{\ell_2}-norms:
#' \deqn{s_i^{(\tilde{F})} = \| \tilde{F}_{i,\cdot} \|_2, \qquad
#'       s_j^{(\tilde{G})} = \| \tilde{G}_{j,\cdot} \|_2.}
#'
#' Selecting the top \eqn{m_x} predictors and \eqn{m_y} responses yields the submatrices of the scaled factor matrices (each with \eqn{r} columns).
#'
#' The reduced coefficient submatrix is then
#' \deqn{\hat{B}_{\mathrm{sub}}(\lambda, r) =
#'       \tilde{F}_{\mathrm{sub}} \, \tilde{G}_{\mathrm{sub}}^{\prime}.}
#'
#' The matrix \eqn{\hat{B}_{\mathrm{sub}}(\lambda, r)} retains the dominant low-rank structure and is visualized as a heatmap (with \eqn{m_x = m_y = 20} by default).
#'
#' @param Y A numeric matrix of response variables.
#' @param X A numeric matrix of explanatory variables.
#' @param nrank Integer rank \eqn{r} of \eqn{\hat{B}} to visualize. If \code{NULL} (default), it is set to \code{min(5, min(dim(X), dim(Y)))}.
#' @param lambda A numeric vector of ridge penalty values. If \code{NULL} (default), it is set to 1.
#' @param mx Integer; number of top \eqn{X}-features (predictors) to display. Defaults to \code{20}.
#' @param my Integer; number of top \eqn{Y}-features (responses) to display. Defaults to \code{20}.
#' @param scaling Character string specifying how to apply the singular values from the compositions of \eqn{\hat{B}(\lambda)} when constructing the biplot factors.
#'Options are:
#'   \code{"symmetric"} (default) distributes singular values evenly to both sides (balanced scaling),
#'   \code{"x"} applies them fully to the X (left) side,
#'   \code{"y"} applies them fully to the Y (right) side,
#'and \code{"none"} removes them (no singular value weighting).
#' @param title Figure title. If \code{TRUE} (default), a formatted title is used. If \code{FALSE} or \code{NULL}, no title is drawn. If a single string, it is passed through to the figure title.
#' @return A list with elements: \code{heatmap} (pheatmap object), \code{B_sub} (mx x my matrix), \code{top_x}, \code{top_y}, \code{b1_sub}, \code{b2_sub}, \code{fit}, \code{scaling}.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' set.seed(10)
#' simdata<-rdasim1(n = 10,p = 50,q = 50,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#' rrda.top(Y=Y,X=X,nrank=5,lambda=1,mx=20,my=20)
#'
#' \dontrun{
#' ### In practice, the parameters nrank and lambda should be selected by CV ###
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' best_lambda<-cv_result$opt_min$lambda
#' best_rank<-cv_result$opt_min$rank
#' rrda.summary(cv_result = cv_result)
#'
#' rrda.top(Y=Y,X=X,nrank=best_rank,lambda=best_lambda,mx=20,my=20)
#' }
#'

rrda.top<-function(Y,X,nrank=NULL,lambda=NULL,mx=20,my=20,scaling = c("symmetric","none","x","y"),title=TRUE){

  scaling <- match.arg(scaling)

  if (!is.numeric(mx) || mx < 1) mx <- 1
  if (!is.numeric(my) || my < 1) my <- 1

  if (mx > ncol(X)) {
    warning("mx is larger than ncol(X). Setting mx = ncol(X).")
    mx <- ncol(X)
  }
  if (my > ncol(Y)) {
    warning("my is larger than ncol(Y). Setting my = ncol(Y).")
    my <- ncol(Y)
  }


  if (is.null(nrank)) {
    nrank <- min(5, min(dim(X), dim(Y)))
    message("nrank is set to default: ", nrank)
  }

  if (is.null(lambda)) {
    lambda <- 1
    message("lambda is set to default: ", lambda)
  }

  if (!is.numeric(nrank)) {
    stop("nrank must be numeric \n")
  }
  else if (max(nrank) > min(dim(X), dim(Y))) {
    stop("rank(B) must be less than or equal to rank(X) and rank(Y) \n")
  }
  if (any(lambda < 0)) {
    stop("lambdas should be non-negative")
  }

  # colnames for X and Y
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("Y", seq_len(ncol(Y)))
  }

  K <- as.integer(nrank)
  L <- lambda

  # Fit the model using the indicated lambda and rank
  b <- rrda.fit(Y = Y, X = X, lambda = L, nrank = K)

  b1<-b$Bhat_comp[[1]][[1]] #UD
  b2<-b$Bhat_comp[[1]][[2]] #V
  D <- b$Bhat_comp[[1]][[3]] #D

  Dvec <- if (is.matrix(D)) diag(D) else as.numeric(D)
  Dvec <- ifelse(Dvec > 0, Dvec, .Machine$double.eps)

  if (scaling == "none") {
    b1 <- sweep(b1, 2, Dvec, "/")
    #b2 = b2
  } else if (scaling == "x") {
    #b1 = b1
    #b2 = b2
  } else if (scaling == "y") {
    b1 <- sweep(b1, 2, Dvec, "/")
    b2 <- sweep(b2, 2, Dvec, "*")
  } else if (scaling == "symmetric") {
    sD <- sqrt(Dvec)
    b1 <- sweep(b1, 2, sD, "/")
    b2 <- sweep(b2, 2, sD, "*")
  }

  rownames(b1)<-colnames(X)
  rownames(b2)<-colnames(Y)

  x_scores <- apply(b1, 1, function(row) sqrt(sum(row^2, na.rm = TRUE)))
  top_x_idx <- order(x_scores, decreasing = TRUE)[1:mx]  # so that it becomes 30 variables afetr filtering

  y_scores <- apply(b2, 1, function(row) sqrt(sum(row^2, na.rm = TRUE)))
  top_y_idx <- order(y_scores, decreasing = TRUE)[1:my] # so that it becomes 30 variables afetr filtering

  b1_sub <- b1[top_x_idx, , drop = FALSE]
  b2_sub <- b2[top_y_idx, , drop = FALSE]

  B_sub <- b1_sub %*% t(b2_sub)

  h<-B_sub
  max_value <- max(h, na.rm = TRUE)
  min_value <- min(h, na.rm = TRUE)

  custom_colors <- colorRampPalette(c("blue", "white", "red"))(200)
  max_abs_value <- max(abs(c(min_value, max_value)))
  breaks <- seq(-max_abs_value, max_abs_value, length.out = 201)

  # --- title handling (no spurious warnings, F/NULL => no title) ---
  if (is.null(title)) {
    main_title <- NA_character_  # no title
  } else if (is.character(title) && length(title) == 1) {
    main_title <- title          # custom title
  } else {
    val <- suppressWarnings(as.logical(title))
    if (isTRUE(val)) {
      main_title <- sprintf("RRDA (rank = %s, lambda = %s)",
                            K, format(L, digits = 6, trim = TRUE))
    } else if (isFALSE(val)) {
      main_title <- NA_character_  # no title
    } else {
      warning("`title` must be TRUE, FALSE, NULL, or a single string. Falling back to default.")
      main_title <- sprintf("RRDA (rank = %s, lambda = %s)",
                            K, format(L, digits = 6, trim = TRUE))
    }
  }

  cluster_rows <- nrow(h) >= 2
  cluster_cols <- ncol(h) >= 2


  hp<-pheatmap::pheatmap(h,
               color = custom_colors,
               breaks = breaks,
               border_color = NA,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               show_rownames = TRUE,
               show_colnames = TRUE,
               cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,
               fontsize_row = 8,
               fontsize_col =8,
               width = 10,
               height = 15,
               main = main_title
  )


  print(hp)

  invisible(list(
    heatmap = hp,
    B_sub = B_sub,
    top_x = rownames(b1_sub),
    top_y = rownames(b2_sub),
    b1_sub = b1_sub,
    b2_sub = b2_sub,
    fit = b,
    scaling = scaling
  ))


}



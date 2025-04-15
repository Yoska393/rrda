
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
#' set.seed(10)
#' simdata<-rdasim1(n = 10,p = 30,q = 30,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' rrda.summary(cv_result = cv_result)
#' rrda.plot(cv_result = cv_result)
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
#' set.seed(10)
#' simdata<-rdasim1(n = 10,p = 30,q = 30,k = 3) # data generation
#' X <- simdata$X
#' Y <- simdata$Y
#'
#' cv_result<- rrda.cv(Y = Y, X = X, maxrank = 5, nfold = 5) # cv
#' rrda.summary(cv_result = cv_result)
#' rrda.heatmap(cv_result=cv_result)

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


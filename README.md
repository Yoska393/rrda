

<img src="rrda.jpg" width="500" height="400">

# Ridge Redundancy Analysis (`rrda`)

## 📖 Please cite:

**Yoshioka, H., Aubert, J., Iwata, H., and Mary-Huard, T. (2025).**  *Ridge Redundancy Analysis for High-Dimensional Omics Data.*

Could you cite when using `rrda` ? Thank you 😊

## 👀 Overview

Hello / Bonjour! 👋

This R package `rrda` provides functions for performing **ridge redundancy analysis (rrda)**, which is useful for modeling the relationship between a matrix of response variables (**Y**; n × q ) and a matrix of explanatory variables (**X**;  n × p ). The method is designed to handle **high-dimensional data efficiently**, allowing computation and storage optimization.


Also, I store all the scripts and my own functions used in our article (Yoshioka et al. 2025). 

## 🛠 Installation

You can install the package from GitHub using the `devtools` package:

```r
devtools::install_github("Yoska393/rrda", dependencies = TRUE)
```

## 📦 Dependencies

- `RSpectra`
- `furrr`
- `dplyr`
- `reshape2`
- `ggplot2`

## 📊 Data

For the application data, the data for breast cancer (Witten et al., 2009) and soybean data (Dang et al., 2023) are stored on this page. Also, refer to Ruiz-Arenas and González 2020) for the methylation data. 

The application data of breast cancer and soybean are stored as .rds file in a folder (RDAdata). For methylation data, you can refer to the R package (MEAL, Ruiz-Arenas and González 2024)

The metabolome data were downloaded from the RIKEN DropMet website (http://prime.psc.riken.jp/menta.cgi/prime/drop_index) ; ID: DM0071, DM0072.

## 🔧 Tutorial

#### Example 1 Fitting:

`rdasim1` function generates rank-restricted matrices X and Y. 

```r
simdata <- rdasim1(n = 100, p = 200, q = 200, k = 5)
X <- simdata×X
Y <- simdata×Y
```
`rrda.fit` function solves the rrda (ridge redundancy) for X and Y. This is equivalent to the prediction from X to Y, where Y = XB + E.

```r
# Sequential
Bhat <- rrda.fit(Y = Y, X = X, nrank = c(1:10))
names(Bhat)
```


#### Example 2 Parameter Tuning by Cross-Validation:
```r
set.seed(123)
simdata <- rdasim1(n = 100, p = 200, q = 200, k = 5)
X <- simdata×X
Y <- simdata×Y

cv_result <- rrda.cv(Y = Y, X = X)
rrda.summary(cv_result = cv_result)

# Complete Example
cv_result <- rrda.cv(Y = Y, X = X, maxrank = 10) # cv
rrda.summary(cv_result = cv_result) # cv result

# Plot the CV result
p <- rrda.plot(cv_result)
print(p)

# Heatmap of the CV result
h <- rrda.heatmap(cv_result)
print(h)

# Extract optimal parameters
estimated_lambda <- cv_result×opt_min×lambda
estimated_rank <- cv_result×opt_min×rank

# Fit the model with the optimal parameters
Bhat <- rrda.fit(Y = Y, X = X, nrank = estimated_rank, lambda = estimated_lambda)
Bhat_mat <- rrda.coef(Bhat)

# Make predictions
Yhat_mat <- rrda.predict(Bhat = Bhat, X = X)
Yhat <- Yhat_mat[[1]][[1]][[1]]

# Correlation
cor_Y_Yhat <- diag(cor(Y, Yhat))
summary(cor_Y_Yhat)
```

## 📚 References
- Dang, T., Fuji, Y., Kumaishi, K., Usui, E., Kobori, S., Sato, T., Toda, Y., Sakurai, K., Yamasaki, Y., Tsujimoto, H. and Hirai, M.Y., 2023. An integrative framework of stochastic variational variable selection for joint analysis of multi-omics microbiome data. bioRxiv, pp.2023-08.
- Ruiz-Arenas C, Gonzalez J (2024). MEAL: Perform methylation analysis. R package version 1.34.0. MEAL: Perform methylation analysis. 
- Witten D, Tibshirani R, Gross S, Narasimhan B (2024). PMA: Penalized Multivariate Analysis. R package version 1.2-4,
- Yoshioka, H., Aubert, J., Iwata, H., and Mary-Huard, T., 2025. Ridge Redundancy Analysis for High-Dimensional Omics Data.

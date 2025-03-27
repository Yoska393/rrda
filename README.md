

<img src="rrda.jpg" width="500" height="400">




# Ridge Redundancy Analysis (RRDA)

## Overview

Hello Bonjour,

Here I store all the scripts for the simulations and applications, and my own functions used for the analysis. Also, all the application data is stored as .rds file in a folder (RDAdata)

This package provides functions for performing **Ridge Redundancy Analysis (RRDA)**, which is useful for modeling the relationship between a matrix of response variables (**Y**; \( n 	\times q \)) and a matrix of explanatory variables (**X**; \( n 	imes p \)). The method is designed to handle **high-dimensional data efficiently**, allowing computation and storage optimization.

## Ridge Redundancy Analysis Model
The model is represented as:

\[
Y = XB + E
\]

where:
- \( Y \) is the response matrix (\( n 	imes q \))
- \( X \) is the predictor matrix (\( n 	imes p \))
- \( B \) is the regression coefficient matrix (\( p 	imes q \))
- \( E \) is the error matrix (\( n 	imes q \))

The regularized estimate of \( B \) is:

\[
\hat{B}(\lambda) = \left(X'X + \lambda P_{X'}
ight)^{-} X'Y
\]

Additionally, the regularized-rank-restricted estimation of \( B \) is:

\[
\hat{B}(\lambda, r) = U_{\hat{B}(\lambda)}^{[r]} D_{\hat{B}(\lambda)}^{[r]} V_{\hat{B}(\lambda)}^{[r]'}
\]

where:
- \( U_{\hat{B}(\lambda)}^{[r]} \) is a \( p 	imes r \) matrix
- \( D_{\hat{B}(\lambda)}^{[r]} \) is a \( r 	imes r \) diagonal matrix
- \( V_{\hat{B}(\lambda)}^{[r]} \) is a \( q 	imes r \) matrix

### Storing Large Coefficients
By default, the estimated coefficient matrix \( \hat{B}(\lambda) \) is stored as:

- **Left component** \( F \) (\( p 	imes r \))
- **Right component** \( G \) (\( q 	imes r \))

For \( i = 1, \dots, r \), the decomposition follows:

\[
F_{.i} = U_{\hat{B}(\lambda)}^{[i]}D_{\hat{B}(\lambda)}^{[i]}, \quad G_{.i} = V_{\hat{B}(\lambda)}^{[i]}
\]

To reconstruct \( \hat{B}(\lambda) \), use the function `rrda.coef()`.

## Functions

### `rrda.fit()`
#### Description:
Fits the Ridge Redundancy Analysis model.

#### Parameters:
- `Y`: Numeric matrix of response variables.
- `X`: Numeric matrix of explanatory variables.
- `nrank`: Vector specifying the rank(s) of \( \hat{B} \). Default: `NULL` (set to `1:min(15, min(dim(X), dim(Y)))`).
- `lambda`: Numeric vector of ridge penalty values. Default: `1`.
- `component`: Logical; if `TRUE`, returns \( \hat{B} \) as component vectors. Default: `TRUE`.
- `center.X`, `center.Y`: Logical; whether to center `X` and `Y`. Default: `TRUE`.
- `scale.X`, `scale.Y`: Logical; whether to scale `X` and `Y`. Default: `FALSE`.

#### Return:
A list containing \( \hat{B} \) as component matrices or full matrices.

#### Example:
```r
simdata <- rdasim1(n = 100, p = 200, q = 200, k = 5)
X <- simdata$X
Y <- simdata$Y

Bhat <- rrda.fit(Y = Y, X = X, nrank = c(1:10))
names(Bhat)
```

---

### `rrda.cv()`
#### Description:
Performs **cross-validation** to evaluate different **ranks** and **ridge penalty values** by computing Mean Squared Error (MSE).

#### Lambda Calculation:
By default, the range of \( \lambda \) is set automatically:

\[
\lambda_{	ext{max}} = rac{\max_{j \in \{1, 2, \dots, p\}} \sqrt{\sum_{k=1}^{q} \left( \sum_{i=1}^{n}  (x_{ij}\cdot y_{ik})  
ight)^2}}{N 	imes 10^{-3}}
\]

and:

\[
\lambda_{\min} = 10^{-4} \lambda_{\max}
\]

#### Parameters:
- `Y`: Numeric response matrix.
- `X`: Numeric predictor matrix.
- `maxrank`: Maximum rank of \( \hat{B} \). Default: `NULL` (set to `min(15, min(dim(X), dim(Y)))`).
- `lambda`: Ridge penalty values. Default: `NULL` (automatically set).
- `nfold`: Number of folds in cross-validation. Default: `10`.
- `folds`: Custom fold assignments. Default: `NULL` (random).
- `sample.X`, `sample.Y`: Number of variables sampled for lambda estimation. Default: `1000`.
- `center.X`, `center.Y`, `scale.X`, `scale.Y`: Centering/scaling options.

#### Return:
A list containing:
- Cross-validation **MSE matrix**
- **Optimal lambda** and **rank** values
- **Lambda sequence** and **rank sequence**

#### Example:
```r
simdata <- rdasim1(n = 100, p = 200, q = 200, k = 5)
X <- simdata$X
Y <- simdata$Y

cv_result <- rrda.cv(Y = Y, X = X)
rrda.summary(cv_result = cv_result)
```

---

## Example Workflow
```r
# Simulated Data
simdata <- rdasim1(n = 100, p = 200, q = 200, k = 5)
X <- simdata$X
Y <- simdata$Y

# Cross-validation
cv_result <- rrda.cv(Y = Y, X = X, maxrank = 10)
rrda.summary(cv_result)

# Plot Results
p <- rrda.plot(cv_result)
print(p)

h <- rrda.heatmap(cv_result)
print(h)

# Extract Optimal Parameters
estimated_lambda <- cv_result$opt_min$lambda
estimated_rank <- cv_result$opt_min$rank

# Fit Model with Optimal Parameters
Bhat <- rrda.fit(Y = Y, X = X, nrank = estimated_rank, lambda = estimated_lambda)

# Get Regression Coefficients
Bhat_mat <- rrda.coef(Bhat)

# Predict Y
Yhat_mat <- rrda.predict(Bhat = Bhat, X = X)
Yhat <- Yhat_mat[[1]][[1]][[1]]

# Compute Correlation
cor_Y_Yhat <- diag(cor(Y, Yhat))
summary(cor_Y_Yhat)
```

## Dependencies
The package depends on:
- `stats`
- `furrr`
- `RSpectra`
- `dplyr`

## Installation
To install this package, use:

```r
devtools::install_github("your-repo/RRDA")
```

## License
This package is open-source and distributed under the MIT license.

---
*For more details, refer to Yoshioka et al. (2025).*



Merci beaucoup, Hayato 24/Mar/2024

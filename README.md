# rrda: Ridge Redundancy Analysis for High-Dimensional Omics Data <img src="src/man/figures/logo.png" align="right" width="170"/>

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rrda)](https://cran.r-project.org/package=rrda)
[![CRAN Downloads (Last Month)](https://cranlogs.r-pkg.org/badges/last-month/rrda)](https://cran.r-project.org/package=rrda)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rrda)](https://cran.r-project.org/package=rrda)

<!-- badges: end -->

> The R package `rrda` provides functions for performing **ridge redundancy analysis (rrda)** 
> for high-dimensional datasets. It is useful for modeling the relationship between a matrix
>  of response variables (**Y**; n × q ) and a matrix of explanatory variables (**X**;  n × p )
>  with ridge penalty and rank restraint. The method is designed to handle **high-dimensional
> data**, allowing efficient computation and storage optimization.


## Overview
Hello / Bonjour / Konichiwa

- rrda/script_rrda:  Scripts and my own functions used in our article.
- rrda/RDAdata:  Application data used in our article.
- rrda/src:  Source code for the package "rrda"

**The link to our article** -> https://doi.org/10.1101/2025.04.16.649138

## 🛠 Installation

You can install the package from CRAN.

```{}
install.packages("rrda")
```

## 💡 Tutorial

```{r}
# rrda is updated if the version is old
required_version <- "0.2.3"  

if (!requireNamespace("rrda", quietly = TRUE) ||
    packageVersion("rrda") < required_version) {
  message("rrda will be updated")
  install.packages("rrda", repos = "https://cloud.r-project.org", type = "source")
}
```


#### Example 1: Fitting

`rdasim1` function generates rank-restricted matrices X and Y. 

```{r}
library(rrda)
set.seed(10)
simdata<-rdasim1(n = 50,p = 100,q = 100,k = 5)
X <- simdata$X
Y <- simdata$Y
```
`rrda.fit` function solves the rrda (ridge redundancy) for X and Y. 
This is equivalent to the prediction from X to Y, where Y = XB + E.

`nrank` indicates the rank restrictions for the model. Here, it is the value of 1 to 5.

`lambda` indicates the ridge penalty for the model. Here, it is the value of 0.1, 1, 10.   

The model solves several ranks and lambdas efficiently. In the default setting, the model returns all the combinations of 15 ranks and 50 lambda grid.

```{r}
Bhat <- rrda.fit(Y = Y, X = X)
names(Bhat)
```

When you see the Bhat, you will see the list composed of each lambda. In each lambda value, you have the coefficient `B` according to each rank.

(Note! The Bhat is stored in a decomposed form. This is because the function is designed for high-dimensional settings.)

#### Example 2: Parameter Tuning by Cross-Validation

Here we illustrate the parameter tuning process (regularization path), which helps identify the optimal parameter for maximizing prediction accuracy from **X** to **Y**.

How do we know the best lambda and rank for the model??
-> Cross-validation by `rrda.cv` function

```{r}
cv_result<- rrda.cv(Y = Y, X = X)
rrda.summary(cv_result = cv_result)

p <- rrda.plot(cv_result) # cv result plot
print(p)
```


`rrda.summary` tells you the parameters suggested via CV. 

```
=== opt_min ===
MSE: 
[1] 3.179695
rank: 
[1] 5
lambda: 
[1] 22.43
```

Also, `rrda.plot` and `rrda.heatmap` show you the figures to select the parameters.

<div style="display: flex; align-items: center; gap: 10px;">
  <img src="image/path.png" width="420" >
</div>

```{r}

# Choose the best parameter sets which gives the minimum MSE

best_lambda<-cv_result$opt_min$lambda  # selected parameter
best_rank<-cv_result$opt_min$rank # selected parameter

# Fitting with the best parameters
Bhat <- rrda.fit(Y = Y, X = X, nrank = best_rank,lambda = best_lambda) 

# Prediction
Yhat_mat <- rrda.predict(Bhat = Bhat, X = X) 
Yhat<-Yhat_mat[[1]][[1]][[1]] # predicted values

plot(Yhat, Y)
abline(0, 1, col = "red") 
```

<div style="display: flex; align-items: center; gap: 10px;">
  <img src="image/fit.png" width="500" >
</div>

#### Visualize and Select the Best Parameter

If you want to plot X and Y matrix in two-dimensional space (like classic RDA approach) :

```{r}
ud<-Bhat$Bhat_comp[[1]][[1]] # SVD component of B (UD) for lambda=0.1
v <-Bhat$Bhat_comp[[1]][[2]] # SVD component of B (V). for lambda=0.1

ud12 <- ud[, 1:2]
v12  <- v[, 1:2]

# Base plot: ud (e.g., site scores)
plot(v12, 
     xlab = "RRDA1", ylab = "RRDA2", 
     xlim = range(c(ud12[,1], v12[,1])) * 1.1, 
     ylim = range(c(ud12[,2], v12[,2])) * 1.1, 
     pch = 19, col = "darkgreen", 
     main = "RRDA")

# Add v (e.g., species scores) as arrows from origin
arrows(0, 0, ud12[,1], ud12[,2], col = "blue3", length = 0.1)

# Optionally add text labels
text(ud12, labels = paste0("X", 1:nrow(ud12)), pos = 3, col = "blue3", cex = 0.6)
text(v12, labels = paste0("Y", 1:nrow(v12)), pos = 3, col = "darkgreen", cex = 0.6)
```
<img src="image/rrda_2dim.png" width="500" >
However, this classic two-dimensional visualization is hard to interpret because of too many features..

For better interpretations, we visualize the **feature–feature matrix** using a selected dimensionality, highlighting the most informative features based on L2 norm.

```{r}
best_lambda<-cv_result$opt_min$lambda  
best_rank<-cv_result$opt_min$rank
rrda.top(Y=Y,X=X,nrank=best_rank,lambda=best_lambda,mx=20,my=20)
```
<img src="image/rrda_heat.png" width="500" >


## For more exercises with application data

Go to Rpubs (https://rpubs.com/Yoska393/1351133).


## 📊 Application Data
The application data of breast cancer and soybean are stored as .rds file in a folder `rrda/RDAdata`.

The TGCA methylation dataset is too large to be hosted on GitHub.  
Please refer to the script `rrda_script/Meth.rmd` for details on downloading and preprocessing.

In short, methylation and gene expression data are available via the **brgedata**:

```r
library(brgedata)

data(brge_methy)
brge_methy

data(brge_gexp)
brge_gexp
```


## 📚 References 

Please cite :)

- Yoshioka, H., Aubert, J., Iwata, H., and Mary-Huard, T., 2025. Ridge Redundancy Analysis for High-Dimensional Omics Data. *bioRxiv*, doi: 10.1101/2025.04.16.649138

- Yoshioka H, Aubert J, and Mary-Huard T (2025). rrda: Ridge Redundancy Analysis for High-Dimensional Omics Data.  https://CRAN.R-project.org/package=rrda (CRAN R Package) 

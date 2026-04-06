# RRDA Data

This repository contains datasets used in the article:  
**“Ridge Redundancy Analysis for High-Dimensional Omics Data” (2025)**  
https://doi.org/10.1101/2025.04.16.649138

## Data availability

The following datasets are included in this repository:

- **Soybean data**
  - `SoyData.RDS`
  - `soymet.csv`
  - `soymicro.csv`

- **Breast cancer data**
  - `breast.RDS`

## Methylation data

The TCGA methylation dataset is too large to be hosted on GitHub.  
Please refer to the script `rrda_script/Meth.rmd` for details on downloading and preprocessing.

The TCGA methylation and gene expression data, along with CpG annotation, can be obtained as follows:

```r
# install if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "brgedata",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
))

library(brgedata)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(brge_methy)
data(brge_gexp)
```

---

Merci beaucoup,  
Hayato  
19/11/2024

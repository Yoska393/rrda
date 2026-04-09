# 📊 Application Data

This repository contains datasets used in the article:  
**“Ridge Redundancy Analysis for High-Dimensional Omics Data” (2025)**  
https://doi.org/10.1101/2025.04.16.649138

## Data availability

The following datasets are included in this repository:

| Dataset | Files | Format | Location | Notes |
|--------|------|--------|----------|------|
| Soybean | `SoyData.RDS`, `soymet.csv`, `soymicro.csv` | RDS / CSV | `rrda/RDAdata` | Multi-omics dataset (genome, metabolome, microbiome) |
| Breast cancer | `breast.RDS` | RDS | `rrda/RDAdata` | Integrated omics dataset |
| TCGA methylation | Not included | — | — | See script `rrda_script/Meth.rmd` for download and preprocessing |

---

## TCGA data - how to access

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

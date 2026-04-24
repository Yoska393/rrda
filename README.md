# Reproducibility Repository: Ridge Redundancy Analysis for High-Dimensional Omics Data <img src="src/man/figures/logo.png" align="right" width="170"/>


> This repository provides the **datasets, analysis scripts, and supplementary materials** used in our study on Ridge Redundancy Analysis (rrda).  
> It is intended for reproducing the analyses in the paper and exploring practical applications.

---

## 🔗 Related Resources

| Resource | Link | Description |
|----------|------|-------------|
| 📦 R Package (CRAN) | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rrda)](https://cran.r-project.org/package=rrda) | Stable release |
| 🔧 Package Development (FORGE) | [forge.inrae.fr/mia-paris/rrda](https://forge.inrae.fr/mia-paris/rrda) | Source code, dev version, issue tracker |
| 📖 Package Documentation | [rrda-00f094.pages-forge.inrae.fr](https://rrda-00f094.pages-forge.inrae.fr) | Function reference & vignettes |
| 📄 Paper | [doi: 10.1101/2025.04.16.649138](https://doi.org/10.1101/2025.04.16.649138) | bioRxiv preprint |
| 🧪 Tutorial (RPubs) | [rpubs.com/Yoska393/1351133](https://rpubs.com/Yoska393/1351133) | Application exercises |

> **For package installation and usage**, please refer to the [FORGE repository](https://forge.inrae.fr/mia-paris/rrda) and its [documentation site](https://rrda-00f094.pages-forge.inrae.fr).

---

## 📂 Repository Structure

| Folder / File | Description |
|----------------|-------------|
| [`RDAdata/`](./RDAdata/) | Application datasets (soybean, breast cancer) |
| [`script_rrda/`](./script_rrda/) | Analysis workflows and preprocessing scripts |
| [`lecture/`](./lecture/) | Supplementary materials and lecture notes |
| [`image/`](./image/) | Figures used in documentation |

---

## 📊 Application Data

The application datasets used in this study are summarized below.

| Dataset | Files | Format | Location | Notes |
|--------|------|--------|----------|------|
| Soybean | `SoyData.RDS`, `soymet.csv`, `soymicro.csv` | RDS / CSV | [`RDAdata/`](./RDAdata/) | Multi-omics dataset (genome, metabolome, microbiome) |
| Breast cancer | `breast.RDS` | RDS | [`RDAdata/`](./RDAdata/) | Integrated omics dataset |
| TCGA methylation | Not included | — | Open Source | See [`script_rrda/Meth.rmd`](./script_rrda/Meth.rmd) for download instructions |

### Obtaining the TCGA Methylation Dataset

The TCGA methylation dataset is too large to be hosted on GitHub.  
The data along with CpG annotation can be obtained as follows:

```r
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

## 🔁 Reproducing the Analyses

Scripts in [`script_rrda/`](./script_rrda/) correspond to the analyses presented in the paper.  
To reproduce results, first install the `rrda` package from CRAN:

```r
install.packages("rrda")
```

Then run the scripts in [`script_rrda/`](./script_rrda/) with the data from [`RDAdata/`](./RDAdata/).  
Detailed instructions are provided within each script.

---

## 📚 Citation

If you use this repository or the `rrda` package in your research, please cite:

- Yoshioka, H., Aubert, J., Iwata, H., and Mary-Huard, T., 2025. Ridge Redundancy Analysis for High-Dimensional Omics Data. *bioRxiv*, doi: [10.1101/2025.04.16.649138](https://doi.org/10.1101/2025.04.16.649138)

- Yoshioka H, Aubert J, and Mary-Huard T (2025). *rrda: Ridge Redundancy Analysis for High-Dimensional Omics Data*. https://CRAN.R-project.org/package=rrda

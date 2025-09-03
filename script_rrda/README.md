## 📂 script_rrda  
All scripts for simulations, applications, and supporting functions used in the analysis.  

The file names correspond to the sections in our article (Yoshioka et al.).  

---

### 📂 Functions  

| File | Description |
|------|-------------|
| **functions_rrda_yoshioka.R** | Supporting original functions. |
| **functions_stars_rrr_yoshioka.R** | Supporting functions for `rrpack` (Chen et al. 2022) and `StARS` (Wen et al. 2024). |

---

### 📂 Simulation  

Scripts under the `Simulation` folder.  

| File | Description |
|------|-------------|
| **scenario1.Rmd** | Computational performance (sRDA, rrpack, and rrda). |
| **scenario2.Rmd** | Rank estimation. |
| **scenario3.Rmd** | Prediction performance. |
| **visualization.Rmd** | Visualization of simulation results. |
| **lambda_max_computation.Rmd** | Sampling of matrices to compute λ<sub>max</sub> (Section 2.4). |
| **rank_estimate_track.Rmd** | Detailed analysis of rank estimation for all criteria and StARS (Fig 2S). |

---

### 📂 Application  

Scripts under the `Application` folder.  

| File | Description |
|------|-------------|
| **Breastcancer.Rmd** | Breast cancer dataset analysis (Section 3.1). |
| **Breastcancer_each_chr.Rmd** | Breast cancer analysis by chromosome (Section 3.1). |
| **Breastcancer_each_chr_rev.Rmd** | Breast cancer analysis with reversed chromosome scenario (Section 3.1). |
| **Soybean.Rmd** | Soybean multi-omics data (Section 3.2). |
| **Meth.Rmd** | Methylation and expression data (Section 3.3). |
| **App_RandomSplit.Rmd** | Random Split analysis for all applications (Sections 3.1–3.3). |
| **App_RandomSplit_rrda_range.Rmd** | Random Split analysis with rrda λ range (Sections 3.1–3.3). |

---

### 📂 sRDA  

Scripts under the `sRDA` folder.  

| File | Description |
|------|-------------|
| **sRDA.Rmd** | Comparison of sRDA vs rrda on breast cancer data (Section 3). |

---
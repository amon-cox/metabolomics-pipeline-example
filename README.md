# metabolomics-pipeline-example
A basic analysis pipeline in R for untargeted LC-MS metabolomics, including quality control, normalization, PCA overview of sample behavior, and differential feature analysis with linear models.

### Case
The (anonymized and modified) data represent panels of untargeted LC–MS features from positive and negative ionization mode runs, such as metabolites measured from tissue, stool samples, or microbial cultures. In this hypothetical scenario, samples belong to control subjects, one of two treatments (tr1 or tr2), or a dual‑treatment group (tr1+tr2). The pipeline demonstrates a conventional approach to characterizing metabolomic differences among these groups.

### Methods Summary
- Preprocessing: merge duplicate LC-MS peaks, sample-wise median normalization, log2 transformation.
- Sample analysis: Principal Components Analysis (PCA), PCA-based outlier detection, and Permutational Analysis of Variance (PERMANOVA) for treatment effects.
- Feature analysis: linear models with cohort blocking and multiple contrasts (limma).
- Reproducibility: environment managed by renv.

### Repository Structure
- run_all.r
- R/
  - 00_setup.r
  - 01_preprocess.r
  - 02_sample_analysis.r
  - 03_feature_analysis.r
- data_raw/
- data_processed/
- output/
  - plots/
  - tables/

### How to Run
In R, from the project root:

```
install.packages("renv")
renv::restore()
source("run_all.r")
```

### Key Results
- PCA demonstrates clear separation of treatments along PC1 in both ion modes.
  - tr1 profiles appear indistinguishable from controls, whereas tr2 and tr1+tr2 profiles demonstrate marked contrasts from controls.
  - No discernible separation by cohort.
- PERMANOVA indicates significant differences amongst metabolomic profiles by treatment group.
- Pairwise PERMANOVA identifies significant differences which corroborate PCA observations.
- limma highlights over 1,000 features differentially abundant between tr2 and control groups.

### Upcoming
- Putative feature annotation using KEGG (example workflow).
- R Markdown or Quarto project report summarizing methods and results.
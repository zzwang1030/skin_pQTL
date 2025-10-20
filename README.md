# Skin pQTL Study

This repository contains analysis scripts for the **skin protein quantitative trait loci (pQTL)** study.  
The project integrates skin proteomics, genetics, and transcriptomics to uncover molecular mechanisms underlying leprosy susceptibility and skin immune regulation.

---

##  Overview

The scripts in this repository cover the complete analysis workflow:

1. **Proteome processing and QC**
2. **Age-associated proteome analysis**
3. **pQTL discovery and fine-mapping**
4. **QTL post-processing and annotation**
5. **Transcriptome-wide and proteome-wide association analyses (TWAS/PWAS)**
6. **Integration with disease GWAS (e.g., leprosy)**

---

##  File Description

| File | Description |
|------|--------------|
| `01_proteome_processing.R` | Preprocessing, normalization, and quality control of proteomics data. |
| `02_proteome_aging.R` | Statistical modeling of age-related changes in protein abundance. |
| `03_pQTL_formal_code.sh` | Main pipeline for pQTL mapping using genotype and proteome data (e.g., via QTLtools or TensorQTL). |
| `04_QTL_process.R` | Downstream processing: variant annotation, LD pruning, and visualization of QTL results. |
| `05_TWAS_fusion_code_leprosy.sh` | Transcriptome-wide association study (TWAS) and cross-trait integration with leprosy GWAS using FUSION. |
| `06_PWAS_process.R` | Proteome-wide association (PWAS) analysis and colocalization with disease loci. |

---

## ⚙️ Dependencies

- **R** (≥ 4.2.0) with packages:  
  `tidyverse`, `data.table`, `ggplot2`, `MatrixEQTL`, `plink2R`, `coloc`, `qvalue`
- **Python** (optional, for specific utilities)
- **QTLtools** or **TensorQTL** for pQTL mapping
- **FUSION** for TWAS
- **PLINK 1.9/2.0**
- **LocusZoom** or **ggplot2** for visualization

---

Modify input/output paths in each script as needed.

📄 License

This project is released under the MIT License.
You are free to reuse and adapt the code with appropriate citation.

📚 Citation
(Manuscript in preparation.)

📬 Contact

For questions or collaborations, please contact:
[Yuanqiang Sun] (sunyuanqiang[AT]pku.edu.cn)

# skin_pQTL Analysis

This repository contains the code and essential data used for the analysis in the skin pQTL study.  
The project includes R scripts, shell workflows, and  datasets.  

---

##  Overview

The scripts in this repository cover the complete analysis workflow:

1. **Proteome processing and QC**
2. **Age-associated proteome analysis**
3. **pQTL discovery and fine-mapping**
4. **QTL post-processing and annotation**
5. **Cross-tissue pQTL comparison**
6. **Integration with disease GWAS (e.g., leprosy)**


---

## Data Availability

The `data/` folder contains **essentail files** to illustrate input formats required to reproduce all results  .  
** Some full datasets are not completely included** due to size constraints (>50–100 MB).


To fully reproduce the results, please download or generate the required data.

Links for downloading full data or instructions to generate them may be provided in R or shell scripts.

---

## Usage

### 1. Running the R scripts
Navigate to the `R/` folder.  
Scripts are named following the analysis workflow (e.g., preprocessing → pQTL mapping → integration → visualization).

Example:
```bash
Rscript Fig1.R
Rscript Fig2.R
```


### 2. Running the shell scripts

Shell scripts in shell/ folder are used for:

batch processing

submitting jobs to HPC clusters

running QTL tools

Example:
```bash
bash 03_pQTL_formal_code.sh
```

Ensure required software and environment modules are installed.


---

## ⚙️ Dependencies

- **R** (≥ 4.2.0) with packages:  
  `susieR`, `data.table`, `ggplot2`, `BisqueRNA`, `plink`, `coloc`, `qvalue`
- **Python** (optional, for specific utilities)
- **QTLtools** or **FastQTL** for pQTL mapping
- **FUSION** for TWAS
- **PLINK 1.9/2.0**
- **LocusZoom** or **ggplot2** for visualization

---

Modify input/output paths in each script as needed.

## 📄 License

This project is released under the MIT License.
You are free to reuse and adapt the code with appropriate citation.

## 📚 Citation
(Manuscript in preparation.)

## 📬 Contact

For questions or collaborations, please contact:
[Yuanqiang Sun] (sunyuanqiang[AT]pku.edu.cn)

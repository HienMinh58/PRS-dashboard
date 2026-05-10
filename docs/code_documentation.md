# PRS Dashboard Technical Documentation

This document provides a technical overview of the PRS Dashboard codebase, detailing the architecture, core modules, and key functionalities.

## 1. Architecture Overview

The PRS Dashboard is a Streamlit-based web application designed for Polygenic Risk Score (PRS) analysis. The architecture follows a modular approach:

- **Frontend (`ui.py`, `app.py`)**: Handles user interaction, file uploads, and configuration.
- **Processing Pipeline (`prs_methods.py`, `qc.py`)**: Executes the heavy lifting, including genetic quality control and PRS calculation.
- **Validation Layer (`validation.py`)**: Ensures data integrity before processing.
- **Analysis Layer (`ml_models.py`, `utils.py`)**: Provides downstream analysis and performance evaluation.

---

## 2. Core Modules

### 🧬 `prs_methods.py`
The heart of the PRS calculation logic. It orchestrates the execution of various state-of-the-art PRS tools.

- **`run_prs_csx(...)`**: 
    - Executes the PRS-CSx pipeline.
    - Handles multi-ancestry GWAS data.
    - Performs per-chromosome processing and merges results.
    - Includes an optional validation step to fit meta-score weights using Linear/Logistic regression.
- **`execute_prs_pipeline(...)`**: 
    - The main entry point for running multiple methods (PRS-CSx, PROSPER, etc.) in sequence.
    - Manages progress bars and result merging.

### 🛡️ `qc.py`
Handles Quality Control for GWAS summary statistics and PLINK binary files.

- **`run_qc_v1(...)`**: 
    - Matches GWAS SNPs to target BIM files.
    - Filters by Minor Allele Frequency (MAF).
    - Removes ambiguous SNPs (A/T, C/G).
- **`clean_plink_invalid_alleles(...)`**: 
    - Uses PLINK2 to remove variants with invalid or duplicate alleles.
- **`standardize_gwas_columns(...)`**: 
    - Maps diverse GWAS column names (e.g., `MARKERNAME` -> `SNP`) to a standard internal format.

### ✅ `validation.py`
Ensures that uploaded files are correctly formatted and suitable for analysis.

- **`validate_gwas(...)`**: Validates GWAS statistics, cleans alleles, and saves a standardized version.
- **`validate_plink(...)`**: Checks for the existence of BED/BIM/FAM files and verifies basic SNP integrity.
- **`validate_ld_ref(...)`**: Verifies that the required LD reference panels for specific populations and chromosomes are present in the system.

### 🤖 `ml_models.py`
Provides a machine learning pipeline to evaluate the predictive power of the calculated PRS.

- **`train_ml_models(...)`**: 
    - Supports **SVM**, **GLM**, and **Random Forest**.
    - Handles both regression (continuous traits) and classification (binary traits).
    - Automatically scales features and splits data into training/testing sets.

---

## 3. Technical Stack

- **UI Framework**: Streamlit
- **Data Analysis**: Pandas, NumPy
- **Machine Learning**: Scikit-learn
- **Scientific Computing**: SciPy
- **Bioinformatics Tools**: PLINK2, PRS-CSx (Python tool), PROSPER (R script)

## 4. Directory Structure

```text
.
├── app.py                 # Main Streamlit application
├── src/                   # Source code
│   ├── ui.py              # UI components and sidebar
│   ├── prs_methods.py     # PRS algorithm wrappers
│   ├── qc.py              # Quality control logic
│   ├── validation.py      # Input validation
│   ├── ml_models.py       # ML evaluation pipeline
│   └── utils.py           # General helper functions
├── ld_reference/          # LD Reference storage (Local)
└── results/               # Output directory for scores
```

---
*Documentation generated for PRS-Dashboard developers.*

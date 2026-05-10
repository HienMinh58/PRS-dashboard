<p align="center">
  <img src="assets/icon_v2.png" width="120" alt="PRS Dashboard Icon">
</p>

# PRS Dashboard

A powerful, user-friendly web interface for Polygenic Risk Score (PRS) analysis and Machine Learning evaluation.

![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-1.20%2B-FF4B4B?logo=streamlit&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED?logo=docker&logoColor=white)
![PRS-CSx](https://img.shields.io/badge/Method-PRS--CSx-orange)
![Status](https://img.shields.io/badge/Status-Active-success)

## Overview

### The Problem
Polygenic Risk Score (PRS) analysis is a cornerstone of modern genomic research, yet the barrier to entry remains high. Most state-of-the-art tools require complex manual setups, including:
- Configuring intricate command-line environments.
- Manually managing GWAS summary statistics and large PLINK genotype files.
- Handling diverse LD reference panels.
- Performing multi-step Quality Control (QC) filtering and data alignment.

### The Solution
The **PRS Dashboard** simplifies this entire workflow into a streamlined, browser-based experience. Built with Streamlit and containerized with Docker, it provides researchers and students with a professional interface to configure, execute, and visualize PRS pipelines without the CLI headache.

## Features

- **📊 Comprehensive GWAS Support:** Easy upload and validation of GWAS summary statistics with automated parsing.
- **🧬 PLINK Integration:** Full support for target genotype files (`.bed`, `.bim`, `.fam`) with automatic SNP matching.
- **🌍 Multi-Ancestry Workflows:** Support for both single-ancestry and complex multi-ancestry (PRS-CSx) configurations.
- **🛡️ Quality Control Pipeline:** Automated QC filtering, including ambiguous SNP removal and MAF thresholding.
- **🧪 Validation & Evaluation:** Optional upload for validation phenotypes and covariates to test score performance.
- **🤖 ML Evaluation Pipeline:** Integrated Machine Learning module (SVM, Random Forest, Regression) to evaluate PRS predictive power.
- **🐳 Reproducible Setup:** Entirely Docker-based environment ensuring consistency across different operating systems.

## 🚀 Recent Optimizations & Updates

- **⚡ High-Performance MCMC:** Optimized core PRS-CSx MCMC iterations by replacing standard Python sum functions with NumPy vectorized operations (`np.sum`), significantly reducing computation time during posterior estimation.
- **🧬 Target-Based LD Computation (Beta):** Integrated local LD computation directly from user-provided target genotype data using PLINK2, providing more flexibility beyond pre-specified reference panels.
- **🛠️ Development Workflow:** Improved Docker configuration with local tool mounting for seamless development and live-code optimization of core algorithms.

## Demo

<!-- Add a screenshot or demo GIF here -->
![PRS Dashboard Demo](assets/demo.gif)

> [!TIP]
> **Video Walkthrough:** Check out our [video demonstration](https://github.com/user-attachments/assets/5c0c26ff-106b-4961-8565-99df306ef613) to see the dashboard in action!

---

## Getting Started

### Prerequisites
- [Docker](https://www.docker.com/get-started) and [Docker Compose](https://docs.docker.com/compose/install/)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/HienMinh58/PRS-dashboard.git
   cd PRS-dashboard
   ```

2. **Prepare Reference Data**
   Place your LD reference panels in the `ld_reference/` directory (e.g., 1000 Genomes HDF5 files).

3. **Launch with Docker**
   ```bash
   docker-compose up --build
   ```

4. **Access the Dashboard**
   Open [http://localhost:8501](http://localhost:8501) in your browser.

## Project Structure

```text
.
├── app.py                 # Streamlit application entry point
├── src/                   # Core logic (QC, PRS methods, ML)
├── ld_reference/          # LD Reference panels (Local storage)
├── data/                  # Input datasets (GWAS, Genotypes)
├── results/               # Output scores and reports
└── Dockerfile             # Container configuration
```

---
*Developed for research and educational purposes. Designed for bioinformatics students and genomic researchers.*

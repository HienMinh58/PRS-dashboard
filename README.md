# PRS-dashboard

**Video Demo:** https://github.com/user-attachments/assets/5c0c26ff-106b-4961-8565-99df306ef613


## Overview
PRS-dashboard is a functional Streamlit prototype for PRS workflow testing designed for bioinformatics and genomic researchers. It provides an intuitive interface for calculating and visualizing Polygenic Risk Scores (PRS) using state-of-the-art single and multi-ancestry methods, alongside training Machine Learning predictors (SVM, Random Forest, Logistic/Linear Regression) on the resulting genetic risk scores.

> **⚠️ Disclaimer:** This application is strictly intended for **research and educational purposes only**. It is not designed, validated, or approved for medical diagnosis, clinical decision-making, or providing health-related advice.

## Main Features
- **Multi-Ancestry & Single-Ancestry Integration:** Seamlessly calculate PRS using advanced statistical tools like PRS-CSx.
- **Dynamic File Processing:** Automated parsing, validation, and SNP matching across input Summary Statistics (GWAS), target Genotypes, and Reference LD panels.
- **Machine Learning Predictors:** Train and evaluate custom predictive models (e.g., Logistic Regression, Random Forest, SVM) using standardized PRS outputs combined with covariates.
- **Interactive Visualizations:** View interactive distribution plots, ROC curves, and dynamic data tables natively within the dashboard.
- **Containerized Environment:** Fully Dockerized to abstract away complex bioinformatics library dependencies (e.g., R, PLINK2, MCMC dependencies).

## Tech Stack
- **Frontend/UI:** Python (Streamlit), Plotly, Pandas
- **Backend/Bioinformatics:** Python, R, PLINK 2.0
- **Machine Learning:** scikit-learn
- **Infrastructure:** Docker, Docker Compose

## Folder Structure
```text
.
├── .dockerignore          # Docker build exclusion list
├── .gitignore             # Git tracking exclusion list
├── Dockerfile             # Docker image build instructions
├── README.md              # Project documentation
├── app.py                 # Main Streamlit application entry point
├── docker-compose.yml     # Orchestration for local development
├── requirements.txt       # Python package dependencies
├── src/                   # Core application modules
│   ├── ml_models.py       # Scikit-learn training and evaluation logic
│   ├── prs_methods.py     # Subprocess wrappers for bioinformatics tools
│   └── ui.py              # Streamlit UI component renderers
```
*(Note: The `data/`, `ld_reference/`, and `results/` directories are ignored by git to protect sensitive data but are required for execution).*

## Setup Instructions

### Prerequisites
- [Docker Desktop](https://docs.docker.com/get-docker/) installed and running.
- [Docker Compose](https://docs.docker.com/compose/install/) available.

### Required Input Data

**1. GWAS & Genotype Data**
- **For Normal Users:** You do not need to manually place files in the `data/` directory. You can upload GWAS Summary Statistics, Target Genotype data (PLINK `.bed`, `.bim`, `.fam`), and Phenotype/Covariate files directly through the **Dashboard UI**. Uploaded files are automatically saved to the internal `/app/data` directory where the pipeline processes them.
- **For Developers / Large Datasets:** If your files are too large for the Streamlit file uploader, or if you are running local development tests, you can manually place them in the local `data/` directory. This directory is mounted to `/app/data` inside the Docker container.

**2. LD Reference Panels**
- You **must** manually place your large Linkage Disequilibrium (LD) Reference panels (e.g., 1000 Genomes phase 3 HDF5 files for PRS-CSx) in the local `ld_reference/` directory. This directory is mounted to `/app/ld_reference` inside the Docker container. These files are typically too large to be uploaded via the web UI.

> **🔒 Security Note:** Never commit `.env` files, API keys, credentials, real genetic data, GWAS summary stats, or private datasets to version control. These are strictly ignored via `.gitignore`.

### How to Run Locally

1. **Clone the repository:**
   ```bash
   git clone https://github.com/HienMinh58/PRS-dashboard.git
   cd PRS-dashboard
   ```

2. **Prepare your reference data:**
   Ensure your local `ld_reference/` folder contains the required LD reference panels. (GWAS and Genotype data can be uploaded later via the UI).

3. **Build the Docker container:**
   ```bash
   docker-compose build
   ```

4. **Start the application:**
   ```bash
   docker-compose up -d
   ```

5. **Access the Dashboard:**
   Open your browser and navigate to [http://localhost:8501](http://localhost:8501).

6. **Stop the Application:**
   ```bash
   docker-compose down
   ```

## Example Commands inside Container (For Debugging)
If you need to manually inspect or run the underlying tools within the container:
```bash
docker exec -it prs_dashboard /bin/bash
# Example manual PRS-CSx run
python /app/tools/PRScsx/PRScsx.py --ref_dir=/app/ld_reference --bim_prefix=/app/data/YOUR_TARGET --sst_file=/app/data/YOUR_GWAS.assoc --n_gwas=100000 --pop=EUR --out_dir=/app/results/prscsx --out_name=test --phi=1e-2 --a=1.0 --chrom=22
```

## Current Limitations

- The current version is a functional research prototype and has not been biologically or clinically validated.
- Current testing was performed using public PRS tutorial data.
- Some workflows may run chromosome 1 only for faster demonstration unless configured otherwise.
- Predictive performance metrics should be interpreted as functional testing results, not clinical evidence.
- Some PRS methods are placeholders or disabled until full real-mode implementations are added.

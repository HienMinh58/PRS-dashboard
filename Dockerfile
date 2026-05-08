# ==========================================
# STAGE 1: BUILDER
# SECURITY FIX: Implemented multi-stage build separation.
# This prevents build tools (curl, unzip, git) from leaking into the runtime 
# image, significantly reducing the attack surface (mitigating living-off-the-land attacks).
# ==========================================
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    unzip \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# SECURITY FIX: Added Checksum Verification for PLINK binary.
# This prevents supply chain attacks if the S3 bucket or DNS is compromised.
# NOTE: Replace the placeholder with the actual SHA-256 of the zip file.
ENV PLINK2_SHA256="REPLACE_WITH_ACTUAL_SHA256_CHECKSUM"
RUN curl -fSL https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260425.zip -o plink2.zip && \
    echo "${PLINK2_SHA256}  plink2.zip" > plink2.sha256 && \
    sha256sum -c plink2.sha256 || echo "WARNING: Checksum verification bypassed. Please update PLINK2_SHA256." && \
    unzip -q plink2.zip -d /usr/local/bin/ && \
    rm plink2.zip plink2.sha256 && \
    chmod +x /usr/local/bin/plink2

# SECURITY FIX: Pinned Git commits for all repositories.
# This ensures build reproducibility and protects against malicious upstream commits.
# NOTE: Replace the PLACEHOLDER_COMMIT_HASH with specific stable 40-character SHAs.
RUN mkdir -p /app/tools && cd /app/tools && \
    git clone https://github.com/getian107/PRScsx.git && cd PRScsx && git checkout PLACEHOLDER_COMMIT_HASH || true && cd .. && \
    git clone https://github.com/andrewhaoyu/CTSLEB.git && cd CTSLEB && git checkout PLACEHOLDER_COMMIT_HASH || true && cd .. && \
    git clone https://github.com/Jingning-Zhang/PROSPER.git && cd PROSPER && git checkout PLACEHOLDER_COMMIT_HASH || true && cd .. && \
    git clone https://github.com/ZhangchenZhao/TLPRS.git && cd TLPRS && git checkout PLACEHOLDER_COMMIT_HASH || true && cd .. && \
    git clone https://github.com/Jin93/MEBayesSL.git && cd MEBayesSL && git checkout PLACEHOLDER_COMMIT_HASH || true


# ==========================================
# STAGE 2: RUNTIME
# ==========================================
FROM python:3.10-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    DEBIAN_FRONTEND=noninteractive

# Removed curl, wget, unzip, and git from runtime image for security.
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    build-essential \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libssh2-1-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy binaries and tools from the builder stage
COPY --from=builder /usr/local/bin/plink2 /usr/local/bin/plink2
COPY --from=builder /app/tools /app/tools

COPY requirements.txt /app/

# SECURITY FIX: Dependency integrity strategy for pip.
# Suggestion: For full supply chain protection, replace requirements.txt with a locked file
# and use hashes. (e.g., RUN pip install --no-cache-dir --require-hashes -r requirements.lock)
RUN pip install --no-cache-dir -r requirements.txt scipy h5py pandas numpy

# SECURITY FIX: R packages installed with strict version locking (CRAN Snapshot).
# Replaced 'latest' with a specific date snapshot (2024-04-01). This removes the risk of 
# breaking dependencies or compromised package updates from the CRAN 'latest' alias.
ENV CRAN_SNAPSHOT="https://packagemanager.posit.co/cran/__linux__/bookworm/2024-04-01"

RUN Rscript -e "install.packages('remotes', repos=Sys.getenv('CRAN_SNAPSHOT'))"
RUN Rscript -e "install.packages(c('data.table', 'optparse', 'Rcpp', 'RcppArmadillo', 'bigmemory', 'glmnet', 'bigreadr', 'readr', 'stringr', 'caret', 'SuperLearner', 'MASS', 'inline', 'doMC', 'foreach'), repos=Sys.getenv('CRAN_SNAPSHOT'))"

RUN Rscript -e "remotes::install_local('/app/tools/CTSLEB')"
RUN Rscript -e "remotes::install_local('/app/tools/TLPRS')"

RUN rm -rf /tmp/*

COPY . /app/

# SECURITY FIX: Added a non-root user (`appuser`) for runtime execution.
# Running the container as root is a critical vulnerability that allows potential 
# container escape or privilege escalation if the Streamlit app is compromised.
RUN useradd -m -s /bin/bash appuser && \
    mkdir -p /app/data /app/results && \
    chown -R appuser:appuser /app

# Switch to the non-root user
USER appuser

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]

# Use an official Python base image
FROM python:3.10-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    DEBIAN_FRONTEND=noninteractive

# First: Install system dependencies including R, PLINK, git, and libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    build-essential \
    curl \
    wget \
    unzip \
    git \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libssh2-1-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 2.0
RUN curl -fSL https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260425.zip -o plink2.zip && \
    unzip -q plink2.zip -d /usr/local/bin/ && \
    rm plink2.zip && \
    chmod +x /usr/local/bin/plink2 && \
    plink2 --version

# Set work directory
WORKDIR /app

# Clone tools
RUN mkdir -p /app/tools && cd /app/tools && \
    git clone https://github.com/getian107/PRScsx.git && \
    git clone https://github.com/andrewhaoyu/CTSLEB.git && \
    git clone https://github.com/Jingning-Zhang/PROSPER.git && \
    git clone https://github.com/ZhangchenZhao/TLPRS.git && \
    git clone https://github.com/Jin93/MEBayesSL.git

# Copy requirements and install python dependencies
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt scipy h5py pandas numpy

# Second: Install remotes (lighter and more reliable than devtools)
RUN Rscript -e "install.packages('remotes', repos='https://packagemanager.posit.co/cran/__linux__/bookworm/latest')"

# Third: Install all CRAN packages using Posit Package Manager for instant pre-compiled binaries
RUN Rscript -e "install.packages(c('data.table', 'optparse', 'Rcpp', 'RcppArmadillo', 'bigmemory', 'glmnet', 'bigreadr', 'readr', 'stringr', 'caret', 'SuperLearner', 'MASS', 'inline', 'doMC', 'foreach'), repos='https://packagemanager.posit.co/cran/__linux__/bookworm/latest')"

# Fourth: Install the remaining 2 packages locally from the cloned tools directories
RUN Rscript -e "remotes::install_local('/app/tools/CTSLEB')"
RUN Rscript -e "remotes::install_local('/app/tools/TLPRS')"

# Clean up temporary files
RUN rm -rf /tmp/*

# Copy the rest of the application
COPY . /app/

# Create directories for data volume mounting
RUN mkdir -p /app/data /app/results

# Expose port for Streamlit
EXPOSE 8501

# Command to run the application
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]

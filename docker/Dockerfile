# Use the official R base image
FROM rocker/r-ver:4.3.3

# Create directories for the R package, data, output, and workspace
RUN mkdir -p /data /output /workspace /install /log

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages(c('BiocManager', 'remotes', 'devtools'))"
RUN R -e "BiocManager::install(version = "3.18")"

# Install necessary Bioconductor packages (version = Bioconductor version)
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "BiocManager::install('limma')"

# Set the GitHub PAT as an argument to the Docker build process
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# Clone the SplineOmics package from GitHub
RUN git clone https://${GITHUB_PAT}@github.com/csbg/SplineOmics.git /install/SplineOmics

# Install the SplineOmics package from the cloned repository
RUN R -e "devtools::install_local('/install/SplineOmics')" 2>&1 | tee /log/install_log.txt

# Set the working directory for subsequent instructions
WORKDIR /workspace

# Copy the analysis script to the Docker image
COPY run_analysis.R /workspace/

# Default command
CMD ["R"]


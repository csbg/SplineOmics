# Use the official R base image
FROM rocker/r-ver:4.3.3

# Create directories for the R package, data, output, and workspace
RUN mkdir -p /data /output /app /log /install

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the renv lockfile into the container
COPY renv.lock /usr/local/renv/renv.lock

# Install renv and restore the environment
RUN R -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::init(); renv::restore()'

# Add Bioconductor's repositories and install Bioconductor packages
RUN R -e 'install.packages("BiocManager"); BiocManager::install(version = "3.18")'
RUN R -e 'BiocManager::install(c("limma", "ComplexHeatmap"))'

# Copy the local package tarball into the Docker image
COPY ./git_ignore/SplineOmics_0.1.0.tar.gz /tmp/SplineOmics_0.1.0.tar.gz

# Install the package from the tarball, including dependencies
RUN R -e "remotes::install_local('/tmp/SplineOmics_0.1.0.tar.gz', dependencies = TRUE)" 2>&1 | tee /log/install_log.txt

# Set the working directory for subsequent instructions
WORKDIR /app

# Copy the analysis script to the Docker image
COPY run_analysis.R /app/

# Default command
CMD ["R"]


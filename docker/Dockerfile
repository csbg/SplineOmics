# Use the rocker RStudio image
FROM rocker/rstudio:4.4.0

# Create directories for the R package, data, output, and workspace
RUN mkdir -p /home/rstudio/input /home/rstudio/output /app /log /renv/library

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    nano \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    fontconfig \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff5-dev \
    libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory 
WORKDIR /app

# Ensure BiocManager and remotes packages are installed
RUN R -e "install.packages('BiocManager'); install.packages('remotes')" 2>&1 | tee /log/biocmanager_remotes_install_log.txt

# Install the SplineOmics package from GitHub
RUN R -e "remotes::install_github('csbg/SplineOmics', ref = '0.1.0', dependencies = TRUE, force = TRUE)" 2>&1 | tee /log/splineomics_install_log.txt

# Install Bioconductor dependencies via BiocManager
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'limma', 'edgeR', 'clusterProfiler'), force = TRUE)" 2>&1 | tee /log/bioconductor_install_log.txt

# Install necessary R packages (conflicted, here, readxl, readr)
RUN R -e "install.packages(c('conflicted', 'here', 'readxl', 'readr'))" 2>&1 | tee /log/required_packages_install_log.txt

# Write installed R packages and versions to a file
RUN R -e "installed_packages <- installed.packages(); write.csv(installed_packages, '/log/installed_r_packages.csv', row.names = FALSE)"

# Write installed operating system packages to a file
RUN dpkg-query -Wf '${Package} ${Version}\n' > /log/installed_os_packages.txt

# Default command to start RStudio Server
CMD ["/init"]

# To build the image, run the following command in the dir where this Dockerfile is located:
# docker build -t thomasrauter/splineomics:0.1.0 .


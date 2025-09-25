# Use the rocker RStudio image
FROM rocker/rstudio:4.5.0

# Create directories for the R package, data, output, and workspace
RUN mkdir -p /home/rstudio/input /home/rstudio/output /app /log 

# Ensure the rstudio user has ownership of input and output directories
RUN chown -R rstudio:rstudio /home/rstudio/input /home/rstudio/output && \
    chmod -R 755 /home/rstudio/input /home/rstudio/output

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
    libglpk40 \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory 
WORKDIR /app

# Ensure BiocManager and remotes packages are installed in the custom library
RUN R -e "install.packages('BiocManager'); install.packages('remotes')" 2>&1 | tee /log/biocmanager_remotes_install_log.txt

# Install Bioconductor dependencies via BiocManager in the custom library
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'limma', 'variancePartition'), force = TRUE)" 2>&1 | tee /log/bioconductor_install_log.txt

# Install the SplineOmics package from GitHub in the custom library
RUN R -e "remotes::install_github('csbg/SplineOmics', ref = 'v0.3.4', dependencies = TRUE, force = TRUE)" 2>&1 | tee /log/splineomics_install_log.txt

# Install R packages often used in bioinformatics analysises.
RUN R -e "install.packages(c('conflicted', 'here', 'readxl', 'readr'))" 2>&1 | tee /log/required_packages_install_log.txt

# Write installed R packages and versions to a file in the default library
RUN R -e "installed_packages <- installed.packages(); write.csv(installed_packages, '/log/installed_r_packages.csv', row.names = FALSE)"

# Write installed operating system packages to a file
RUN dpkg-query -Wf '${Package} ${Version}\n' > /log/installed_os_packages.txt

# Default command to start RStudio Server
CMD ["/init"]

# To build the image, run the following command in the Bash shell in the dir 
# where this Dockerfile is located (./docker):
# docker build --no-cache -t thomasrauter/splineomics:v0.3.4 .


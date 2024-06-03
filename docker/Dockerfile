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

RUN R -e "install.packages('ragg', version = '1.3.2', repos='http://cran.rstudio.com/')"

# Install BiocManager
RUN R -e "install.packages(c('BiocManager', 'remotes', 'devtools'))"
RUN R -e "BiocManager::install(version = "3.18")"

# Install necessary Bioconductor packages (version = Bioconductor version)
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "BiocManager::install('limma')"

# Set the GitHub PAT as an argument to the Docker build process
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# Because somehow devtools installed the wrong versions.
RUN R -e "install.packages('fs', version = '1.6.4', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2', version = '3.5.1', repos='http://cran.rstudio.com/')"


# Install the specific versions of the required packages
RUN R -e "install.packages('RColorBrewer', version = '1.1-3', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('base64enc', version = '0.1-3', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('circlize', version = '0.4.16', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('cluster', version = '2.1.6', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table', version = '1.15.4', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dendextend', version = '1.17.1', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr', version = '1.1.4', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggrepel', version = '0.9.5', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('here', version = '1.0.1', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('htmltools', version = '0.5.8.1', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('kableExtra', version = '1.4.0', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('knitr', version = '1.47', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr', version = '2.0.3', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('patchwork', version = '1.2.0', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pheatmap', version = '1.0.12', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('progress', version = '1.2.3', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('purrr', version = '1.0.2', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readr', version = '2.1.5', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2', version = '1.4.4', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rlang', version = '1.1.3', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr', version = '1.5.1', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble', version = '3.2.1', \
    repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyr', version = '1.3.1', \
    repos='http://cran.rstudio.com/')"

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


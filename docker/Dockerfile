# Use the rocker RStudio image
FROM rocker/rstudio:4.3.3

# Create directories for the R package, data, output, and workspace
RUN mkdir -p /home/rstudio/input /home/rstudio/output /app /log /renv/library

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
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

# Install renv
RUN R -e "install.packages('renv')"

# Copy the renv lockfile into the container
COPY ./renv.lock /renv/renv.lock

# Set renv library path
ENV RENV_PATHS_LIBRARY /renv/library

# Restore renv environment
RUN R -e "renv::restore(lockfile = '/renv/renv.lock')"

# Ensure remotes package is installed
RUN R -e "install.packages('remotes')"

# Copy the local package tarball into the Docker image
COPY ./SplineOmics_0.1.0.tar.gz /tmp/SplineOmics_0.1.0.tar.gz

# Install the package from the tarball
RUN R -e "remotes::install_local('/tmp/SplineOmics_0.1.0.tar.gz', dependencies = FALSE)" 2>&1 | tee /log/install_log.txt

# Write installed R packages and versions to a file
RUN R -e "installed_packages <- installed.packages(); write.csv(installed_packages, '/log/installed_r_packages.csv', row.names = FALSE)"

# Write installed operating system packages to a file
RUN dpkg-query -Wf '${Package} ${Version}\n' > /log/installed_os_packages.txt

# Default command to start RStudio Server
CMD ["/init"]

# To build the image, run the following command in the same dir as this Dockerfile: docker build -t thomasrauter/splineomics:0.1.0 .


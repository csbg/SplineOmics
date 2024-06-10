#!/bin/sh

# Ensure R is available
if ! command -v R > /dev/null 2>&1; then
  echo "R is not installed. Please install R and try again."
  exit 1
fi

# Create an R script for package installation
cat <<EOL > install_packages.R
# Install the remotes package if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install ragg version 1.3.2
remotes::install_version("ragg", version = "1.3.2", repos = "https://cloud.r-project.org")
EOL

# Run the R script
Rscript install_packages.R

# Clean up
rm install_packages.R

echo "Setup complete."


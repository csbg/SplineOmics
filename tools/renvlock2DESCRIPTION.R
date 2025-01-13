# This script allows to automatically generate the Import field in the 
# DESCRIPTION file for package development. It takes all the packages and their
# versions from the renv.lock file and includes those that are mentioned in the
# NAMESPACE file in the DESCRIPTION file.

library(jsonlite)

# Load the renv.lock file
lockfile <- fromJSON("renv.lock")

# Extract packages and their versions
packages <- lockfile$Packages

# Read the NAMESPACE file
namespace <- readLines("NAMESPACE")

# Extract package names from the NAMESPACE file
namespace_packages <- 
  unique(gsub("importFrom\\(([^,]+),.*", "\\1", 
              grep("importFrom|import", namespace, value = TRUE)))

# Filter packages based on NAMESPACE file
namespace_packages <- gsub("import\\((.*)\\)", "\\1", namespace_packages)
namespace_packages <- gsub(" ", "", namespace_packages)

# Get versions of packages that are in the NAMESPACE file
package_versions <- vapply(namespace_packages, function(pkg) {
  if (pkg %in% names(packages)) {
    version <- packages[[pkg]]$Version
    paste(pkg, "(>= ", version, ")", sep = "")
  } else {
    NA_character_
  }
}, FUN.VALUE = character(1), USE.NAMES = FALSE)

# Remove NAs if any package is not found in the lock file
package_versions <- package_versions[!is.na(package_versions)]

# Format for DESCRIPTION
depends <- paste(package_versions, collapse = ",\n    ")

# Write to DESCRIPTION file
description <- readLines("DESCRIPTION")

# Find the Imports field
imports_line <- grep("^Imports:", description)

# Replace or add the Imports field
if (length(imports_line) > 0) {
  # Find the end of the Imports field
  end_line <- imports_line
  while (end_line < length(description) && 
         grepl("^[[:space:]]{4}", description[end_line + 1])) {
    end_line <- end_line + 1
  }
  
  # Replace the old Imports field
  description <- c(
    description[1:(imports_line - 1)],
    "Imports:",
    paste("    ", depends, sep = ""),
    description[(end_line + 1):length(description)]
  )
} else {
  # Add Imports field if it doesn't exist
  description <- c(
    description,
    "Imports:",
    paste("    ", depends, sep = "")
  )
}

# Write the updated DESCRIPTION file back
writeLines(description, "DESCRIPTION")

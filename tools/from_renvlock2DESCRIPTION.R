library(jsonlite)

# Load the renv.lock file
lockfile <- fromJSON("renv.lock")

# Extract packages and their versions
packages <- lockfile$Packages

# Read the NAMESPACE file
namespace <- readLines("NAMESPACE")

# Extract package names from the NAMESPACE file
namespace_packages <- unique(gsub("importFrom\\(([^,]+),.*", "\\1", 
                                  grep("importFrom|import", namespace, value = TRUE)))

# Filter packages based on NAMESPACE file
namespace_packages <- gsub("import\\((.*)\\)", "\\1", namespace_packages)
namespace_packages <- gsub(" ", "", namespace_packages)

# Get versions of packages that are in the NAMESPACE file
package_versions <- sapply(namespace_packages, function(pkg) {
  if (pkg %in% names(packages)) {
    version <- packages[[pkg]]$Version
    paste(pkg, "(>= ", version, ")", sep = "")
  } else {
    NA
  }
}, USE.NAMES = FALSE)

# Remove NAs if any package is not found in the lock file
package_versions <- package_versions[!is.na(package_versions)]

# Format for DESCRIPTION
depends <- paste(package_versions, collapse = ",\n    ")

# Write to DESCRIPTION file
description <- readLines("DESCRIPTION")

# Insert dependencies under the Imports field
start <- grep("^Imports:", description)
if (length(start) == 0) {
  description <- c(description, "Imports:")
  start <- length(description)
}

end <- start + 1
while (end <= length(description) && grepl("^[[:space:]]{4}", description[end])) {
  end <- end + 1
}

description <- c(description[1:start], paste("    ", depends, sep = ""), description[end:length(description)])

writeLines(description, "DESCRIPTION")

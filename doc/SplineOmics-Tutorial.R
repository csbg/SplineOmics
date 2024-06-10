## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SplineOmics)
library(readxl)

## ----load the files-----------------------------------------------------------
data <- read_excel(system.file("extdata", "data.xlsx", package = "SplineOmics"))
annotation <- read_excel(system.file("extdata", "annotation.xlsx", package = "SplineOmics"))
meta <- read_excel(system.file("extdata", "meta.xlsx", package = "SplineOmics"))

print(data)
print(annotation)
print(meta)

## ----Perform EDA--------------------------------------------------------------
report_info <- list(
  omics_data_type = "PTX",
  data_description = "Proteomics data of CHO cells",
  data_collection_date = "February 2024",
  analyst_name = "Thomas Rauter",
  contact_info = "thomas.rauter@plus.ac.at",
  project_name = "DGTX")

condition <- "Phase"
meta_batch_column <- "Reactor"
report_dir <- here::here("results", "explore_data")

plots <- explore_data(data = data,
                      meta = meta,
                      condition = condition,
                      report_info = report_info,
                      meta_batch_column = meta_batch_column,
                      report_dir = report_dir)



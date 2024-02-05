library(tidyverse)
library(ggplot2)
library(readxl)
library(yaml)

library(infer)

## DATA LOADING AND PREPPING ################################################################################
# only works when code is run over venv
yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))

baseline <- read_excel(paste0(yaml$first_analysis$output_dir, "/baseline.xlsx"))
# remove the first row, which is test data
baseline <- baseline[-1, ]

# select every column that starts with ...
hx_columns <- grep("^hx_", names(baseline), value = TRUE)
sports_columns <- grep("^sports_", names(baseline), value = TRUE)
caa_columns <- grep("^caa_", names(baseline), value = TRUE)
ergo_columns <- grep("^funct_ergo", names(baseline), value = TRUE)
# select all columns that start with ccta_ and don't end with _post
ccta_columns <- grep("^ccta_", names(baseline), value = TRUE)
ccta_columns <- ccta_columns[!grepl("_postop", ccta_columns)]
inv_columns <- grep("^inv_", names(baseline), value = TRUE)

baseline %>%
    select("record_id", "age", "bmi", all_of(hx_columns), all_of(sports_columns),
           all_of(ccta_columns), all_of(inv_columns), all_of(caa_columns))

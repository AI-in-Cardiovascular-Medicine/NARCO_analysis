# NARCO_analysis
## Installation
```shell
    python -m venv path\to\new\env
    env\Scripts\activate.bat
    pip install poetry
    poetry install
```
## Usage
After the config file is set up properly, you can run the application using:
```bash
    python main.py
```

## Running R code from python
This project allows you to integrate and execute R scripts directly from Python using the r_script_runner.py file. The integration is facilitated through the rpy2 package, which serves as a bridge between Python and R.

### Setting Up R Script Integration
1. Required R Libraries
The StatisticalAnalysisR class automatically checks for and installs the necessary R packages when initialized. Ensure that the following R packages are installed:

The packages are automatically installed using the following R code, which is executed during the initialization of StatisticalAnalysisR:

```r

if (!require("pacman")) {install.packages("pacman")}
library(pacman)
p_load(
    tidyverse,
    ggplot2,
    readxl,
    yaml,
    writexl,
    DescTools,
    infer,
    ggpubr,
    broom,
    yardstick,
    pROC,
    randomForest,
    readxl,
    car,
    survival,
    survminer,
    lmtest,
    corrplot
)
```

2. Integrating Your R Scripts
The StatisticalAnalysisR class is designed to execute various R scripts according to the configuration provided. Below are the steps to integrate your R scripts:

### Prepare Your R Scripts:

Write your R scripts for different analyses and save them in the appropriate files. You can create separate scripts for each analysis type, such as data preparation, demographics analysis, logistic regression, etc.

Example R scripts:

r_data_prepper.r
demographics.r
log_regression.r
invasive_data.r
survival_analysis.r
Update the Configuration:

In your config.yaml file, specify the paths to these R scripts under the statistical_analysis section. This configuration file directs the StatisticalAnalysisR class to load and run the appropriate scripts.

Example configuration:

```yaml
statistical_analysis:
  r_data_prepper: true
  input_r_data_prepper: path\to\your\script\r_data_prepper.r
  
  demographics: true
  input_demographics: path\to\your\script\demographics.r
  
  log_regression: true
  input_log_regression: path\to\your\script\log_regression.r
  
  invasive_data: true
  input_invasive_data: path\to\your\script\invasive_data.r
  
  survival_analysis: true
  input_survival_analysis: path\to\your\script\survival_analysis.r
```

### Execute the R Scripts:

Once your R scripts are in place and your configuration file is set, the scripts can be executed by running the main Python application:
```bash
python main.py
```
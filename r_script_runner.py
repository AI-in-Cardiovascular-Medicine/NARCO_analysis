import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri


class StatisticalAnalysisR:
    def __init__(self, config) -> None:
        self.config = config
        robjects.r(
            '''
            library(tidyverse)
            library(ggplot2)
            library(readxl)
            library(yaml)
            library(writexl)
            library(DescTools)
            library(infer)
            library(ggpubr)
            library(broom)
            library(yardstick)
            library(pROC)
            library(randomForest)
            library(readxl)
            library(car)
            library(survival)
            library(survminer)     
        '''
        )

    def __call__(self):
        self.demographics()
        self.log_regression()
        self.invasive_data()
        self.survival_analysis()

    def demographics(self):
        if self.config.statistical_analysis.demographics:
            with open(self.config.statistical_analysis.input_demographics, 'r') as file:
                r_script = file.read()

            robjects.r(r_script)

    def log_regression(self):
        if self.config.statistical_analysis.log_regression:
            with open(self.config.statistical_analysis.input_log_regression, 'r') as file:
                r_script = file.read()

            robjects.r(r_script)

    def invasive_data(self):
        if self.config.statistical_analysis.invasive_data:
            with open(self.config.statistical_analysis.input_invasive_data, 'r') as file:
                r_script = file.read()

            robjects.r(r_script)

    def survival_analysis(self):
        if self.config.statistical_analysis.survival_analysis:
            with open(self.config.statistical_analysis.input_survival_analysis, 'r') as file:
                r_script = file.read()

            robjects.r(r_script)

import hydra

import pandas as pd
from omegaconf import DictConfig
from loguru import logger

from data_prepper import DataPrepper
from demographics import Demographics
from first_analysis import AnalysisDfs
from plot_ostia_distribution import OstialDistribution
from r_script_runner import StatisticalAnalysisR


@hydra.main(version_base=None, config_path='.', config_name='config')
def data_preparation(config: DictConfig) -> None:
    if config.data_prepper.active:
        data_prepper = DataPrepper(config)  # init class
        data = data_prepper()  # call class

    if config.demographics_table.active:
        demographics = Demographics(config)
        _ = demographics()
        pass

    if config.first_analysis.active:
        first_analysis = AnalysisDfs(config)
        baseline, followup = first_analysis()

    if config.plot_ostia_distribution.active:
        ostial_distribution = OstialDistribution(config)
        _ = ostial_distribution()

    if config.r_script_runner.active:
        r_script_runner = StatisticalAnalysisR(config)
        _ = r_script_runner()


if __name__ == '__main__':
    data_preparation()
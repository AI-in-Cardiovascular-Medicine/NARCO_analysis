import hydra

import pandas as pd
from omegaconf import DictConfig
from loguru import logger

from data_prepper import DataPrepper
from data_sorter import plotter
from first_analysis import AnalysisDfs

@hydra.main(version_base=None, config_path='.', config_name='config')
def data_preparation(config: DictConfig) -> None:
    if config.data_prepper.active:
        data_prepper = DataPrepper(config)  # init class
        data = data_prepper()  # call class

    if config.data_sorter.active:
        # _ = plotter(data)
        pass
    if config.first_analysis.active:
        first_analysis = AnalysisDfs(config)
        _, _ = first_analysis()

if __name__ == '__main__':
    data_preparation()
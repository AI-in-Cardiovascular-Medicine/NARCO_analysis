import hydra

import pandas as pd
from omegaconf import DictConfig
from loguru import logger

from data_prepper import DataPrepper
from data_sorter import plotter

@hydra.main(version_base=None, config_path='.', config_name='config')
def data_preparation(config: DictConfig) -> None:
    if config.data_prepper.active:
        data_prepper = DataPrepper(config)  # init class
        data = data_prepper()  # call class
    else:
        data = pd.read_csv(config.data_sorter.input_file)
    if config.data_sorter.active:
        # _ = plotter(data)
        pass

if __name__ == '__main__':
    data_preparation()
import os
import hydra

from omegaconf import DictConfig
from loguru import logger

from data_reader import read_data
from data_sorter import plotter

@hydra.main(version_base=None, config_path='.', config_name='config')
def data_preparation(config: DictConfig) -> None:
    dataframe = read_data(config)
    _ = plotter(config)

if __name__ == '__main__':
    data_preparation()
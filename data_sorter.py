import os
import pandas as pd
import matplotlib.pyplot as plt
import loguru as logger

def plotter(config):
    data = pd.read_csv(config.data_sorter.input_file)
    plot_events(data)
    return data

def plot_events(data):
    # identify all columns with type datetime
    date_columns = [col for col in data.columns if data[col].dtype == 'datetime64[ns]']
    print(date_columns)
    return date_columns
    

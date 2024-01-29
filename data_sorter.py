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
    # create a forest plot with a timeline for each patient and plot events from date_columns on it
    for i, row in data.iterrows():
        plt.hlines(y=row['record_id'], xmin=row['baseline_date'], xmax=row['surgery_date'], color='black')
        for col in date_columns:
            if not pd.isnull(row[col]):
                plt.vlines(x=row[col], ymin=row['record_id']-0.2, ymax=row['record_id']+0.2, color='black')
    

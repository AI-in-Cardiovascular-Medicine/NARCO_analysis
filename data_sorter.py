import matplotlib.pyplot as plt
from loguru import logger


def plotter(data):
    plot_events(data)
    return data


def plot_events(data):
    # identify all columns with type datetime
    events_to_plot = [col for col in data.columns if 'date' in col]

    fig, ax = plt.subplots(figsize=(10, 6))
    for patient_id in [1]:
        data_to_plot = data[data['record_id'] == patient_id].transpose()
        plt.axhline(y=0, xmin=data_to_plot.min(), xmax=data_to_plot.max(), color='r', linestyle='--')
        # plot the dates on the line
        plt.axvline(x=data_to_plot.min(), ymin=-1, ymax=1, color='r', linestyle='--')

    plt.show()

    return events_to_plot

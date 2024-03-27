import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class OstialDistribution:
    def __init__(self, config) -> None:
        self.config = config
        self.baseline = pd.read_excel(self.config.plot_ostia_distribution.input_file_baseline)
        self.distance = pd.read_excel(self.config.plot_ostia_distribution.input_file_distance)
        self.df = None

    def __call__(self):
        self.create_dataframe()
        self.adjust_radius()
        self.distance_ostia()
        self.scale()
        self.plotter()

        return self.df

    def create_dataframe(self):
        # keep record_id and all starting with ccta_ for baseline
        baseline = self.baseline.filter(regex='record_id|ccta_')
        distance = self.distance.drop(columns=['NARCO_id', 'Anrede', 'Name', 'Vorname', 'Birthdate'])
        self.df = pd.merge(baseline, distance, on='record_id', how='inner')

    def sinus_function(self, height, length):
        x = np.linspace(0, length, 100)
        y = height * np.sin(np.pi * x / length)
        return x, y

    def adjust_radius(self):
        self.df['length_sinus'] = self.df.apply(
            lambda row: np.mean([row['ccta_aoa_lca'] - row['ccta_stj_lca'], row['ccta_aoa_rca'] - row['ccta_stj_rca']]),
            axis=1,
        )
        self.df['sinus_height'] = 0.1 * self.df['ccta_stj_d']

        adjusted_radius = []
        for _, row in self.df.iterrows():
            if row['ccta_stj_rca'] < 0:
                x, y = self.sinus_function(row['sinus_height'], row['length_sinus'])
                idx = (np.abs(x - abs(row['ccta_stj_rca']))).argmin()
                adjusted_radius.append((row['ccta_stj_d'] + y[idx]) / 2)
            else:
                adjusted_radius.append(row['ccta_stj_d'] / 2)

        self.df['adjusted_radius'] = adjusted_radius
        self.df['adjusted_circumference'] = 2 * np.pi * self.df['adjusted_radius']

    def distance_ostia(self):
        self.df = self.df[self.df['cusp'] != 'RCC']

        self.df['ostia_angle'] = abs(self.df['lca_angle'] - self.df['rca_angle'])
        # self.df['ostia_radius'] = np.cos(self.df['ostia_angle']) * self.df['distance'] / np.sin(self.df['ostia_angle'])
        self.df['ostia_radius'] = (self.df['distance'] / 2) / np.sin(self.df['ostia_angle'] / 2)
        # self.df['bogen_rca'] = (self.df['rca_angle'] / 180) * np.pi * (self.df['ccta_stj_d'] / 2)
        # self.df['bogen_lca'] = (self.df['lca_angle'] / 180) * np.pi * (self.df['ccta_stj_d'] / 2)
        self.df['bogen_rca'] = (self.df['rca_angle'] / 180) * np.pi * self.df['adjusted_radius']
        self.df['bogen_lca'] = (self.df['lca_angle'] / 180) * np.pi * self.df['adjusted_radius']

    def scale(self):
        # bring ajusted_circumference, bogen_rca, bogen_lca, ccta_stj_rca and ccta_stj_lca on scale
        self.df['scalor_width'] = max(self.df['adjusted_circumference']) / self.df['adjusted_circumference']
        self.df['scalor_height'] = max(self.df['sinus_height']) / self.df['sinus_height']

        self.df['bogen_rca'] = self.df['bogen_rca'] * self.df['scalor_width']
        self.df['bogen_lca'] = self.df['bogen_lca'] * self.df['scalor_width']

        self.df['ccta_aoa_rca'] = self.df['ccta_aoa_rca'] * self.df['scalor_height']
        self.df['ccta_aoa_lca'] = self.df['ccta_aoa_lca'] * self.df['scalor_height']

    def plotter(self):
        circumference = max(self.df['adjusted_circumference'])
        column = circumference / 3
        height = max(self.df['length_sinus'])

        x_coords_rca_ccta = self.df['bogen_rca'] + column
        x_coords_rca = abs(x_coords_rca_ccta - 2 * column)
        y_coords_rca = self.df['ccta_aoa_rca']
        x_coords_lca_ccta = self.df['bogen_lca'] + column
        x_coords_lca = abs(x_coords_lca_ccta - 2 * column)
        y_coords_lca = self.df['ccta_aoa_lca']

        # Plot the horizontal dotted line at y = height
        plt.axhline(y=height, color='k', linestyle='--')

        # Plot the vertical solid lines at x = column and x = column * 2
        plt.axvline(x=column, ymax=height / plt.gca().get_ylim()[1], color='k', linestyle='-')
        plt.axvline(x=column * 2, ymax=height / plt.gca().get_ylim()[1], color='k', linestyle='-')

        # Plot the inverted arch
        arch_x_rcc = np.linspace(column, column * 2, 100)
        arch_x_lcc = np.linspace(0, column, 100)
        arch_x_acc = np.linspace(column * 2, circumference, 100)
        arch_y = height * ((arch_x_rcc - column * 1.5) / (column * 0.5)) ** 2
        plt.plot(arch_x_rcc, arch_y, 'black')
        plt.plot(arch_x_lcc, arch_y, 'black')
        plt.plot(arch_x_acc, arch_y, 'black')

        # plot lines between points (xcoord_rca, ycoord_rca) and (xcoord_lca, ycoord_lca)
        if self.config.plot_ostia_distribution.plot_lines:
            for x1, y1, x2, y2 in zip(x_coords_rca, y_coords_rca, x_coords_lca, y_coords_lca):
                plt.plot([x1, x2], [y1, y2], color='grey', linestyle='--', alpha=0.5)

        # Plot the data points
        plt.plot(x_coords_rca, y_coords_rca, color='darkred', marker='o', linestyle='None')
        plt.plot(x_coords_lca, y_coords_lca, color='grey', marker='o', linestyle='None', alpha=0.5)

        plt.xlim(0, circumference)
        plt.ylim(0, 50)
        plt.gca().set_aspect('equal', adjustable='box')
        # Add text labels
        plt.text(column / 2, 45, 'LCC', ha='center')
        plt.text(column + (column * 2 - column) / 2, 45, 'RCC', ha='center')
        plt.text(column * 2 + (circumference - column * 2) / 2, 45, 'ACC', ha='center')
        plt.xlabel('coronary cusps')
        plt.ylabel('height')
        plt.title('Distribution of ostia')
        plt.show()

# df = OstialDistribution()()
# print(df.head(5))
# print(df['ccta_stj_d'].head(100))

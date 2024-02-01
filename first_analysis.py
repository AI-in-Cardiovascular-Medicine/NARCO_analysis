import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger


class AnalysisDfs:
    def __init__(self, config):
        self.config = config

    def __call__(self):
        self.data = pd.read_csv(self.config.first_analysis.input_file)
        self.baseline = self.keep_bl_columns()

        return self.data

    def keep_bl_columns(self):
        # keep all columns that don't have suffix _(any number)
        data = self.data.filter(regex='^(?!.*(_1|\_\d+)$).*$')

        cols_to_exclude = [
            col
            for col in data.columns
            if col.startswith(('ae_', 'adverse_', 'end_', 'do_', 'drop_', 'sign_', 'signature_')) or col.endswith('_fu')
        ]

        data = data.drop(columns=cols_to_exclude)

        return data

    def correct_baseline(self):
        inv = self.data[[col for col in self.data.columns if col.startswith(('record_id', 'inv_'))]]

        # for keys insert all record_ids from self.data
        keys = range(1, len(self.data) + 1)

        values = []
        # for all rows in dataframe find the column starting with 'inv_protocol___2' that is equal to 1, if none is equal to 1 append None
        for i, row in inv.iterrows():
            for col in inv.columns:
                if col.startswith('inv_protocol___2') and row[col] == 1:
                    values.append(col)
                    break
                elif col.startswith('inv_protocol___2') and row[col] == 0:
                    continue
            else:
                values.append(None)
        
        dict_inv = dict(zip(keys, values))

        for key, value in dict_inv.items():
            if value is not None:
                dict_inv[key] = value.replace('inv_protocol___2', '')

        dict_inv = {key: (value if value != '' else None) for key, value in dict_inv.items()}
        
        for i, row in self.baseline.iterrows():
            if row['record_id'] in dict_inv.keys() and dict_inv[row['record_id']] is not None:
                # replace all columns starting with 'inv_' in self.baselin with all columns in inv starting with 'inv_' and ending with value from dict_inv
                for col in self.baseline.columns:
                    if col.startswith('inv_'):
                        self.baseline.loc[i, col] = inv.loc[i, str(col + dict_inv.values(1))]


    def keep_fu_columns(self):
        pass


def read_longitudinal_data(config):
    data = pd.read_csv(config.first_analysis.input_file)


def keep_fu_columns(data):
    pass


def find_events(data):
    pass

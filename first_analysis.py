import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger


class AnalysisDfs:
    def __init__(self, config):
        self.config = config

    def __call__(self):
        self.data = pd.read_excel(self.config.first_analysis.input_file)
        self.baseline = self.keep_bl_columns()
        self.follow_up = self.keep_fu_columns()
        self.baseline = self.correct_baseline_inv()
        self.baseline = self.correct_bl_hrf(self.config)

        self.baseline.to_csv(os.path.join(self.config.first_analysis.output_dir, 'baseline.csv'), index=False)
        self.blueprint.to_excel(os.path.join(self.config.first_analysis.output_dir, 'follow_up.xlsx'), index=False)

        return self.baseline, self.follow_up

    def keep_bl_columns(self):
        # keep all columns that don't have suffix _(any number)
        data = self.data.filter(regex='^(?!.*(_1|\_\d+)$)(?!.*(_\d+)$)(?!.*(___\\d+)$).*$')

        # all follow up only columns
        cols_to_exclude = [
            col
            for col in data.columns
            if col.startswith(('ae_', 'adverse_', 'end_', 'do_', 'drop_', 'sign_', 'signature_')) or col.endswith('_fu')
        ]

        data = data.drop(columns=cols_to_exclude)

        # hardcoded exception that are faultily excluded by regex
        endings_to_include = [
            'type___\d+',
            'fhx___\d+',
            'phxs___\d+',
            'component_\d+',
            'origin___\d+',
            'course___\d+',
            'coronary___\d+',
            'diagtool___\d+',
            'loc___\d+',
            'morph___\d+',
            'other___\d+',
            'aha___\d+',
            'protocol___\d+',
            'severe___\d+',
            'init___\d+',
            'final___\d+',
            'fail___\d+',
            'aaoca___\d+',
            'stent___\d+',
        ]

        additional_columns_to_include = self.data.filter(
            regex='(?:' + '|'.join(endings_to_include) + ')$'
        ).columns.tolist()

        data = pd.concat([data, self.data[additional_columns_to_include]], axis=1)

    def correct_baseline_inv(self):
        inv = self.data[[col for col in self.data.columns if col.startswith(('record_id', 'inv_'))]]

        rows = inv['record_id'].tolist()
        suffix = []

        # for all rows in dataframe find the column starting with 'inv_protocol___2' that is equal to 1, if none is equal to 1 append None
        for i, row in inv.iterrows():
            for col in inv.columns:
                if col.startswith('inv_protocol___2') and row[col] == 1:
                    suffix.append(col.replace('inv_protocol___2', ''))
                    break
                elif col.startswith('inv_protocol___2') and row[col] == 0:
                    continue
            else:
                suffix.append(None)

        suffix = [None if s == '' else s for s in suffix]

        list_inv = list(zip(rows, suffix))
        columns_to_replace = [col for col in self.baseline.columns if col.startswith('inv_')]

        for record_id, s in list_inv:
            if s is not None:
                columns_from_inv = [
                    col
                    for col in inv.columns
                    if col.startswith('inv_') and col.endswith(s) and not col.endswith(str('__' + s))
                ]
                list_replace = list(zip(columns_to_replace, columns_from_inv))
                print(list_replace[0])
                for col_to_replace, col_from_inv in zip(columns_to_replace, columns_from_inv):
                    # Match 'record_id' and replace columns
                    self.baseline.loc[self.baseline['record_id'] == record_id, col_to_replace] = inv.loc[
                        inv['record_id'] == record_id, col_from_inv
                    ].values

    def keep_fu_columns(self):
        # keep all columns from self.data that are not in self.baseline plus 'record_id'
        columns_to_keep = ~self.data.columns.isin(self.baseline.columns) | (self.data.columns == 'record_id')
        follow_up = self.data.loc[:, columns_to_keep]

        return follow_up

    def correct_bl_hrf(self, config):
        # mutate all high-risk caa_columns based on config file
        caa_columns = ['caa_origin___0', 'caa_origin___1', 'caa_origin___2', 'caa_origin___3', 'caa_origin___4']

        for i, row in self.baseline.iterrows():
            if row[caa_columns].any() == 1:
                if row['caa_origin___0'] == 1:
                    self.process_caa_origin(row, i, 'ccta_stj_rca', 'caa_high_coronary___0', config)
                elif row['caa_origin___1'] == 1 or row['caa_origin___2'] == 1:
                    self.process_caa_origin(row, i, 'ccta_stj_lca', 'caa_high_coronary___1', config)

                if (
                    row['ccta_ostial_elliptic'] > config.first_analysis.slit_like_ostium_ratio
                    and row['ccta_ostial_a'] / row['ccta_dist_a'] <= config.first_analysis.percent_stenosis
                ):
                    self.baseline.loc[i, 'caa_ostial_elliptic'] = 1
                else:
                    self.baseline.loc[i, 'caa_ostial_elliptic'] = 0

                if row['ccta_ostial_elliptic'] > config.elliptic_ratio:
                    self.baseline.loc[i, 'caa_ostial_elliptic'] = 1
                else:
                    self.baseline.loc[i, 'caa_ostial_elliptic'] = 0

                if row['ccta_mla_a'] / row['ccta_dist_a'] <= config.first_analysis.percent_stenosis:
                    self.baseline.loc[i, 'caa_pn'] = 1
                else:
                    self.baseline.loc[i, 'caa_pn'] = 0

                if row['ccta_aa_degree'] < config.first_analysis.acute_take_off:
                    self.baseline.loc[i, 'caa_angle'] = 1
                else:
                    self.baseline.loc[i, 'caa_angle'] = 0

                if row['ccta_imc_length'] >= config.first_analysis.imc_length_cutoff:
                    self.baseline.loc[i, 'caa_im'] = 1
                else:
                    self.baseline.loc[i, 'caa_im'] = 0

                if (
                    row['caa_im'] == 1
                    and row['ccta_mla_c'] / row['ccta_dist_c'] <= config.first_analysis.percent_stenosis
                ):
                    self.baseline.loc[i, 'caa_hypoplasia'] = 1
                else:
                    self.baseline.loc[i, 'caa_hypoplasia'] = 0

    def process_caa_origin(self, row, i, stenosis_column, high_coronary_column, config):
        """Helper function for correct_bl_hrf"""
        if row[stenosis_column] >= 10:
            self.baseline.loc[i, 'caa_high'] = 1
            self.baseline.loc[i, high_coronary_column] = 1
        else:
            self.baseline.loc[i, 'caa_high'] = 0
            self.baseline.loc[i, high_coronary_column] = 0

    def find_events(data):
        pass

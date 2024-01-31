import os
from typing import Any
import pandas as pd
from loguru import logger


class DataPrepper:
    def __init__(self, config):
        self.config = config
    
    def __call__(self):
        data = pd.read_csv(self.config.data_reader.redcap_database)

        for col in data.columns:
            if 'date' in col or 'year' in col:
                data[col] = pd.to_datetime(data[col], format = r'%Y-%m-%d', errors='coerce')

        # structure is built by redcap_event_name baseline_arm_1 and redcap_repeat_instance is missing
        blueprint = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'].isnull())]

        list_bl = ['diagnostic_exams', 'ccta', 'invasive_testing', 'cmr']
        for i in list_bl:
            blueprint = baseline_arm_1(data, blueprint, redcap_repeat_instrument=i)
        blueprint = surgery_arm_1(data, blueprint)
        blueprint = adverse_events(data, blueprint)
        list_fu = ['1_year_follow_up_arm_1', '3_year_follow_up_arm_1', '5_year_follow_up_arm_1', 'other_follow_up_arm_1']
        for i in list_fu:
            blueprint = follow_up(data, blueprint, redcap_event_name=i)

        blueprint.to_csv(os.path.join(self.config.data_reader.output_dir, 'complete_dataframe.csv'), index=False)
        blueprint.to_excel(os.path.join(self.config.data_reader.output_dir, 'complete_dataframe.xlsx'), index=False)

        return data




def read_data(config):
    """
    Read original csv file, keep original baseline data and save it as a new csv to build
    the structure for the rest of the project.
    """
    data = pd.read_csv(config.data_reader.redcap_database)

    # change every variable containing 'date' to datetime
    for col in data.columns:
        if 'date' in col or 'year' in col:
            data[col] = pd.to_datetime(data[col], format = r'%Y-%m-%d', errors='coerce')

    # structure is built by redcap_event_name baseline_arm_1 and redcap_repeat_instance is missing
    blueprint = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'].isnull())]

    list_bl = ['diagnostic_exams', 'ccta', 'invasive_testing', 'cmr']
    for i in list_bl:
        blueprint = baseline_arm_1(data, blueprint, redcap_repeat_instrument=i)
    blueprint = surgery_arm_1(data, blueprint)
    blueprint = adverse_events(data, blueprint)
    list_fu = ['1_year_follow_up_arm_1', '3_year_follow_up_arm_1', '5_year_follow_up_arm_1', 'other_follow_up_arm_1']
    for i in list_fu:
        blueprint = follow_up(data, blueprint, redcap_event_name=i)

    blueprint.to_csv(os.path.join(config.data_reader.output_dir, 'complete_dataframe.csv'), index=False)
    blueprint.to_excel(os.path.join(config.data_reader.output_dir, 'complete_dataframe.xlsx'), index=False)

    return data


def baseline_arm_1(data, dataframe, redcap_repeat_instrument=None):
    bl = data[
        (data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == redcap_repeat_instrument)
    ]
    
    dict_bl = {'diagnostic_exams': 'dgn_', 'ccta': 'ccta_', 'invasive_testing': 'inv_', 'cmr': 'cmr_'}

    bl_columns = [col for col in bl.columns if col.startswith(dict_bl[redcap_repeat_instrument])]
    add_columns(bl, dataframe, bl_columns, '_')

    return dataframe


def surgery_arm_1(data, dataframe):
    surgery = data[(data['redcap_event_name'] == 'surgery_arm_1')]

    surgery_columns = [col for col in surgery.columns if col.startswith(('surgery_', 'ccta_'))]

    for i, row in surgery.iterrows():
        for col in surgery_columns:
            if col.startswith('ccta_'):
                dataframe.loc[dataframe['record_id'] == row['record_id'], col + '_postop'] = row[col]
            elif col.startswith('surgery_'):
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]

    return dataframe


def adverse_events(data, dataframe):
    ae = data[(data['redcap_event_name'] == 'adverse_event_arm_1')]

    ae_columns = [col for col in ae.columns if col.startswith('ae_')]
    add_columns(ae, dataframe, ae_columns, '_')

    return dataframe


def follow_up(data, dataframe, redcap_event_name = None):
    fu = data[(data['redcap_event_name'] == redcap_event_name)]
    if 'year_follow_up' in redcap_event_name:
        fu['redcap_repeat_instance'] = int(redcap_event_name.split('_')[0])
        decorator = 'y'
        treat_1_different = True
    else:
        decorator = 'other_year_'
        treat_1_different = False

    fu_columns = [col for col in fu.columns if 'fu' in col and not col.startswith('funct_')]
    fu_funct_columns = [col for col in fu.columns if col.startswith('funct_')]

    add_columns(fu, dataframe, fu_columns, '_' + decorator, treat_1_different)
    add_columns(fu, dataframe, fu_funct_columns, '_fu_' + decorator, treat_1_different)

    return dataframe


def add_columns(data, dataframe, columns, suffix_separator, treat_1_different=True):
    for i, row in data.iterrows():
        for col in columns:
            if treat_1_different and data['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif not treat_1_different or data['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + suffix_separator + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe
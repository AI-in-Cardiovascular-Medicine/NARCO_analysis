import os
import pandas as pd
from loguru import logger

def read_data(config):
    """
    Read original csv file, keep original baseline data and save it as a new csv to build
    the structure for the rest of the project.
    """
    data = pd.read_csv(config.data_reader.redcap_database)

    # change every variable containing 'date' to datetime
    for col in data.columns:
        if 'date' in col or 'year' in col:
            data[col] = pd.to_datetime(data[col], format=r'%Y-%m-%d', errors='coerce')

    # structure is built by redcap_event_name baseline_arm_1 and redcap_repeat_instance is missing
    blueprint = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'].isnull())]

    blueprint = ccta_data(data, blueprint)
    blueprint = diagnostics_data(data, blueprint)
    blueprint = invasive_data(data, blueprint)
    blueprint = cmr_data(data, blueprint)
    blueprint = adverse_events(data, blueprint)
    # blueprint = follow_up_1y(data, blueprint)

    blueprint.to_csv(os.path.join(config.data_reader.output_dir, 'complete_dataframe.csv'), index=False)

    return data


def ccta_data(data, dataframe):
    ccta = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == 'ccta')]

    ccta_postop = data[data['redcap_event_name'] == 'surgery_arm_1']

    # make a list with all columns that start with ccta_ in ccta dataframe
    ccta_columns = [col for col in ccta.columns if col.startswith('ccta_')]

    for i, row in ccta.iterrows():
        for col in ccta_columns:
            if ccta['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif ccta['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    ccta_postop_columns = [col for col in ccta_postop.columns if col.startswith(('ccta_', 'surgery_'))]
    for i, row in ccta_postop.iterrows():
        for col in ccta_postop_columns:
            if col.startswith('ccta_'):
                dataframe.loc[dataframe['record_id'] == row['record_id'], col + '_postop'] = row[col]
            elif col.startswith('surgery_'):
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]

    return dataframe

def diagnostics_data(data, dataframe):
    dgn = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == 'diagnostic_exams')]

    # make a list with all columns that start with ccta_ in ccta dataframe
    dgn_columns = [col for col in dgn.columns if col.startswith('dgn_')]

    for i, row in dgn.iterrows():
        for col in dgn_columns:
            if dgn['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif dgn['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe

def invasive_data(data, dataframe):
    inv = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == 'invasive_testing')]

    # make a list with all columns that start with ccta_ in ccta dataframe
    inv_columns = [col for col in inv.columns if col.startswith('inv_')]

    for i, row in inv.iterrows():
        for col in inv_columns:
            if inv['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif inv['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe

def cmr_data(data, dataframe):
    cmr = data[(data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == 'cmr')]

    # make a list with all columns that start with ccta_ in ccta dataframe
    cmr_columns = [col for col in cmr.columns if col.startswith('cmr_')]

    for i, row in cmr.iterrows():
        for col in cmr_columns:
            if cmr['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif cmr['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe

def adverse_events(data, dataframe):
    ae = data[(data['redcap_event_name'] == 'adverse_event_arm_1')]

    # make a list with all columns that start with ccta_ in ccta dataframe
    ae_columns = [col for col in ae.columns if col.startswith('ae_')]

    for i, row in ae.iterrows():
        for col in ae_columns:
            if ae['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif ae['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe

# def follow_up_1y(data, dataframe):
#     fu1 = data[(data['redcap_event_name'] == '1_year_follow_up_arm_1')]

#     # Iterate through columns of fu1
#     for col in fu1.columns:
#         # Check if the column is not empty
#         if not fu1[col].isnull().all():
#             # Add non-empty columns with the suffix '_fu1' to the dataframe
#             new_col_name = col + '_fu1'
#             dataframe[new_col_name] = fu1[col]

#     return dataframe


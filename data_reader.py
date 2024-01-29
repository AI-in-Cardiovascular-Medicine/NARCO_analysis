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

    list_bl = ['diagnostic_exams', 'ccta', 'invasive_testing', 'cmr']
    for i in list_bl:
        blueprint = baseline_arm_1(data, blueprint, redcap_repeat_instrument=i)
    list_fu = ['1_year_follow_up_arm_1', '3_year_follow_up_arm_1', '5_year_follow_up_arm_1', 'other_follow_up_arm_1']
    blueprint = surgery_arm_1(data, blueprint)
    blueprint = adverse_events(data, blueprint)
    for i in list_fu:
        blueprint = follow_up(data, blueprint, event_name=i)

    blueprint.to_csv(os.path.join(config.data_reader.output_dir, 'complete_dataframe.csv'), index=False)

    return data


def baseline_arm_1(data, dataframe, redcap_repeat_instrument=None):
    bl = data[
        (data['redcap_event_name'] == 'baseline_arm_1') & (data['redcap_repeat_instrument'] == redcap_repeat_instrument)
    ]
    
    dict_bl = {'diagnostic_exams': 'dgn_', 'ccta': 'ccta_', 'invasive_testing': 'inv_', 'cmr': 'cmr_'}

    bl_columns = [col for col in bl.columns if col.startswith(dict_bl[redcap_repeat_instrument])]

    for i, row in bl.iterrows():
        for col in bl_columns:
            if bl['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif bl['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

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

    for i, row in ae.iterrows():
        for col in ae_columns:
            if ae['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif ae['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe


def follow_up(data, dataframe, redcap_event_name = None):
    fu = data[(data['redcap_event_name'] == redcap_event_name)]
    fu_columns = [
        'pf_date_fu',
        'pf_days_fu',
        'pf_assessment_fu',
        'pf_funotes',
        'planned_followup_complete',
        'sym_caa_fu',
        'sym_op_fu',
        'sym_op_postopsym_fu___0',
        'sym_op_postopsym_fu___1',
        'sym_op_postopsym_fu___2',
        'sym_op_postopsym_fu___3',
        'sym_fu___0',
        'sym_fu___1',
        'sym_fu___2',
        'sym_fu___3',
        'sym_fu___4',
        'sym_fu___5',
        'sym_fu___6',
        'sym_fu___7',
        'sym_fu___8',
        'sym_fu___9',
        'sym_fu___10',
        'sym_fu___11',
        'sym_fu___12',
        'sym_other_fu',
        'sym_ap_fu',
        'sym_ap_ccs_fu',
        'sym_ap_stress_fu',
        'sym_sync_exercise_fu',
        'sym_sync_stress_fu',
        'sym_dysp_nyha_fu',
        'sym_div_exercise_fu',
        'sym_bp_measure',
        'sym_bp_change_fu',
        'sym_bp_fu',
        'sym_hr_fu',
        'sym_hemodyn_patientreport',
        'sym_comments',
        'symptoms_follow_up_complete',
        'sports_recom_fu',
        'sports_recom_type_fu',
        'sports_recom_yn_fu',
        'sports_yn_fu',
        'sports_level_fu',
        'sports_type_fu',
        'sports_static_component_fu',
        'sports_dynamic_component_fu',
        'sports_hours_fu',
        'sports_years_fu',
        'sports_competitions_fu',
        'sports_work_fu',
        'sports_work_hours_fu',
        'sports_symptoms_fu',
        'sports_symp_type_fu___0',
        'sports_symp_type_fu___1',
        'sports_symp_type_fu___2',
        'sports_symp_type_fu___3',
        'sports_symp_type_fu___4',
        'sports_symp_other_fu',
        'sports_symp_specific_fu',
        'sports_behaviour_follow_up_complete',
    ]

    fu_funct_columns = [col for col in fu.columns if col.startswith('funct_')]

    for i, row in fu.iterrows():
        for col in fu_columns:
            if fu['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif fu['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    for i, row in fu.iterrows():
        for col in fu_funct_columns:
            if fu['redcap_repeat_instance'][i] == 1:
                dataframe.loc[dataframe['record_id'] == row['record_id'], col] = row[col]
            elif fu['redcap_repeat_instance'][i] > 1:
                dataframe.loc[
                    dataframe['record_id'] == row['record_id'], col + '_fu_' + str(int(row['redcap_repeat_instance']))
                ] = row[col]

    return dataframe
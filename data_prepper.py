import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
from loguru import logger
from tqdm import tqdm


class DataPrepper:
    """Creates a longitudinal dataset, by adding repeat_instance columns with a suffix _+ number of repeat instance, 
    and further translating the different arms also to longitudinal format"""

    def __init__(self, config):
        self.config = config
        self.dict_bl = {'diagnostic_exams': 'dgn_', 'ccta': 'ccta_', 'invasive_testing': 'inv_', 'cmr': 'cmr_'}

    def __call__(self):
        self.data = pd.read_csv(self.config.data_prepper.redcap_database)
        date_columns = [col for col in self.data.columns if 'date' in col or 'patient_year' in col]
        self.data[date_columns] = self.data[date_columns].apply(pd.to_datetime, format=r'%Y-%m-%d', errors='coerce')

        # structure is built by redcap_event_name baseline_arm_1 and redcap_repeat_instance is missing
        self.blueprint = self.data[
            (self.data['redcap_event_name'] == 'baseline_arm_1') & (self.data['redcap_repeat_instrument'].isnull())
        ]

        # save blueprint as excel
        self.blueprint.to_excel(
            os.path.join(self.config.data_prepper.output_dir, 'blueprint.xlsx'), index=False
        )

        for instrument in tqdm(self.dict_bl.keys()):
            self.baseline_arm_1(redcap_repeat_instrument=instrument)
        self.surgery_arm_1()
        self.adverse_events()
        list_fu = [
            '1_year_follow_up_arm_1',
            '3_year_follow_up_arm_1',
            '5_year_follow_up_arm_1',
            'other_follow_up_arm_1',
        ]
        for i in tqdm(list_fu):
            self.follow_up(redcap_event_name=i)

        # self.blueprint.to_csv(os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.csv'), index=False)
        self.blueprint.to_excel(
            os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.xlsx'), index=False
        )

        # # stupid fix since I didn't find the mistake yet
        # self.blueprint = self.update_blueprint_from_complete()
        
        # self.blueprint.to_excel(
        #     os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.xlsx'), index=False
        # )
        return self.blueprint

    def baseline_arm_1(self, redcap_repeat_instrument=None):
        bl = self.data[
            (self.data['redcap_event_name'] == 'baseline_arm_1')
            & (self.data['redcap_repeat_instrument'] == redcap_repeat_instrument)
        ]

        bl_columns = [col for col in bl.columns if col.startswith(self.dict_bl[redcap_repeat_instrument])]
        self.add_columns(bl, bl_columns, '_')

    def surgery_arm_1(self):
        surgery = self.data[(self.data['redcap_event_name'] == 'surgery_arm_1')]

        surgery_columns = [col for col in surgery.columns if col.startswith(('surgery_', 'ccta_'))]

        for i, row in surgery.iterrows():
            for col in surgery_columns:
                if col.startswith('ccta_'):
                    self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col + '_postop'] = row[col]
                elif col.startswith('surgery_'):
                    self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col] = row[col]

    def adverse_events(self):
        ae = self.data[(self.data['redcap_event_name'] == 'adverse_event_arm_1')]

        ae_columns = [col for col in ae.columns if col.startswith('ae_')]
        self.add_columns(ae, ae_columns, '_')

    def follow_up(self, redcap_event_name=None):
        fu = self.data[(self.data['redcap_event_name'] == redcap_event_name)]
        if 'year_follow_up' in redcap_event_name:
            fu['redcap_repeat_instance'] = int(redcap_event_name.split('_')[0])
            decorator = 'y'
            treat_1_different = True
        else:
            decorator = 'other_y'
            treat_1_different = False

        fu_columns = [col for col in fu.columns if 'fu' in col and not col.startswith('funct_')]
        fu_funct_columns = [col for col in fu.columns if col.startswith('funct_')]

        self.add_columns(fu, fu_columns, '_' + decorator, treat_1_different)
        self.add_columns(fu, fu_funct_columns, '_fu_' + decorator, treat_1_different)

    def add_columns(self, data, columns, suffix_separator, treat_1_different=True): 
        for i, row in data.iterrows(): 
            for col in columns: 
                if treat_1_different and data['redcap_repeat_instance'][i] == 1: 
                    self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col] = row[col] 
                else: 
                    self.blueprint.loc[ 
                        self.blueprint['record_id'] == row['record_id'], 
                        col + suffix_separator + str(int(row['redcap_repeat_instance'])), 
                    ] = row[col] 

    # def add_columns(self, data, columns, suffix_separator, treat_1_different=True):
    #     exception_columns = ['hx_', 'sports_', 'caa_', 'funct_', 'synp_']

    #     for i, row in data.iterrows():
    #         for col in columns:
    #             # Check if the column starts with any of the exception columns
    #             if any(col.startswith(prefix) for prefix in exception_columns):
    #                 # For the exception columns, we ALWAYS add a suffix, regardless of the instance
    #                 if pd.notnull(row[col]):
    #                     suffix = suffix_separator + str(int(row['redcap_repeat_instance']))
    #                     self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col + suffix] = row[col]
    #             else:
    #                 # For non-exception columns (including inv_ and others)
    #                 if treat_1_different and data['redcap_repeat_instance'][i] == 1:
    #                     # For the first instance, we do not add a suffix
    #                     self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col] = row[col]
    #                 elif pd.notnull(row[col]):
    #                     # For repeated instances, we add a suffix based on redcap_repeat_instance
    #                     suffix = suffix_separator + str(int(row['redcap_repeat_instance']))
    #                     self.blueprint.loc[self.blueprint['record_id'] == row['record_id'], col + suffix] = row[col]

    def update_blueprint_from_complete(self):
        """Reads blueprint.xlsx and complete_dataframe.xlsx, and fills missing values in blueprint from complete_dataframe"""
        blueprint_path = os.path.join(self.config.data_prepper.output_dir, 'blueprint.xlsx')
        complete_dataframe_path = os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.xlsx')

        # Load the blueprint and complete dataframe
        blueprint_df = pd.read_excel(blueprint_path)
        complete_df = pd.read_excel(complete_dataframe_path)

        # Ensure they align on 'record_id' or the primary key column
        common_cols = [col for col in blueprint_df.columns if col in complete_df.columns and col != 'record_id']

        # Fill missing values in blueprint using values from complete_df
        for col in common_cols:
            blueprint_df[col] = blueprint_df[col].fillna(complete_df[col])

        # Save the updated blueprint
        blueprint_df.to_excel(blueprint_path, index=False)
        logger.info(f"Updated blueprint saved to {blueprint_path}")

        return blueprint_df
        
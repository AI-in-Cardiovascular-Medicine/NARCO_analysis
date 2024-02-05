import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
from loguru import logger


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

        for instrument in self.dict_bl.keys():
            self.baseline_arm_1(redcap_repeat_instrument=instrument)
        self.surgery_arm_1()
        self.adverse_events()
        list_fu = [
            '1_year_follow_up_arm_1',
            '3_year_follow_up_arm_1',
            '5_year_follow_up_arm_1',
            'other_follow_up_arm_1',
        ]
        for i in list_fu:
            self.follow_up(redcap_event_name=i)

        # self.blueprint.to_csv(os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.csv'), index=False)
        self.blueprint.to_excel(
            os.path.join(self.config.data_prepper.output_dir, 'complete_dataframe.xlsx'), index=False
        )

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

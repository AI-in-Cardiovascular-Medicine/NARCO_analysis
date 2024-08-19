import os

import pandas as pd
from loguru import logger


class Demographics:
    """..."""

    def __init__(self, config):
        self.config = config
        # keep column names only
        self.dict = pd.read_csv(config.demographics_table.dict)
        print(self.dict)

    def __call__(self):
        self.data = pd.read_excel(self.config.first_analysis.input_file)
        self.baseline = self.keep_bl_columns()
        self.correct_baseline_inv()
        self.correct_bl_hrf(self.config)

        self.follow_up = self.keep_fu_columns()
        self.find_events()

        self.baseline.to_excel(os.path.join(self.config.first_analysis.output_dir, 'baseline.xlsx'), index=False)
        self.follow_up.to_excel(os.path.join(self.config.first_analysis.output_dir, 'follow_up.xlsx'), index=False)

        return self.baseline, self.follow_up
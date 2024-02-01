import pandas as pd
import numpy as np
from loguru import logger

test = pd.read_csv('L1958Narco_DATA_2024-01-31_1014.csv')

final = pd.read_csv('C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/dataframes/complete_dataframe.csv')

filtered1 = final.filter(regex='^(?!.*_\d+)')
filtered2 = filtered1.filter(regex='^(?!.*(_1|\_\d+)$).*$')

cols_to_exclude = [
            col
            for col in filtered2.columns
            if col.startswith(('ae_', 'adverse_', 'end_', 'do_', 'drop_', 'sign_', 'signature_')) or col.endswith('_fu')
        ]

baseline = filtered2.drop(columns=cols_to_exclude)

inv = final[[col for col in final.columns if col.startswith(('record_', 'inv_'))]]

keys = range(1, len(final) + 1)

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

# replace for values 'inv_protocol___2' with ''
for key, value in dict_inv.items():
    if value is not None:
        dict_inv[key] = value.replace('inv_protocol___2', '')

dict_inv = {key: (value if value != '' else None) for key, value in dict_inv.items()}

# create subset from final with only record_id == 1, 2, 6, 10, 13, 14, 15, 32, 41, 58, 66, 68
test = final[final['record_id'].isin([1, 2, 6, 10, 13, 14, 15, 32, 41, 58, 66, 68])]
test = test[[col for col in test.columns if col.startswith(('record_', 'inv_'))]]

print(test)
print(dict_inv)

for i, row in test.iterrows():
    if row['record_id'] in dict_inv.keys() and dict_inv[row['record_id']] is not None:
        suffix = dict_inv[row['record_id']]
        # replace all columns starting with 'inv_' in self.baselin with all columns in inv starting with 'inv_' and ending with value from dict_inv
        for col in test.columns:
            if col.startswith('inv_'):
                test.loc[i, col] = inv.loc[i, str(col + dict_inv.values())]

print(test)
'''
Created on 29/08/2013

@author: dlawrence
'''

import pandas as pd

PIRNA_MIN_SIZE = 23
PIRNA_MAX_SIZE = 32

WILDTYPE = "WT"
MUTANT = "MUTANT"

CONDITION_NAMES = {"Mut" : "pin/pin"}
SAMPLE_NAMES = {"F2" : "Spermatocytes", "F5" : "Round spermatids"}

SAMPLE_TYPE = "sample_type"
CONDITION = "condition"
GROUP_BY_COLUMNS = [SAMPLE_TYPE, CONDITION]

OTHER = "other"
PRIMARY_PIRNA = "p-piRNA"
SECONDARY_PIRNA = "s-piRNA"
PIRNA_TYPES = [PRIMARY_PIRNA, SECONDARY_PIRNA]
CLASSIFICATION_TYPES = [OTHER] + PIRNA_TYPES

def get_conditions_sample_types_reps(index):
    condition_series = pd.Series(index=index, dtype=str)
    sample_type_series = pd.Series(index=index, dtype=str)
    rep_series = pd.Series(index=index, dtype=str)
         
    for i in index:
        name = i.split('.')[0] # Looks like: 'SL_Mut_F2_Rep2'
        (condition, sample_type, rep) = name.split('_')[1:4]
        
        condition_series[i] = CONDITION_NAMES.get(condition, condition)
        sample_type_series[i] = SAMPLE_NAMES.get(sample_type, sample_type)
        rep_series[i] = rep
    return (condition_series, sample_type_series, rep_series)

def add_conditions_and_sample_types(df):
    (conditions, sample_types, _) = get_conditions_sample_types_reps(df.index)
    df[CONDITION] = conditions
    df[SAMPLE_TYPE] = sample_types
    return df

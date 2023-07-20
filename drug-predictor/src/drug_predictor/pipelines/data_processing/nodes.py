import numpy as np
import pandas as pd
import json
import os
from typing import List, Dict

def drop_columns(x: pd.DataFrame, parameters: List) -> pd.DataFrame:
    x = x.drop(parameters, axis=1)
    return x

def rename_columns(x: pd.DataFrame, parameters: Dict) -> pd.DataFrame:
    x = x.rename(columns=parameters)
    return x

def is_lipinski(x: pd.DataFrame) -> pd.DataFrame:
    hdonor = x['HBondDonorCount'] < 6
    haccept = x['HBondAcceptorCount'] < 10
    mw = x['MolecularWeight'] < 500
    clogP = x['LogP'] < 5
    # Apply rules to dataframe
    x['RuleFive'] = np.where(((hdonor & haccept & mw) | (hdonor & haccept & clogP) | (hdonor & mw & clogP) | (haccept & mw & clogP)), 1, 0)
    return x

def matc_conversion(x: pd.DataFrame, parameters: Dict) -> pd.DataFrame:
    x['MATC_Code_Short'] = x['drug_class'].map(parameters)
    return x

def matc_explanation(x: pd.Series, parameters: Dict) -> pd.Series:
    x = x.map(parameters)
    return x

def process_gitter(gitter: pd.DataFrame, columns: Dict, conversions: Dict, explanations: Dict) -> pd.DataFrame:
    """ Processes the data for gitter dataset.
    Args:
      gitter: Raw data.
    Returns:
      Processed dataset without irrelevant columns, drug_class changed to MATC_Code_Short, and added columns RuleFive and MATC_conversion"""
    gitter = drop_columns(gitter, columns['irrelevant_columns'])
    gitter = rename_columns(gitter, columns['columns_to_rename'])
    gitter = is_lipinski(gitter)
    gitter = matc_conversion(gitter, conversions)
    gitter = drop_columns(gitter, columns['columns_to_drop'])
    gitter['MATC_Code_Explanation'] = matc_explanation(gitter['MATC_Code_Short'], explanations)
    return gitter



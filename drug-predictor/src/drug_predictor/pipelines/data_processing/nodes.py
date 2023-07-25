import requests as rq
import json
import re
import numpy as np
import pandas as pd
from typing import List, Dict
import pubchempy as pcp

def get_inchikey(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'InChIKey':
                return dic['value']  
    except:
        return None
    
def get_h_bond_accept(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'H Bond Acceptor Count':
                return dic['value']
    except:        
        return None
    
def get_h_bond_donor(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'H Bond Donor Count':
                return dic['value']
    except:        
        return None
    
def get_mol_weight(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'Molecular Weight':
                return dic['value']
    except:        
        return None
    
def get_logp(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'logP' and dic['source'] == 'ChemAxon':
                return dic['value']
    except:        
        return None

def get_rule_five(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'Rule of Five':
                return dic['value']
    except:        
        return None
    
def get_isomeric_smiles(x: Dict) -> pd.Series:
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'SMILES':
                return dic['value']
    except:        
        return None
    
def get_atc_code(x: Dict) -> pd.Series:
    try:
        for code in x['atc-codes']['atc-code']:
            return code['@code']
    except:
        return None
    
def shorten_atc_code(x: pd.Series) -> pd.Series:
    x = x.str[0]
    return x

def select_columns(x: pd.DataFrame, parameters: List) -> pd.DataFrame:
    x = x[parameters]
    return x

def get_atc_code(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
    response = rq.get(url)
    pat = r"\"[A-Z]\d{2}[A-Z]{2}\d{2}\""
    atc_code_found = re.search(pat, response.text)
    if atc_code_found:
        return atc_code_found.group(0).strip('"')

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

def process_pubchem(pubchem: pd.DataFrame, columns: Dict, explanations: Dict) -> pd.DataFrame:
    """ Processes the data for pubchem dataset.
    Args:
      pubchem: Raw data.
    Returns:
      Processed dataset without irrelevant columns, ATC_Code changed to MATC_Code_Short, and added columns RuleFive and MATC_conversion"""
    pubchem = select_columns(pubchem, columns['columns_to_select'])
    pubchem = rename_columns(pubchem, columns['columns_to_rename'])
    pubchem = is_lipinski(pubchem)
    pubchem['ATC_Code'] = pubchem['CID'].map(get_atc_code)
    pubchem['MATC_Code_Short'] = shorten_atc_code(pubchem['ATC_Code'])
    pubchem['MATC_Code_Explanation'] = matc_explanation(pubchem['MATC_Code_Short'], explanations)
    pubchem = drop_columns(pubchem, columns['columns_to_drop'])
    return pubchem

def preprocess_drugbank(drugbank: List) -> pd.DataFrame:
    """ Processes the data for drugbank dataset.
    Args:
      drugbank: Raw data in xml format
    Returns:
      Processed dataset without irrelevant columns, ATC_Code changed to MATC_Code_Short, and added columns RuleFive and MATC_conversion"""
    
    drugs_df = pd.DataFrame()
    drugs_df['InChIKey'] = list(map(get_inchikey, drugbank))
    drugs_df['HBondAcceptorCount'] = list(map(get_h_bond_accept, drugbank))
    drugs_df['HBondDonorCount'] = list(map(get_h_bond_donor, drugbank))
    drugs_df['MolecularWeight'] = list(map(get_mol_weight, drugbank))
    drugs_df['LogP'] = list(map(get_logp, drugbank))
    drugs_df['RuleFive'] = list(map(get_rule_five, drugbank))
    drugs_df['IsomericSMILES'] = list(map(get_isomeric_smiles, drugbank))
    drugs_df['ATC_Code'] = list(map(get_atc_code, drugbank))
    return drugs_df

def 


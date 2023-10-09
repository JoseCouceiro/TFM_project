import requests as rq
import re
import numpy as np
import pandas as pd
from typing import List, Dict
import pubchempy as pcp

# Processing Drugbank database

def get_small_molecules(x: list) -> List:
    """
    Function that extracts the small molecules from the Drugbank database and discards the biomolecules.
    Input: list of molecules.
    Output: list of small molecules.
    """
    small_drug_list = []
    for drug in x:
        if drug['@type'] == 'small molecule':
            small_drug_list.append(drug)
    return small_drug_list

def get_inchikey(x: Dict) -> pd.Series:
    """
    Function that searches the InChiKey value of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: string.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'InChIKey':
                return dic['value']  
    except:
        return None
    
def get_h_bond_accept(x: Dict) -> pd.Series:
    """
    Function that searches the number of H bond acceptors of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: integer.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'H Bond Acceptor Count':
                return dic['value']
    except:        
        return None
    
def get_h_bond_donor(x: Dict) -> pd.Series:
    """
    Function that searches the number of H bond donors of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: integer.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'H Bond Donor Count':
                return dic['value']
    except:        
        return None
    
def get_mol_weight(x: Dict) -> pd.Series:
    """
    Function that searches the molecular weight of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: float.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'Molecular Weight':
                return dic['value']
    except:        
        return None
    
def get_logp(x: Dict) -> pd.Series:
    """
    Function that searches the partition coefficient of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: integer.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'logP' and dic['source'] == 'ChemAxon':
                return dic['value']
    except:        
        return None

def get_rule_five(x: Dict) -> pd.Series:
    """
    Function that searches the rule of five value for a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: integer.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'Rule of Five':
                return dic['value']
    except:        
        return None
    
def get_isomeric_smiles(x: Dict) -> pd.Series:
    """
    Function that searches the SMILES of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: string.
    """
    try:
        for dic in x['calculated-properties']['property']:
            if dic['kind'] == 'SMILES':
                return dic['value']
    except:        
        return None
    
def get_atc_code_drugbank(x: Dict) -> pd.Series:
    """
    Function that searches the ATC code of a given molecule in the Drugbank database and returns it.
    Input: dictionary.
    Output: string.
    """
    try:
        for code in x['atc-codes']['atc-code']:
            return code['@code']
    except:
        return None
    
def get_cid_from_smiles(smiles: str) -> int:
    """
    Function uses Pubchempy function "get_compounds" to obtain the CID number of a molecule from its SMILES.
    Input: string.
    Output: integer.
    """
    try:
        mol = pcp.get_compounds(smiles, 'smiles')
        return int(mol[0].cid)
    except:
        print(f'no cid for smiles {smiles}')
        return None
    
def shorten_atc_code(x: pd.Series) -> pd.Series:
    """
    Function that shortens the ATC code to the first letter.
    Input: string.
    Output: string.
    """
    x = x.str[0]
    return x

def select_columns(x: pd.DataFrame, parameters: List) -> pd.DataFrame:
    """
    Function that selects a set of defined columns from a pandas dataframe and returns a shortened dataframe.
    Input: pandas dataframe, list of column names.
    Output: pandas dataframe.
    """
    x = x[parameters]
    return x

# Processing Pubchem database

def get_atc_code(cid: int) -> str:
    """
    Function that gets the ATC code of a given molecule using its cid from the Pubchem database using "requests".
    Input: integer.
    Output: string.
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
    response = rq.get(url)
    pat = r"\"[A-Z]\d{2}[A-Z]{2}\d{2}\""
    atc_code_found = re.search(pat, response.text)
    if atc_code_found:
        return atc_code_found.group(0).strip('"')

def drop_columns(x: pd.DataFrame, parameters: List) -> pd.DataFrame:
    """
    Function that drops a set of defined columns from a pandas dataframe and returns the resulting dataframe.
    Input: pandas dataframe, list o column names.
    Output: pandas dataframe.
    """
    x = x.drop(parameters, axis=1)
    return x

def rename_columns(x: pd.DataFrame, parameters: Dict) -> pd.DataFrame:
    """
    Function that renames a set of defined columns from a pandas dataframe and returns the resulting dataframe.
    Input: pandas dataframem, dictionary.
    Output: pandas dataframe.
    """
    x = x.rename(columns=parameters)
    return x

def is_lipinski(x: pd.DataFrame) -> pd.DataFrame:
    """
    Function that applies a set of rules (Lipinski rules) to several columns of a pandas dataframe and returns \
          a dataframe with a new column that states if said rules were passed or not.
    Input: pandas dataframe.
    Output: pandas dataframe.
    """
    # Lipinski rules
    hdonor = x['HBondDonorCount'] <= 5
    haccept = x['HBondAcceptorCount'] <= 10
    mw = x['MolecularWeight'] < 500
    clogP = x['LogP'] <= 5
    # Apply rules to dataframe
    x['RuleFive'] = np.where(((hdonor & haccept & mw) | (hdonor & haccept & clogP) | (hdonor & mw & clogP) | (haccept & mw & clogP)), 1, 0)
    return x

# Processing Pubchem dataset

def matc_conversion(x: pd.DataFrame, parameters: Dict) -> pd.DataFrame:
    """
    Function that transforms the label data "drug_class" into the standard for every dataframe ("MATC_Code_Short") using a dictionary.
    Input: pandas dataframe.
    Output: pandas dataframe.
    """
    x['MATC_Code_Short'] = x['drug_class'].map(parameters)
    return x

def matc_explanation(x: pd.Series, parameters: Dict) -> pd.Series:
    """
    Function that adds a column with the explanation of the MATC codes to a dataframe using a dictionary.
    Input: pandas dataframe, dictionary.
    Output: pandas dataframe.
    """
    x = x.map(parameters)
    return x

# Combined functions

def process_gitter(gitter: pd.DataFrame, columns: Dict, conversions: Dict, explanations: Dict) -> pd.DataFrame:
    """ Functions that processes the data from the gitter dataset.
    Args:
      gitter: raw Gitter dataset.
      columns: list of columns to drop and rename arranged in a dictionary.
      conversions: a dictionary with label conversions.
      explanations: a dictionary with the explanations of the MATC code.
    Returns:
      Processed dataset without irrelevant columns, "drug_class" changed to "MATC_Code_Short", \
        and added columns "RuleFive" and "MATC_conversion".
    """
    gitter = drop_columns(gitter, columns['irrelevant_columns'])
    gitter = rename_columns(gitter, columns['columns_to_rename'])
    gitter = is_lipinski(gitter)
    gitter = matc_conversion(gitter, conversions)
    gitter = drop_columns(gitter, columns['columns_to_drop'])
    gitter['MATC_Code_Explanation'] = matc_explanation(gitter['MATC_Code_Short'], explanations)
    return gitter

def preprocess_pubchem(pubchem: pd.DataFrame, columns: Dict) -> pd.DataFrame:
    """ Function that processes the data from the pubchem dataset.
    Args:
      pubchem: raw Pubchem data.
      columns: list of columns to keep and rename.
    Returns:
      Processed dataset without irrelevant columns and added columns "RuleFive" and "ATC_Code"
    """
    pubchem = select_columns(pubchem, columns['columns_to_select'])
    pubchem = rename_columns(pubchem, columns['columns_to_rename'])
    pubchem = is_lipinski(pubchem)
    pubchem['ATC_Code'] = pubchem['CID'].map(get_atc_code)
    pubchem = pubchem[pubchem['ATC_Code'].isna() == False]
    return pubchem

def process_pubchem(pubchem: pd.DataFrame, columns: Dict, explanations: Dict) -> pd.DataFrame:
    """ Function that further processes the data from the pubchem dataset.
    Args:
      pubchem: prerocessed Pubchem dataset.
      columns: list of columns to drop.
      explanations: dictionary with the explanations of the MATC_Code.
    Returns:
      Processed dataset without irrelevant columns, "ATC_Code" changed to "MATC_Code_Short", and added column "MATC_Code_Explanation".
    """
    pubchem['MATC_Code_Short'] = shorten_atc_code(pubchem['ATC_Code'])
    pubchem['MATC_Code_Explanation'] = matc_explanation(pubchem['MATC_Code_Short'], explanations)
    pubchem = drop_columns(pubchem, columns['columns_to_drop'])
    return pubchem

def preprocess_drugbank(drugbank: List) -> pd.DataFrame:
    """ Function that processes the data from the drugbank dataset.
    Args:
      drugbank: raw Drugbank data.
    Returns:
      Processed dataset with selected information
    """
    drugbank = get_small_molecules(drugbank)
    drugs_df = pd.DataFrame()
    drugs_df['InChIKey'] = list(map(get_inchikey, drugbank))
    drugs_df['HBondAcceptorCount'] = list(map(get_h_bond_accept, drugbank))
    drugs_df['HBondDonorCount'] = list(map(get_h_bond_donor, drugbank))
    drugs_df['MolecularWeight'] = list(map(get_mol_weight, drugbank))
    drugs_df['LogP'] = list(map(get_logp, drugbank))
    drugs_df['RuleFive'] = list(map(get_rule_five, drugbank))
    drugs_df['IsomericSMILES'] = list(map(get_isomeric_smiles, drugbank))
    drugs_df['ATC_Code'] = list(map(get_atc_code_drugbank, drugbank))
    return drugs_df

def process_drugbank(drugbank: pd.DataFrame, pubchem:pd.DataFrame, columns: Dict, explanations: Dict) -> pd.DataFrame:
    """ Function that further processes the data from the drugbank dataset.
    Args:
      drugbank: preprocessed Drugbank dataset.
      pubchem: processed Pubchen dataset.
      columns: list of columns to drop.
      explanations: dictionary with the explanations of the MATC_Code.
    Returns:
      Processed dataset containing only molecules that have an ATC_Code and SMILES and are not present in the pubchem dataset. \
      Irrelevant columns are dropped and "MATC_Code_Short", "CID" and "MATC_Code_Explantion" are added.
    """
    drugbank = drugbank[drugbank['ATC_Code'].isna()==False]
    drugbank = drugbank[drugbank['IsomericSMILES'].isna()==False]
    drugbank = drugbank[drugbank['ATC_Code'].isin(pubchem['ATC_Code'])==False]
    drugbank['MATC_Code_Short'] = shorten_atc_code(drugbank['ATC_Code'])
    drugbank['CID'] = drugbank['IsomericSMILES'].map(get_cid_from_smiles)
    drugbank = drugbank[drugbank['CID'].isna()==False]
    drugbank = drop_columns(drugbank, columns['columns_to_drop'])
    drugbank['MATC_Code_Explanation'] = matc_explanation(drugbank['MATC_Code_Short'], explanations)
    return drugbank

def join_datasets(pubchem: pd.DataFrame, drugbank: pd.DataFrame, gitter: pd.DataFrame) -> pd.DataFrame:
    """ Function joins three datasets together. To avoid duplications, it only concatenates the frames after duplicated CIDs are\
    eliminated from a given dataset by using the method 'isin'.
    Args:
      drugbank: processed Drugbank dataset.
      pubchem: processed Pubchen dataset.
      gitter: proceseed Gitter dataset.
    Returns:
      Processed dataset containing a non-redundant set of molecules from the three processed datasets.
    """
    drugbank['CID'] = drugbank['CID'].astype('int')
    drugbank = drugbank[drugbank['CID'].isin(pubchem['CID'])==False]
    pubchem_drugbank = pd.concat([pubchem, drugbank])
    gitter = gitter[gitter['CID'].isin(pubchem_drugbank['CID'])==False]
    all_drugs = pd.concat([gitter, pubchem_drugbank])
    return all_drugs
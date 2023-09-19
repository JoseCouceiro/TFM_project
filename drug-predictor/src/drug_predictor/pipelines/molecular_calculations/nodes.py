import pandas as pd
import numpy as np
from sklearn import preprocessing
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import PandasTools as pt
import pubchempy as pcp
from typing import List, Dict

# Functions to obtain fingerprints

def compute_connectivity_invariants(mol) -> np.array:
    """Function that obtains the connectivity invariants of a molecule.
    Input: RDKit molecule object.
    Output: numpy array.
    """
    try:
        con_inv_fp = rdMolDescriptors.GetConnectivityInvariants(mol)
    except:
        print('Something went wrong computing Connectivity Invariants')
        return None
    return np.array(con_inv_fp)

def compute_feature_invariants(mol) -> np.array:
    """Function that obtains the feature invariants of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        inv_fp = rdMolDescriptors.GetFeatureInvariants(mol)
    except:
        print('Something went wrong computing Feature Invariants')
        return None
    return np.array(inv_fp)

def compute_morgan_fp(mol, depth=2, nBits=2048) -> np.array:
    """Function that obtains the Morgan fingerprints of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        mor_fp = AllChem.GetMorganFingerprintAsBitVect(mol,depth,nBits)
    except:
        print('Something went wrong computing Morgan fingerprints')
        return None
    return np.array(mor_fp)

def compute_maccskeys(mol) -> np.array:
    """Function that obtains the MACCSKeys of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        mkeys = MACCSkeys.GenMACCSKeys(mol)   
    except:
        print('Something went wrong computing MACCSKeys')
        return None
    return np.array(mkeys)

def compute_atom_pair_fp(mol, nBits=2048) -> np.array:
    """Function that obtains the atom pair Fingerprints of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        atom_pair_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits)
    except:
        print('Something went wrong computing Atom Pair fingerprints')
        return None
    return np.array(atom_pair_fp)

def compute_topological_torsion_fp(mol) -> np.array:
    """Function that obtains the topological torsion fingerprints of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        tt_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol)
    except:
        print('Something went wrong computing Topological Torsion fingerprints')
        return None
    return np.array(tt_fp)

def compute_avalon_fp(mol, nBits=2048) -> np.array:
    """Function that obtains the Avalon fingerprints of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        av_fp = pyAvalonTools.GetAvalonFP(mol, nBits)
    except:
        print('Something went wrong computing Avalon fingerprints')
        return None
    return np.array(av_fp)

def compute_rdkit_fp(mol, maxPath=5, fpSize=2048) -> np.array:
    """Function that obtains the RDKit fingerprints of a molecule.
    Input: RDKit molecule.
    Output: numpy array.
    """
    try:
        rdkit_fp = AllChem.RDKFingerprint(mol, maxPath, fpSize)
    except:
        print('Something went wrong computing RDKit fingerprints')
        return None
    return np.array(rdkit_fp)

def compute_pubchem_fingerprints(cid: int) -> np.array:
    """Function that obtains the PubChem fingerprints of a molecule.
    Input: molecules's CID.
    Output: numpy array.
    """
    try:
        comp = pcp.Compound.from_cid(int(cid))
        fp_bin = bin(int(comp.fingerprint, 16))[2:]   
    except:
        print('Something went wrong computing Pubchem fingerprints')
        return None
    return np.array(list(fp_bin)).astype('int')

""" def compute_cactvs_fingerprints(cid: int) -> np.array:
    Function that obtains the Cactvs fingerprints of a molecule.
    Input: molecule's CID.
    Output: numpy array.
    
    try:
        comp = pcp.Compound.from_cid(int(cid))
        cactvs_fp_bin = bin(int(comp.fingerprint, 16))[2:]
    except:
        print('Something went wrong computing Cactvs fingerprints')
        return None
    return np.array(list(cactvs_fp_bin)).astype('int') """

# Getting RDKit molecules

def add_molecule_column(all_drugs: pd.DataFrame, columns: Dict) -> pd.DataFrame:
    """
    Function that adds a "Molecule" column to a dataframe contining the RDKit molecule objet of each entry.
    Input: pandas dataframe, dictionary containing columns.
    Output: pandas dataframe.
    """
    base_column = columns['base']
    calculated_column = columns['calculated']
    pt.AddMoleculeColumnToFrame(frame=all_drugs, smilesCol=base_column, molCol=calculated_column)
    return all_drugs

# Feature engineering to obtain the input table

def drop_columns(all_drugs: pd.DataFrame, columns: Dict) -> pd.DataFrame:
    """
    Function that drops unnecessary columns from a dataframe. Columns are set in the parameters.
    Input: pandas dataframe, dictionary containing column names.
    Output: pandas dataframe.
    """
    all_drugs = all_drugs.drop(columns['columns_to_drop'], axis=1)
    return all_drugs

def label_encoder(all_drugs: pd.DataFrame) -> pd.DataFrame:
    """
    Function that replaces the values of the "Label" column, that contains alphabetical characters, for integers.
    Input: pandas dataframe.
    Output: pandas dataframe.
    """
    le = preprocessing.LabelEncoder()
    all_drugs['Label'] = le.fit_transform(all_drugs['MATC_Code_Short'])
    return all_drugs

def extract_code_to_label_dic(all_drugs: pd.DataFrame) -> Dict:
    """
    Function that obtains a dictionary of the label numerical values as keys and the code explanations as values.
    Input: pandas dataframe.
    Output: dictionary.
    """
    code_to_label_group = all_drugs.groupby(['Label', 'MATC_Code_Explanation']).count().reset_index()
    code_to_label_dataset = code_to_label_group[['Label', 'MATC_Code_Explanation']].set_index('Label')
    code_to_label_dic = code_to_label_dataset.to_dict()['MATC_Code_Explanation']
    return code_to_label_dic

# Building a dataframe of fingerprints

def get_fingerprints(all_drugs: pd.DataFrame) -> pd.DataFrame:
    """
    Function that takes RDKit mol objects from a pandas dataframe column and, by mapping the functions described above, \
        generates from it different columns containing fingerprints. 
    Input: pandas dataframe containing a 'Molecule' column.
    Output: pandas dataframe containing several fingerprints columns.
    """
    #all_drugs['FeatInvariants'] = all_drugs['Molecule'].map(compute_feature_invariants)
    #all_drugs['ConnInvariants'] = all_drugs['Molecule'].map(compute_connectivity_invariants)
    all_drugs['Morgan2FP'] = all_drugs['Molecule'].map(compute_morgan_fp)
    all_drugs['MACCSKeys'] = all_drugs['Molecule'].map(compute_maccskeys)
    all_drugs['AtomPairFP'] = all_drugs['Molecule'].map(compute_atom_pair_fp)
    all_drugs['TopTorFP'] = all_drugs['Molecule'].map(compute_topological_torsion_fp)
    all_drugs['AvalonFP'] = all_drugs['Molecule'].map(compute_avalon_fp)
    
    all_drugs['PubchemFP']= all_drugs['CID'].map(compute_pubchem_fingerprints) #This takes over 1 hour in my computer
    #all_drugs['CactvsFP']= all_drugs['CID'].map(compute_cactvs_fingerprints) #This takes over 1 hour in my computer
    #all_drugs['RDKitFP']= all_drugs['Molecule'].map(compute_rdkit_fp) #This takes so long that crashes my computer, but I coudn't find a way around
    
    return all_drugs

def clean_dataset(all_drugs: pd.DataFrame, selected_columns: Dict) -> pd.DataFrame:
    """
    Function that drops unnecessary columns from a pandas dataset and encodes the 'Label' column as an integer. \
        The unnecesary columns are set in the parameters. Finally, drop entries which may have not worked correctly \
        using pandas' funtion "dropna()".
    Input: pandas dataframe, dictionary containing column names.
    Output: pandas dataframe.
    """
    all_drugs = drop_columns(all_drugs, selected_columns)
    all_drugs = label_encoder(all_drugs)
    all_drugs = all_drugs.dropna()
    return all_drugs

# Combined functions

def get_model_input(all_drugs: pd.DataFrame, columns: Dict, selected_columns: Dict) -> pd.DataFrame:
    """
    Function build a dataframe containing molecule descriptors, their corresponding fingerprints and a numerical label.
    Args:
      all_drugs: processed dataframe.
      columns: list of column names convert from descriptor to RDKit molecule.
      selected_columns: list of column names to drop.
    Output: pandas dataframe, dictionary of numerical labels to code explanations.
    """
    all_drugs = add_molecule_column(all_drugs, columns)
    all_drugs = get_fingerprints(all_drugs)
    all_drugs = clean_dataset(all_drugs, selected_columns)
    code_to_label_dic = extract_code_to_label_dic(all_drugs)
    return all_drugs, code_to_label_dic
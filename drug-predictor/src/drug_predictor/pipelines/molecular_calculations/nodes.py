import pandas as pd
import numpy as np
from sklearn import preprocessing
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import PandasTools as pt
import pubchempy as pcp
from typing import List, Dict, AnyStr

#Functions to obtain fingerprints
def compute_connectivity_invariants(mol):
    """Function that obtains the connectivity invariants of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        con_inv_fp = rdMolDescriptors.GetConnectivityInvariants(mol)
    except:
        print('Something went wrong computing Connectivity Invariants')
        return None
    return np.array(con_inv_fp)

def compute_feature_invariants(mol):
    """Function that obtains the feature invariants of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        inv_fp = rdMolDescriptors.GetFeatureInvariants(mol)
    except:
        print('Something went wrong computing Feature Invariants')
        return None
    return np.array(inv_fp)

def compute_morgan_fp(mol, depth=2, nBits=2048):
    """Function that obtains the Morgan fingerprints of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        mor_fp = AllChem.GetMorganFingerprintAsBitVect(mol,depth,nBits)
    except:
        print('Something went wrong computing Morgan fingerprints')
        return None
    return np.array(mor_fp)

def compute_maccskeys(mol):
    """Function that obtains the MACCSKeys of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        mkeys = MACCSkeys.GenMACCSKeys(mol)   
    except:
        print('Something went wrong computing MACCSKeys')
        return None
    return np.array(mkeys)

def compute_atom_pair_fp(mol, nBits=2048):
    """Function that obtains the atom pair Fingerprints of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        atom_pair_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits)
    except:
        print('Something went wrong computing Atom Pair fingerprints')
        return None
    return np.array(atom_pair_fp)

def compute_topological_torsion_fp(mol, nBits=2048):
    """Function that obtains the topological torsion fingerprints of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        tt_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol)
    except:
        print('Something went wrong computing Topological Torsion fingerprints')
        return None
    return np.array(tt_fp)

def compute_avalon_fp(mol, nBits=2048):
    """Function that obtains the Avalon fingerprints of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        av_fp = pyAvalonTools.GetAvalonFP(mol, nBits)
    except:
        print('Something went wrong computing Avalon fingerprints')
        return None
    return np.array(av_fp)

def compute_rdkit_fp(mol, maxPath=5, fpSize=2048):
    """Function that obtains the RDKit fingerprints of a molecule
    Input: RDKit molecule
    Output: Numpy array
    """
    try:
        rdkit_fp = AllChem.RDKFingerprint(mol, maxPath, fpSize)
    except:
        print('Something went wrong computing RDKit fingerprints')
        return None
    return np.array(rdkit_fp)

def compute_pubchem_fingerprints(cid):
    """Function that obtains the PubChem fingerprints of a molecule
    Input: molecules's CID
    Output: Numpy array
    """
    try:
        comp = pcp.Compound.from_cid(int(cid))
        fp_bin = bin(int(comp.fingerprint, 16))[2:]   
    except:
        print('Something went wrong computing Pubchem fingerprints')
        return None
    return np.array(list(fp_bin)).astype('int')

def compute_cactvs_fingerprints(cid):
    """Function that obtains the Cactvs fingerprints of a molecule
    Input: molecule's CID
    Output: Numpy array
    """
    try:
        comp = pcp.Compound.from_cid(int(cid))
        cactvs_fp_bin = bin(int(comp.fingerprint, 16))[2:]
    except:
        print('Something went wrong computing Cactvs fingerprints')
        return None
    return np.array(list(cactvs_fp_bin)).astype('int')

#get rdkit molecules

def add_molecule_column(all_drugs: pd.DataFrame, columns: Dict):
    base_column = columns['base']
    print(base_column)
    calculated_column = columns['calculated']
    print(calculated_column)
    pt.AddMoleculeColumnToFrame(frame=all_drugs, smilesCol=base_column, molCol=calculated_column)
    print(all_drugs.head(3))
    return all_drugs

#feature engineering input table
def drop_columns(all_drugs: pd.DataFrame, columns: Dict):
    all_drugs = all_drugs.drop(columns['columns_to_drop'], axis=1)
    return all_drugs

def label_encoder(all_drugs: pd.DataFrame):
    le = preprocessing.LabelEncoder()
    all_drugs['Label'] = le.fit_transform(all_drugs['MATC_Code_Short'])
    return all_drugs

#dataframe construction

def get_fingerprints(all_drugs: pd.DataFrame) -> pd.DataFrame:
    
    all_drugs['FeatInvariants'] = all_drugs['Molecule'].map(compute_feature_invariants)
    all_drugs['ConnInvariants'] = all_drugs['Molecule'].map(compute_connectivity_invariants)
    all_drugs['Morgan2FP'] = all_drugs['Molecule'].map(compute_morgan_fp)
    all_drugs['MACCSKeys'] = all_drugs['Molecule'].map(compute_maccskeys)
    all_drugs['AtomPairFP'] = all_drugs['Molecule'].map(compute_atom_pair_fp)
    all_drugs['TopTorFP'] = all_drugs['Molecule'].map(compute_topological_torsion_fp)
    all_drugs['AvalonFP'] = all_drugs['Molecule'].map(compute_avalon_fp)
    
    all_drugs['PubchemFP']= all_drugs['CID'].map(compute_pubchem_fingerprints) #This takes over 1 hour in my computer
    #all_drugs['CactvsFP']= all_drugs['CID'].map(compute_cactvs_fingerprints) #This takes over 1 hour in my computer
    #all_drugs['RDKitFP']= all_drugs['Molecule'].map(compute_rdkit_fp) #This takes so long that crashes my computer, but I coudn't find a way around
    
    return all_drugs

def clean_dataset(all_drugs, selected_columns: Dict):
    all_drugs = drop_columns(all_drugs, selected_columns)
    all_drugs = label_encoder(all_drugs)
    return all_drugs

def get_model_input(all_drugs: pd.DataFrame, columns: Dict, selected_columns: Dict):
    all_drugs = add_molecule_column(all_drugs, columns)
    all_drugs = get_fingerprints(all_drugs)
    all_drugs = clean_dataset(all_drugs, selected_columns)
    return all_drugs
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp



query=st.text_input("Insert molecule's CID or SMILES: ")

def get_molecule_from_smiles(query):
    """
    Function that obtains an RDKit molecule object from a SMILES.
    Input: SMILES as a string
    Output: RDKit molecule
    """
    mol = Chem.MolFromSmiles(query)
    st.write(mol)
    return mol

get_molecule_from_smiles(query)
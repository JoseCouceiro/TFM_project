import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp



st.title('This is how you print a molecule')

cid=st.text_input('a cid, please: ')

if cid:
    comp_mol = pcp.Compound.from_cid(cid)
    smiles_mol = comp_mol.canonical_smiles
    mol = Chem.MolFromSmiles(smiles_mol)

    drawing = Draw.MolToImage(mol)


    st.image(drawing)

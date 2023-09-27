import streamlit as st
import os
import numpy as np
import json
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from tensorflow.keras.models import load_model
from compute_fp import Compute_FP

class Display:
    def show_display(self):
        st.title('Drug predictor')
        st.markdown('A drug discovery tool by José R. Couceiro')
        with st.sidebar:
            st.write("Get your molecule's CID [here](https://pubchem.ncbi.nlm.nih.gov/)")
            st.image(os.path.join('res', 'images', 'pubchem.png'), width=200)
            st.write("Draw your molecule and get its SMILES [here](https://web.chemdoodle.com/demos/smiles#customise-template)")
            st.image(os.path.join('res', 'images', 'chemdoodleweb.png'), width=200)

class Calcs:

    def __init__(self):

        self.computer = Compute_FP()
        
        self.CNN_model = load_model(os.path.join('..', '..', 'data', '06_models', 'def_model.hd5'))

        with open(os.path.join('..', '..', 'data', '03_primary', 'code_to_label_dic.json'), 'r') as file:
            self.class_codes_dict = json.load(file)

        with open(os.path.join('..', '..', 'data', '05_model_input', 'selected_fp.txt')) as file:
            self.selected_fp = file.readline()

    def get_molecule_from_cid(self, cid):
        """
        Function that obtains an RDKit molecule object from a Pubchem "CID" number.
        Input: Pubchem CID as an integer
        Output: RDKit molecule
        """
        comp_mol = pcp.Compound.from_cid(cid)
        smiles_mol = comp_mol.canonical_smiles
        mol = Chem.MolFromSmiles(smiles_mol)
        return mol

    def get_molecule_from_smiles(self, smiles):
        """
        Function that obtains an RDKit molecule object from a SMILES.
        Input: SMILES as a string
        Output: RDKit molecule
        """
        mol = Chem.MolFromSmiles(smiles)
        return mol

    def get_fp(self, mol):
        """
        Function that obtains the fingerprints from an RDKit molecule and returns them as a Numpy array.
        Input: RDKit Molecule
        Output: Numpy array
        """
        mc_fp = self.computer.relate_fp_functions(self.selected_fp, mol)
        return mc_fp

    def reshape_fp_array(self, fp_chain):
        """
        Function that reshapes a Numpy array to pass it into a 1D CNN model.
        Input: Numpy array with shape (n,)
        Output: Numpy array with shape (1, n, 1)
        """
        reshaped = np.array(fp_chain.reshape(1, fp_chain.shape[0], 1))
        return reshaped

    def return_predictions(self, arr, query):
        """
        Function that takes an array and predicts its label using a CNN model.
        Input: Numpy array
        Output: Prints predictions of which class a molecule is assigned to and the probability with which it is assigned to that class.
        """
        y_pred = self.CNN_model.predict(arr)
        st.markdown(f"##### The compound with {self.is_cid_or_smiles(query)} '{query}' is predicted as {self.class_codes_dict[str(np.argmax(y_pred))]} with probability: {y_pred.max():.2f}**")

    def return_output(self, mol, query):
        """
        Function that returns all the output of the application, giving the structure and the prediction fort the query molecule. It uses the functions \
        'get_fp', 'reshape_fp_array' and 'return_predictions' to achieve this.
        Input: molecule RDKit object
        Output: draws the molecule structure and prints predictions.
        """
        st.write("This is your compound's structure")
        mol_image = Draw.MolToImage(mol)
        st.image(mol_image)
        fp_chain = self.get_fp(mol)
        arr = self.reshape_fp_array(fp_chain)
        self.return_predictions(arr, query)

    def pred_from_cid(self, query):
        """
        Checks the input and runs all the functions in order to obtain a prediction from a CID identifier.
        Input: user input
        Output: prints prediction
        """
        mol = False
        while not mol:
            try:
                mol = self.get_molecule_from_cid(query)
                if mol:
                    self.return_output(mol, query)
            except:
                st.write(':red[Please, enter a valid CID]')
                break

    def pred_from_smiles(self, query):
        """
        Checks the input and runs all the functions in order to obtain a prediction from a CID identifier.
        Input: user input
        Output: prints prediction
        """
        mol = False
        while not mol:
            mol = self.get_molecule_from_smiles(query)
            if mol:
                self.return_output(mol, query)
            else:
                st.write(':red[Please, enter a valid SMILES]')
                break
            
    def is_cid_or_smiles(self, user_input):
        """
        This functions return what kind of input the user has entered
        Input: user input
        Output: a string
        """
        if user_input == '':
            return None
        try:
            if int(user_input):
                return 'CID'
        except:
            return('SMILES')
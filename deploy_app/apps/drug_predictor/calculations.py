import streamlit as st
import os
import numpy as np
import json
import pubchempy as pcp
from tensorflow import keras
from keras.models import load_model
from keras.layers import TFSMLayer
from sklearn import metrics
from rdkit import Chem
from rdkit.Chem import Draw
from compute_fp import Compute_FP

class Display:
    """
    This class generates the layout of Drug Predictor's application using Streamlit functions. It comprises 3 functions:
        'show_tab1': creates tab1 for Drug Predictor.
        'show_tab2': creates tab2 for Drug Predictor High Throughput.
        'show_sidebar': creates the sidebar with pictures and links to resources.
    """
    def show_tab1(self):
        st.title('Drug predictor')
        st.markdown('A drug discovery tool by José R. Couceiro')
    
    def show_tab2(self):
        st.title('Drug predictor high throughput')
        st.markdown('A drug discovery tool by José R. Couceiro')

    def show_sidebar(self):
        with st.sidebar:
            st.write("Get your molecule's CID [here](https://pubchem.ncbi.nlm.nih.gov/)")
            pubchem_image_path = os.path.join('deploy_app', 'apps', 'drug_predictor', 'res', 'images', 'pubchem.png')
            if os.path.isfile(pubchem_image_path):
                try:
                    st.image(pubchem_image_path, width=200)
                except Exception as e:
                    st.error(f"Error opening the image file 'pubchem.png': {e}")
            else:
                st.error(f"File not found: {pubchem_image_path}")
            
            st.write("Draw your molecule and get its SMILES [here](https://web.chemdoodle.com/demos/smiles#customise-template)")
            chemdoodle_image_path = os.path.join('res', 'images', 'chemdoodleweb.png')
            if os.path.isfile(chemdoodle_image_path):
                try:
                    st.image(chemdoodle_image_path, width=200)
                except Exception as e:
                    st.error(f"Error opening the image file 'chemdoodleweb.png': {e}")
            else:
                st.error(f"File not found: {chemdoodle_image_path}")
    
class Calcs:
    """
    This class loads the necessary files for the app Drug Predictor to run. See the documentation of each function for more details.
    """

    def __init__(self):

        self.computer = Compute_FP()
        
        self.CNN_model = load_model(os.path.join('deploy_app', 'data', '06_models', 'def_model.h5'))

        with open(os.path.join('deploy_app', 'data', '03_primary', 'code_to_label_dic.json'), 'r') as file:
            self.class_codes_dict = json.load(file)

        with open(os.path.join('deploy_app', 'data', '05_model_input', 'selected_fp.txt')) as file:
            self.selected_fp = file.readline()

    def get_molecule_from_cid(self, cid):
        """
        Function that obtains an RDKit molecule object from a Pubchem "CID" number.
        Input: Pubchem CID as an integer.
        Output: RDKit molecule.
        """
        try:
            comp_mol = pcp.Compound.from_cid(cid)
        except:
            st.error(f'Wrong cid format! {cid} is not a valid cid')

        smiles_mol = comp_mol.canonical_smiles
        mol = Chem.MolFromSmiles(smiles_mol)
        return mol
        

    def get_molecule_from_smiles(self, smiles):
        """
        Function that obtains an RDKit molecule object from a SMILES.
        Input: SMILES as a string.
        Output: RDKit molecule.
        """
        #try:
        mol = Chem.MolFromSmiles(smiles)
        #except:
        if not mol:
            st.error(f'Wrong SMILES! {smiles} is not a valid SMILES')
        return mol

    def get_fp(self, mol):
        """
        Function that obtains the fingerprints from an RDKit molecule and returns them as a Numpy array.
        Input: RDKit Molecule.
        Output: Numpy array.
        """
        mc_fp = self.computer.relate_fp_functions(self.selected_fp, mol)
        return mc_fp

    def reshape_fp_array(self, fp_chain, entry_type):
        """
        Function that reshapes a Numpy array to pass it into a 1D CNN model.
        Input: Numpy array with shape (n,) or (m,n,), entry_type either 'list' or 'single'.
        Output: Numpy array with shape (1, n, 1) or (m,n,1).
        """
        
        if entry_type == 'list':
            reshaped = np.array(fp_chain.reshape((fp_chain.shape[0], fp_chain.shape[1], 1)))
        else:
            reshaped = np.array(fp_chain.reshape(1, fp_chain.shape[0], 1))
        return reshaped
    
    def return_predictions_dataframe(self, df):
        """
        Function that takes a dataframe that contains a column with CID numbers or SMILES descriptors as values, trasforms those values to \
        feed them to a CNN model, and then obtains predictions and the probabilities of those predictions. It uses several functions to \
        achieve this.
        Args: a pandas dataframe.
        Output: two lists, one of predictions and one of probabilities.
        """
        try:
            if 'smiles' in df.columns:
                df['molecule'] = df['smiles'].map(self.get_molecule_from_smiles)
            elif 'cid' in df.columns:
                df['molecule'] = df['cid'].map(self.get_molecule_from_cid)
            else:
                raise KeyError             
            df['fingerprints'] = df['molecule'].map(self.get_fp)
            fingerprints = np.array(list(df['fingerprints']))
            reshaped_fps = self.reshape_fp_array(fingerprints, 'list')
            preds, probs = self.return_predictions(reshaped_fps, df, 'list')
            return preds, probs
        except KeyError:
            st.error(':red[Wrong format! Please, make sure there is at least one header called either "cid" or "smiles"]')
        

    def return_predictions(self, arr, query, query_type):
        """
        Function that takes an array and predicts its label using a CNN model.
        Input: Numpy array, a cid number or array of cids, the type of query (either 'list' or 'single']
        Output: Prints predictions of which class a molecule is assigned to and the probability with which it is assigned to that class.
        """
        if query_type == 'list':
            probs = self.CNN_model.predict(arr)
            preds = [np.argmax(x) for x in probs]
            max_probs = [np.max(x) for x in probs]
            return preds, max_probs
            
        else:
            y_pred = self.CNN_model.predict(arr)
            st.markdown(f"##### The compound with \
                {self.is_cid_or_smiles(query)} '{query}' is predicted as \
                {self.class_codes_dict[str(np.argmax(y_pred))]} with probability: \
                {y_pred.max():.2f}"
                )
    

    def return_output_dataframe(self, preds, probs, df):
        """
        Function that adds the columns 'prediction', 'probability', 'description', and 'url' to a dataframe.\
        Values for 'prediction' and 'probability' are inputted. 'description' is mapped from a dictionary, \
        and 'url' is calculated from the column 'cid', if there is one, otherwise a 'cid' is calculated from \
        a 'smiles' column.
        Args: one list of predictions and one of probabilities, a pandas dataframe.
        Output: none. Output data is handled by streamlit.
        """
        try:
            df['prediction'] = preds
            df['probability'] = probs
            df['description'] = df['prediction'].astype('str').map(self.class_codes_dict)
            if 'cid' not in df.columns:
                df['cid'] = df['smiles'].apply(lambda x: pcp.get_compounds(x, 'smiles')[0].cid)
                df['cid'] = df['cid'].replace({np.nan: None})
                df['cid'] = df['cid'].apply(lambda x: int(x) if x!=None else 0)
            df['url'] = df['cid'].apply(lambda x: f'https://pubchem.ncbi.nlm.nih.gov/compound/{x}' if x !=0 else None)
                
            output_df = df[['cid', 'prediction', 'probability', 'description', 'url']]
            st.dataframe(output_df)
            st.download_button(label='Download as a CSV file', data=output_df.to_csv(), file_name='drug_prediction.csv')
        except:
            st.error('Something went wrong while mounting the output dataframe :(')

    def return_output(self, mol, query):
        """
        Function that returns all the output of the application Drug Predictor: the structure and the prediction for the query molecule. \
              It uses the functions 'get_fp', 'reshape_fp_array' and 'return_predictions' to achieve this.
        Input: molecule RDKit object, a query (a CID number or SMILES string).
        Output: draws the molecule structure and prints predictions.
        """
        st.write("This is your compound's structure")
        mol_image = Draw.MolToImage(mol)
        st.image(mol_image)
        fp_chain = self.get_fp(mol)
        arr = self.reshape_fp_array(fp_chain, 'single')
        self.return_predictions(arr, query, 'single')

    def pred_from_cid(self, query):
        """
        Checks a query for its validity. If the query is valid, it runs it into the function 'return_output'.
        Input: user input.
        Output: none, it generates an RDKit molecule object that passes to 'return_output'.
        """
        mol = False
        while not mol:
            try:
                mol = self.get_molecule_from_cid(query)
                if mol:
                    self.return_output(mol, query)
            except:
                st.error(':red[Please, enter a valid CID]')
                break

    def pred_from_smiles(self, query):
        """
        Checks a query for its validity. If the query is valid, it runs it into the function 'return_output'.
        Input: user input.
        Output: none, it generates an RDKit molecule object that passes to 'return_output'.
        """
        mol = False
        while not mol:
            mol = self.get_molecule_from_smiles(query)
            if mol:
                self.return_output(mol, query)
            else:
                st.error(':red[Please, enter a valid SMILES]')
                break
    
    def save_evaluation_dataframe(self, df):
        """
        This function takes a dataframe and saves a column named 'prediction' as a pickle and a column named 'label as a csv.
        It also calculates a classification_report with these values and saves it as text.
        Input: a pandas dataframe.
        Output: saves to disc a pickle, a csv file and a text file.
        """
        df['prediction'].to_pickle(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_predictions.pickle'))
        df['label'].to_csv(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_true_values.csv'))

        class_report = metrics.classification_report(df['label'], df['prediction'])
        with open(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_class_report.txt'), 'w') as f:
            f.write(class_report)

            
    def is_cid_or_smiles(self, user_input):
        """
        This functions return what kind of input the user has entered.
        Input: user input.
        Output: a string.
        """
        if user_input == '':
            return None
        try:
            if int(user_input):
                return 'CID'
        except:
            return('SMILES')
        
        
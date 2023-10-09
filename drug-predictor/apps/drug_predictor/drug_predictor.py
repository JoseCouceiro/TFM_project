import os
import streamlit as st
from calculations import Display, Calcs
import pandas as pd


# Instantiate classes
displayer = Display()
calculator = Calcs()

tab1, tab2 = st.tabs(['Drug predictor', 'Drug predictor high throughput'])
displayer.show_sidebar()

with tab1:
    displayer.show_tab1()

    # Introduce input
    query = st.text_input("Insert molecule's CID or SMILES: ")

    # Calculate and return predictions
    calculator = Calcs()
    query_type = calculator.is_cid_or_smiles(query)
    if query_type == 'SMILES':
        calculator.pred_from_smiles(query)
    if query_type == 'CID':
        calculator.pred_from_cid(query)
    
with tab2:
    displayer.show_tab2()

    # Introduce input
    query = st.file_uploader("Please upload a CSV file. One molecule CID or SMILES per line. \
                             'cid' or 'smiles' as header for molecules. 'label' as header for labels (optional).")

    # Calculate and return predictions
    if query is not None:
        cids = pd.read_csv(query)     
        st.write(f'Analysing {len(cids)} molecules')
        try:
            preds, probs = calculator.return_predictions_dataframe(cids)
            calculator.return_output_dataframe(preds, probs, cids)
        except:
            st.image(os.path.join('..','..', 'apps', 'app', 'res', 'images', 'monkey_mistake.jpg'))

        if 'label' in cids.columns:
            calculator.save_evaluation_dataframe(cids)


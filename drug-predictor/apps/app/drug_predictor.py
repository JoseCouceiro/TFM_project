import streamlit as st
from calculations import Display, Calcs

# Show web display
displayer = Display()
displayer.show_display()

# Introduce input
query = st.text_input("Insert molecule's CID or SMILES: ")

# Calculate and return predictions
calculator = Calcs()
query_type = calculator.is_cid_or_smiles(query)
if query_type == 'SMILES':
    calculator.pred_from_smiles(query)
if query_type == 'CID':
    calculator.pred_from_cid(query)
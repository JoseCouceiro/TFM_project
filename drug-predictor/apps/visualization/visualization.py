import streamlit as st
from motor import Display

st.set_option('deprecation.showPyplotGlobalUse', False)
displayer = Display()
displayer.show_display()
import streamlit as st
import plotly.express as px
import yaml
from kedro.framework.session import KedroSession

#session = KedroSession('build_model')

context = KedroSession.load_context('./')
data = context.catalog.load('gitter')

found_max_trials = int(data['tune_model']['max_trials'])

st.title('Kedro model')

max_trials = st.slider('Maximum number of trials')


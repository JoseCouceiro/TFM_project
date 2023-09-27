import os
import pandas as pd
import streamlit as st
from draw_graphs import ClassificationReportGraph

class_graph = ClassificationReportGraph()
with open(os.path.join('..', '..', 'data', '08_reporting', 'classiffication_report.txt')) as file:
    clas_rep = file.read()

lines = clas_rep.split('\n')
accuracy_score = float(lines[-4:-3][0][-15:-5].strip())

st.write('Accuracy Score: ', accuracy_score)

fig = class_graph.plot_classification_report(clas_rep)
st.set_option('deprecation.showPyplotGlobalUse', False)
st.pyplot(fig)
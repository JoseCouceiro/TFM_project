import os
import pandas as pd
import streamlit as st
from draw_graphs import Graphs

graph = Graphs()
with open(os.path.join('..', '..', 'data', '08_reporting', 'classiffication_report.txt')) as file:
    train_clas_rep = file.read()
with open(os.path.join('..', '..', 'data', '08_reporting', 'clas_rep_val.txt')) as file:
    val_clas_rep = file.read()
with open(os.path.join('..', '..', 'data', '08_reporting', 'clas_rep_stage3.txt')) as file:
    stage3_clas_rep = file.read()
train_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'evaluation.pickle'))
train_splits = pd.read_pickle(os.path.join('..', '..', 'data', '05_model_input', 'split_col.pickle'))
val_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'evaluation.pickle'))
val_splits = pd.read_pickle(os.path.join('..', '..', 'data', '05_model_input', 'split_col.pickle'))
stage3_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'evaluation.pickle'))
stage3_splits = pd.read_pickle(os.path.join('..', '..', 'data', '05_model_input', 'split_col.pickle'))


def get_acc_score(clas_rep):
    lines = clas_rep.split('\n')
    return float(lines[-4:-3][0][-15:-5].strip())

class_rep_dic = {'training data': [train_clas_rep, train_preds, train_splits],
                 'validation data': [val_clas_rep, val_preds, val_splits],
                 'stage 3 examples': [stage3_clas_rep, stage3_preds, stage3_splits]
                                      }

st.title('Evaluation of drug_predictor model')

selected_clas_rep = st.selectbox('Please, choose the dataset',
                list(class_rep_dic.keys()))

accuracy_score = get_acc_score(class_rep_dic[selected_clas_rep][0])

st.write('Accuracy Score: ', accuracy_score)

st.write('Classification report')
fig1 = graph.plot_classification_report(class_rep_dic[selected_clas_rep][0])
st.pyplot(fig1)

st.write('Confusion matrix')
fig2 = graph.plot_confusion_matrix(class_rep_dic[selected_clas_rep][1], class_rep_dic[selected_clas_rep][2][3])
st.pyplot(fig2)

st.set_option('deprecation.showPyplotGlobalUse', False)
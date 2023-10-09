import os
import pandas as pd
import streamlit as st
from draw_graphs import Graphs

class Loads():
    def load_class_reports(self):
        with open(os.path.join('..', '..', 'data', '08_reporting', 'train_classification_report.txt')) as file:
                train_clas_rep = file.read()
        with open(os.path.join('..', '..', 'data', '08_reporting', 'validation_classification_report.txt')) as file:
                val_clas_rep = file.read()
        try:
            with open(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_class_report.txt')) as file:
                ht_clas_rep = file.read()    
        except:
            ht_clas_rep = None
        return train_clas_rep, val_clas_rep, ht_clas_rep

    def load_predictions(self):
        train_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'train_predictions.pickle'))
        val_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'validation_predictions.pickle'))
        try:
            drug_predictor_ht_preds = pd.read_pickle(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_predictions.pickle'))
        except:
            drug_predictor_ht_preds = None
        return train_preds, val_preds, drug_predictor_ht_preds

    def load_true_values(self):
        train_true_values = pd.read_pickle(os.path.join('..', '..', 'data', '05_model_input', 'split_col.pickle'))
        val_true_values = pd.read_pickle(os.path.join('..', '..', 'data', '05_model_input', 'input_table_val.pickle'))
        try:
            drug_predictor_ht_true_values = pd.read_csv(os.path.join('..', '..', 'data', '08_reporting', 'drug_predictor_ht_true_values.csv'))
        except:
            return train_true_values[3], val_true_values['Label'], None
        return train_true_values[3], val_true_values['Label'], drug_predictor_ht_true_values['label']
    
class Organizer():
    def __init__(self):
        self.loader = Loads()
        class_reports = self.loader.load_class_reports()
        predictions = self.loader.load_predictions()
        true_values = self.loader.load_true_values()
        self.class_rep_dic = {
                        'training data': [class_reports[0], predictions[0], true_values[0]],
                        'validation data': [class_reports[1], predictions[1], true_values[1]],
                        'drug predictor high throughput': [class_reports[2], predictions[2], true_values[2]]
                        }

    def get_acc_score(self, clas_rep):
        lines = clas_rep.split('\n')
        return float(lines[-4:-3][0][-15:-5].strip())


class Display():

    def __init__(self):
        self.organizer = Organizer()
        self.graph = Graphs()

    def show_display(self):
        tab1, tab2 = st.tabs(['evaluate predictions', 'evaluate model training'])
    
        with tab1:
            st.title('Evaluation of drug_predictor model')
            selected_clas_rep = st.selectbox('Please, choose a dataset',
                            list(self.organizer.class_rep_dic.keys()))
            try:
                accuracy_score = self.organizer.get_acc_score(self.organizer.class_rep_dic[selected_clas_rep][0])
            except AttributeError:
                st.error('You must run some predrictions in Drug Predictor High Throuput to see the analysis')
            else:     
                st.write('**:orange[Accuracy Score: ]**', f'**:blue[{accuracy_score}]**')
                st.markdown('**:orange[Classification report]**')
                try:
                    fig1 = self.graph.plot_classification_report(self.organizer.class_rep_dic[selected_clas_rep][0])
                    st.pyplot(fig1)

                    st.markdown('**:orange[Confusion matrix]**')
                
                    fig2 = self.graph.plot_confusion_matrix(self.organizer.class_rep_dic[selected_clas_rep][1],
                                                            self.organizer.class_rep_dic[selected_clas_rep][2])
                    st.pyplot(fig2)
                except KeyError:
                    st.error('A label was not recognized by the model! Please, make sure that all your labels were used to train the model')

        with tab2:
            st.title("Model's training history")
            st.image(os.path.join('..', '..', 'data', '08_reporting', 'training_fig.png'))


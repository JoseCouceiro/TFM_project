import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from tensorflow import keras
from tensorflow.keras.models import Sequential

def get_predictions(model: Sequential, x_test: pd.Series, y_test: pd.Series):
    """
    Function that obtains the predictions of a CNN model and a classification report.
    Args: 
      model: a CNN model.
      y_train: pandas series containing y_train data.
      y_test: pandas series containing y_test.
    Output: predictions, classification report.
    """
    y_pred = model.predict(x_test)
    y_pred_list = [np.argmax(x) for x in y_pred]
    class_rep = metrics.classification_report(y_test,y_pred_list)
    print(class_rep)
    return pd.Series(y_pred_list), class_rep

def evaluate_model(model: Sequential, model_data:tuple, split_col:tuple, val_dataset: pd.DataFrame, selected_fp: str):
    """
    Function that takes the data of x_test from both a training dataset and a validation dataset, and\
    makes the predictions using CNN model. The function also takes the true values (y_test) of both datasets and uses them \
    to make a classification report. The function returns both predictions and classification reports.
    Args: model: CNN model.
          model_data: tuple containing x_test.
          split_col: tuple containing y_test.
          val_dataset: a pandas dataframe with the validation data.
          selected_fp: a string indicating the selected fingerprints.
    Output: tuple containing predictions and classification reports for both training data and validation data.
    """
    # Reshaping input data
    fingerprints = np.array(list(val_dataset[selected_fp]))
    reshaped_fps = fingerprints.reshape((fingerprints.shape[0], fingerprints.shape[1], 1))
    
    dic = {'training': {'x_test': model_data[1], 'true_value': split_col[3]},
           'validation': {'x_test': reshaped_fps, 'true_value': val_dataset['Label']}}
    
    val_predictions, val_class_rep = get_predictions(model, dic['validation']['x_test'], dic['validation']['true_value'])
    train_predictions, train_class_rep = get_predictions(model, dic['training']['x_test'], dic['training']['true_value'])
    return train_predictions, train_class_rep, val_predictions, val_class_rep

def visualize_training(history_dic):
    """
    Function that plots the training history of a CNN model showing values of training and validation loss and accuracy for each epoch.
    Args: a training history.
    Output: a matplotlib figure.
    """
    fig,ax = plt.subplots()
    plt.style.use('ggplot')
    epochs = len(history_dic['loss'])
    epoch_values = list(range(epochs))

    ax.plot(epoch_values, history_dic['loss'], label='Training loss')
    ax.plot(epoch_values, history_dic['val_loss'], label='Validation loss')
    ax.plot(epoch_values, history_dic['accuracy'], label='Training accuracy')
    ax.plot(epoch_values, history_dic['val_accuracy'], label='Validation accuracy')

    ax.set_title('Training loss and accuracy')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Loss/Accuracy')
    ax.legend()
    return fig
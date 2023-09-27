import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import metrics
from numpy import unique, argmax
from keras_tuner import RandomSearch
from tensorflow.keras import layers
from tensorflow import keras
from tensorflow.keras.models import Sequential #load_model
from tensorflow.keras.layers import Dense, Conv1D, MaxPool1D, Flatten, Dropout
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

# Functions to prepare data input

def train_test_split_column(input_pickle: pd.DataFrame, column: str, label: dict) -> tuple:
    """
    Function that makes a train/test split of a selected column of a pandas dataframe. Y values are taken from the label column.
    Args:
      input_pickle: input dataframe.
      column: name of the column selected as X values.
      label: name of the column selected as y values.
    Output: train/test tuple.
    """
    print(len(input_pickle.shape))
    X = input_pickle[column]
    y = input_pickle[label]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    return X_train, X_test, y_train, y_test

def train_test_split_columns(input_pickle: pd.DataFrame, column_list: dict, label: dict) -> dict:
    """
    Function that takes a list of columns and makes a train/test split for each of them with the column containing the labels as y.
    Then returns a dictionary with the column name as a key and the tuple with X_train, X_test, y_train, y_test as value.
    Args:
      input_pickle: input pandas dataframe.
      columns_list: list of column names.
      label: name of the column selected as y values.
    Output: dictionary.
    """
    splits_dic = {}
    for column in column_list['list_of_fingerprints']:
        X_train, X_test, y_train, y_test = train_test_split_column(input_pickle, column, label['label'])
        splits_dic[column] = X_train, X_test, y_train, y_test
        print('X_train_shape: ', X_train.shape)
    return splits_dic
    

def reshape_input(tuple_split: tuple) -> tuple:
    """
    Function that reshapes the input data to fit the CNN model and determines both the number of classes from the label data \
        and the shape of the reshaped data.
    Input: train/test split tuple.
    Output: tuple with the data necessary to run the CNN model.
    """
    # Transform arrays into lists
    X_train = np.array(list(tuple_split[0]))
    X_test = np.array(list(tuple_split[1]))
    # Determine number of classes
    n_classes = len(np.unique(tuple_split[2]))
    # Reshape X
    reshaped_X_train= X_train.reshape((X_train.shape[0], X_train.shape[1], 1))
    reshaped_X_test= X_test.reshape((X_test.shape[0], X_test.shape[1], 1))
    # Determine in_shape
    in_shape = reshaped_X_train.shape[1:]
    return reshaped_X_train, reshaped_X_test, in_shape, n_classes

# Functions to build a model to select the best fingerprints

def build_array_dic(splits_dic: dict) -> dict:
    """
    Function that takes the dictionary generated by "train_test_split_columns" and reshapes X_train to fit a CNN model. \
         Then builds the CNN model using the function 'build_selection_model'. Returns both the reshaped arrays and the models.
    Input: dictionary with column names as keys and a tuple of arrays as values.
    Output: dictionary with column names as keys and a tuple with both reshaped arrays and CNN models.
    """
    arrays_models_dic = {}
    for column, tup in splits_dic.items():                                                                                                                                                                                                                                                                              
        model_data = reshape_input(tup)  
        model = build_selection_model(model_data[2], model_data[3])
        arrays_models_dic[column] = model_data[0], model_data[1], tup[2], tup[3], model
    return arrays_models_dic

def build_selection_model(inshape: tuple, nclasses: int) -> Sequential:
    """
    Function that builds a CNN model with non optimal conditions.
    Input: tuple with inner shape of arrays, integer with number of classes.
    Output: a CNN model.
    """
    model = Sequential()
    model.add(Conv1D(100, 9, activation='relu', kernel_initializer='he_uniform', input_shape=inshape))
    model.add(MaxPool1D(2))
    model.add(Flatten())
    model.add(Dense(100, activation='relu', kernel_initializer='he_uniform'))
    model.add(Dropout(0.5))
    model.add(Dense(nclasses, activation='softmax'))
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

def fit_selection_model(arrays_models_dic: dict) -> dict:
    """
    Function that takes the dictionary generated by 'reshape_and_build' and fits the CNN model contained in the values \
        with the X_train and y_train arrays also in the values. Then evaluates the model against X_test e y_test \
        and returns the accuracy of each model.
    Input: dictionary with column names as keys and a tuple with reshaped arrays and CNN models as values.
    Output: a dictionary with column names as keys and the accuracy obtained by the CNN model for each column as values.
    """
    accuracies_dic = {}
    es = EarlyStopping(monitor='val_loss', patience=1)
    for column, tup in arrays_models_dic.items():
        print(f"Analysing {column}")
        tup[4].fit(tup[0], tup[2], epochs=10, batch_size=128, verbose=1, validation_split = 0.2, callbacks = [es])
        loss, acc = tup[4].evaluate(tup[1], tup[3], verbose=1)
        accuracies_dic[column] = f'{acc:.3f}'
    return accuracies_dic

# Functions to build the definitive model

def tune_hp_def_model(model_data: tuple, y_train: pd.Series, y_test: pd.Series, tune_params: dict) -> Sequential:
    """
    Funtion that uses RandomSearch from keras.tuner to obtain the best hyperparameters for a CNN model.
    Args:
      model_data: tuple containing pandas series for X_train, X_test and the values for the in_shape and number of classes.
      y_train: pandas series with y_train data.
      y_test: pandas series with y_test.
      tune_params: a set of optional parameters.
    Output: a tuned CNN model.
    """
    print('y_train type: ', type(y_train))
    print('X_test type: ', type(y_test))
    tuner = RandomSearch(build_def_model,
                    objective=tune_params['objective'],
                    max_trials =tune_params['max_trials'],
                    directory = os.path.join('temp', 'tuner', tune_params['date']))
    tuner.search(model_data[0],y_train,epochs=3,validation_data=(model_data[1],y_test))
    tuned_model = tuner.get_best_models(num_models=1)[0]
    tuned_model.summary() 
    return tuned_model

def build_def_model(hp) -> Sequential:
    """
    Function that builds a CNN model with sets of hyperparameters that will be chosen by the function "tune_hp_def_model".
    Input: a set of hyperparamneters
    Output: a cnn model
    """
    # Create model object
    model = keras.Sequential()
    # Choose number of layers
    for i in range(hp.Int("num_layers", 1, 5)):
        model.add(
            layers.Conv1D(
        filters=hp.Int('conv_1_filter', min_value=16, max_value=128, step=16),
        kernel_size=hp.Choice('conv_1_kernel', values = [3,5]),
        activation='relu',
        input_shape=(2048, 1),
        padding='valid')) #no padding
        model.add(layers.MaxPool1D(hp.Int('pool_size', min_value=2, max_value=6)))
        if hp.Boolean("dropout"):
            model.add(layers.Dropout(rate=0.25))
    model.add(layers.Flatten())
    model.add(layers.Dense(
        units=hp.Int('dense_1_units', min_value=32, max_value=128, step=16),
        activation='relu', kernel_initializer = 'he_uniform'
    ))
    model.add(layers.Dropout(0.5))
    model.add(layers.Dense(16, activation='softmax'))
    # Compilation of model
    model.compile(optimizer=keras.optimizers.Adam(hp.Choice('learning_rate', values=[1e-2, 1e-3])),
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
    return model

def train_def_model(tuned_model: Sequential, model_data: tuple, y_train: pd.Series, y_test: pd.Series) -> tuple:
    """
    Function that trains a CNN model.
    Args: 
      tuned_mode: a CNN model.
      model_data: tuple containing pandas series for X_train, X_test and the values for the in_shape and number of classes.
      y_train: pandas series with y_train data.
      y_test: pandas series with y_test.
    Output: a cnn model, the training history.
    """
    # Configure early stopping
    es = EarlyStopping(monitor='val_loss', patience=10)
    mc = ModelCheckpoint(filepath = os.path.join('temp', 'compiled_model', 'checkpoints', '{epoch:02d}-{val_accuracy:.3f}.hdf5'), monitor = 'val_loss', save_best_only = True)
    # Fit the model
    history = tuned_model.fit(model_data[0], y_train, epochs=200, batch_size=128, verbose=1, validation_split = 0.3, callbacks = [es,mc])
    # Evaluate the model
    loss, acc = tuned_model.evaluate(model_data[1], y_test, verbose=1)
    print('history type: ', type(history))
    
    return tuned_model, history

def get_predictions(tuned_model: Sequential, x_test: pd.Series, y_test: pd.Series):
    """
    Function that obtains the predictions of a CNN model and determines the accuracy.
    Args: 
      tuned_mode: a CNN model.
      model_data: tuple containing pandas series for X_train, X_test and the values for the in_shape and number of classes.
      y_train: pandas series with y_train data.
      y_test: pandas series with y_test.
    Output: a cnn model, the training history.
    """
    y_pred = tuned_model.predict(x_test)
    y_pred_list = [argmax(x) for x in y_pred]
    class_rep = metrics.classification_report(y_test,y_pred_list)
    print(class_rep)
    return pd.Series(y_pred_list), class_rep

def visualize_training(history):
    """
    Function that plots the training history of a CNN model showing values of training and validation loss and accuracy.
    Args: the training history
    Output: nothing
    """
    fig,ax = plt.subplots()
    plt.style.use('ggplot')

    epochs = len(history.history['loss'])
    epoch_values = list(range(epochs))

    ax.plot(epoch_values, history.history['loss'], label='Training loss')
    ax.plot(epoch_values, history.history['val_loss'], label='Validation loss')
    ax.plot(epoch_values, history.history['accuracy'], label='Training accuracy')
    ax.plot(epoch_values, history.history['val_accuracy'], label='Validation accuracy')

    ax.set_title('Training loss and accuracy')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Loss/Accuracy')
    ax.legend()
    plt.show()

# Combined functions

def determine_best_fp(input_pickle: pd.DataFrame, column_list: list, label: dict) -> dict:
    """
    Function that splits and reshape the data of several columns (determined by input) so a CNN model can be build for each of these columns \
    (each corresponding to a different type of fingerprint). The accuracy of each model is calculated and returned in a dictionary.
    Args: 
      input_pickle: input pandas dataframe
      column_list: list of column names (determined in parameters).
      label: name of the column selected as y values (determined in parameters).
    Output: a dictionary with column names as keys and the accuracy obtained by the CNN model for each column as values.
    """
    splits_dic = train_test_split_columns(input_pickle, column_list, label)
    array_model_dic= build_array_dic(splits_dic)
    accuracies_dic = fit_selection_model(array_model_dic)
    return accuracies_dic

def prepare_data(input_pickle: pd.DataFrame, fingerprints_accuracies: dict, label: dict) -> tuple:
    """
    Function that selects the fingerprints with wich the best accuracy is obtained, train/test splits the data of those fingerprints, \
        reshapes X_train and X_test to fit the model, and obtains the number of classes and the in_shape.
    Args: 
      input_pickle: input pandas dataframe
      fingerprints_accuracies: a dictionary with column names as keys and the accuracy obtained by the CNN model for each column as values.
      label: name of the column selected as y values (determined in parameters).
    Output:
      split_col: train/test tuple.
      model_data: tuple containing pandas series for X_train, X_test and the values for the in_shape and number of classes.
      selected_fp: the name of the fingerprints with a higher accuracy value as a string.
    """
    selected_fp = max(fingerprints_accuracies, key=fingerprints_accuracies.get)
    print("selected fingerprints: ", selected_fp)
    split_col = train_test_split_column(input_pickle, selected_fp, label['label'])
    model_data = reshape_input(split_col)
    return split_col, model_data, selected_fp

def obtain_model(split_col: tuple, model_data: tuple, tune_params: dict) -> tuple:
    """
    Function that builds a CNN model, chooses the best hyperparameters, trains it and return the tuned model and the predictions obtained during the training.
    Args: 
      split_col: train/test tuple.
      model_data: tuple containing pandas series for reshaped X_train, X_test and the values for the in_shape and number of classes.
    Output:
      def_model: a tuned and trained CNN_model.
      predictions: the predictions obtained in the training.
    """
    tuned_model = tune_hp_def_model(model_data, split_col[2], split_col[3], tune_params)
    def_model, history = train_def_model(tuned_model, model_data, split_col[2], split_col[3])
    predictions, class_rep = get_predictions(tuned_model, model_data[1], split_col[3])
    #visualize_training(history)
    return def_model, history, predictions, class_rep
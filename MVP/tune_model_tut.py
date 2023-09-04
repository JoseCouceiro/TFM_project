import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import metrics
from numpy import unique
from numpy import argmax
from keras_tuner import RandomSearch
from tensorflow.keras import layers
from tensorflow import keras
from tensorflow.keras.datasets.mnist import load_data
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Conv1D, MaxPool1D, Flatten, Dropout
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint


input_pickle = pd.read_pickle('C:/Users/josin/GitRepositories/TFM_project/drug-predictor/data/05_model_input/input_table.pickle/2023-08-25T16.54.10.335Z/input_table.pickle')
fingerprints_accuracies = {
  "Morgan2FP": "0.484",
  "MACCSKeys": "0.387",
  "AtomPairFP": "0.387",
  "TopTorFP": "0.403",
  "AvalonFP": "0.419",
  "PubchemFP": "0.355"
}

def train_test_split_column(input_pickle, column):
    X = input_pickle[column]
    y = input_pickle['Label'] #Label drug_class_code
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    return X_train, X_test, y_train, y_test

def reshape_input(tuple_split):
    X_train = np.array(list(tuple_split[0]))
    print('Shape X_train: ', X_train.shape, type(X_train))
    X_test = np.array(list(tuple_split[1]))
    n_classes = len(np.unique(tuple_split[2]))
    print('Number of classes: ', n_classes)
    X_train= X_train.reshape((X_train.shape[0], X_train.shape[1], 1))
    X_test= X_test.reshape((X_test.shape[0], X_test.shape[1], 1))
    print('Reshaped X_train: ', X_train.shape)
    in_shape = X_train.shape[1:]
    print('In_shape: ', in_shape)
    print('Shape y_train: ', tuple_split[2].shape)
    return X_train, X_test, in_shape, n_classes

def tune_hp_def_model(model_data, y_train, y_test):
    print('X_train shape: ', model_data[0].shape)
    print('y_train shape: ', y_train.shape)
    print('X_test shape: ', model_data[1].shape)
    print('X_test shape: ', y_test.shape)
    tuner = RandomSearch(build_def_model,
                    objective='val_loss',
                    max_trials = 5,
                    directory = os.path.join('tempHERE', 'tuner', 'RS_tuned_model', '230714_08'))
    tuner.search(model_data[0],y_train,epochs=3,validation_data=(model_data[1],y_test))
    tuned_model = tuner.get_best_models(num_models=1)[0]
    tuned_model.summary()
    return tuned_model

def train_def_model(tuned_model, model_data, y_train, y_test):
    # Configure early stopping
    es = EarlyStopping(monitor='val_loss', patience=10)
    mc = ModelCheckpoint(filepath = os.path.join('compiled_modelsHERE','checkpoints', '{epoch:02d}-{val_accuracy:.3f}.hdf5'), monitor = 'val_loss', save_best_only = True)
    # Fit the model
    history = tuned_model.fit(model_data[0], y_train, epochs=200, batch_size=128, verbose=1, validation_split = 0.3, callbacks = [es,mc])
    # Evaluate the model
    loss, acc = tuned_model.evaluate(model_data[1], y_test, verbose=1)
    return mc

def build_def_model(hp):
    """
    Function that chooses the best hyperparameters for a CNN model and then compiles it.
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


def obtain_model(input_pickle, fingerprints_accuracies):
    selected_fp = max(fingerprints_accuracies, key=fingerprints_accuracies.get)
    print('selected_fp: ', selected_fp)
    split_col = train_test_split_column(input_pickle, selected_fp)
    model_data = reshape_input(split_col)
    print('split_col: ', split_col)
    print('model_data[0] shape: ', model_data[0].shape)
    tuned_model = tune_hp_def_model(model_data, split_col[2], split_col[3])
    def_model = train_def_model(tuned_model, model_data, split_col[2], split_col[3])
    return def_model

obtain_model(input_pickle, fingerprints_accuracies)
# TFM_project José Rodríguez Couceiro
### Drug Predictor: Predicting Medicinal Application From Molecular Structure
#### A drug discovery tool

This repository contains the code to run the app Drug Predictor,
a python application capable of predicting the therapeutic potential of molecules from their chemical structures.

It serves as Master thesis for the Master's degree in Data Science 2022/23 at KSCHOOL together with the Memoria that can be found in the folder 'memoria'.

The project comprises a Kedro application as backend, which generates a convolutional neural network model capable of predicting the medicinal application of molecules from their molecular descriptors.

As frontend, two applications can be found in the folder \drug-predictor\apps:

1. drug_predictor: application where the final user can make queries for medicinal application of drugs using PubChem CIDs or SMILES codes.
![image](https://github.com/JoseCouceiro/TFM_project/assets/118387556/84bb072e-b7b4-4ec6-b367-d0c1e8681a17)


2. visualization: application for the user to analyse model training and performance.
![image](https://github.com/JoseCouceiro/TFM_project/assets/118387556/26ddd129-1009-4339-a6dd-9cb0a5592c83)


To run the project, the following steps are required:

-> Mount the environment:

A requirements text for a conda environmnet can be found at \drug-predictor\src\requirements.txt.

Create a new conda environment with the command:

    conda create --name <YOUR_ENV>> python==<YOUR_VERSION>

Run from inside the project's main folder:

    ~\anaconda3\envs\<YOUR_ENV>\python.exe  -m pip install -r src/requirements.txt

 -> Run the Kedro application:

The raw data can be found at https://1drv.ms/f/s!Aorqmaz_NWu2j9FQhoGSSgid9dQVfg?e=PQWhxx.

First, the raw data must be moved into the folder \drug-predictor\data\01_raw

(Many more details about the data are provided in the README_evaluator.md file inside \drug-predictor)

Then run:

    kedro run

-> This will generate a model whose performance can be analysed with the Visualization app.

To run it, move to \drug-predictor\apps\visualization, then run:

    streamlit run visualization.py

-> Finally, the user interface of the model, where a potential user can send queries to the model, can be found at \drug-predictor\apps\drug_predictor

To run it, use the command:

    streamlit run drug_predictor.py

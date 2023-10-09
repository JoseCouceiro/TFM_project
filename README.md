# TFM Project: José Rodríguez Couceiro
### Drug Predictor: Predicting Medicinal Application From Molecular Structure
#### A drug discovery tool

This repository contains the code to run the app Drug Predictor,
a python application capable of predicting the therapeutic potential of molecules from their chemical structures.

It serves as the Master thesis for the Master's degree in Data Science 2022/23 at KSCHOOL together with the Memoria that can be found in the folder 'memoria'.

To run the project, the following steps are required:
 -> Mount the environment:
     A requirements text for a conda environmnet can be found at \drug-predictor\src\requirements.txt.
     Activate conda and run from inside the project's folder:
          pip install -r src/requirements.txt
 -> Run the Kedro application:
     First, the data must be moved into the folder \drug-predictor\data\01_raw
     Many more details about the data are provided in the README_evaluator.md file inside \drug-predictor
     Then, inside the projects folder, run:
         kedro run
 -> This will generate a model whose performance can be analysed with the Visualization app.
      To run it move to \drug-predictor\apps\visualization, then use:
          streamlit run visualization.py
 -> Finally, the user interface of the model, where a potential user can send queries to the model, can be found at \drug-predictor\apps\drug_predictor
     To run it, use the command:
         streamlit run drug_predictor.py


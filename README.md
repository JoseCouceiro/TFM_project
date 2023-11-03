# TFM_project José Rodríguez Couceiro
### Drug Predictor: Predicting Medicinal Application From Molecular Structure
#### A drug discovery tool

This repository contains the code to run the app Drug Predictor,
a python application capable of predicting the therapeutic potential of molecules from their chemical structures.

It serves as Master thesis for the Master's degree in Data Science 2022/23 at KSCHOOL together with the Memoria that can be found in the folder 'memoria'.

### Project Structure

The folder 'drug-predictor' comprises a Kedro application organized in the following manner:

### Running the project

To run the project, the following steps are required:

-> Mount the environment:

    A requirements text for a conda environmnet can be found at \drug-predictor\src\requirements.txt.
    Create a new conda environment with the command:
        conda create --name <YOUR_ENV>> python==<YOUR_VERSION>
    and run from inside the project's main folder:
        ~\anaconda3\envs\<YOUR_ENV>\python.exe  -m pip install -r src/requirements.txt

 -> Run the Kedro application:

    The raw data can be found at https://1drv.ms/f/s!Aorqmaz_NWu2j9FQhoGSSgid9dQVfg?e=PQWhxx
    First, the raw data must be moved into the folder \drug-predictor\data\01_raw

    (Many more details about the data are provided in the README_evaluator.md file inside \drug-predictor)

    Then run:

        kedro run

-> This will generate a model whose performance can be analysed with the Visualization app.

    To run it move to \drug-predictor\apps\visualization, then run:

        streamlit run visualization.py

-> Finally, the user interface of the model, where a potential user can send queries to the model, can be found at \drug-predictor\apps\drug_predictor

    To run it, use the command:

        streamlit run drug_predictor.py

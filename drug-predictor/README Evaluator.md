DRUG PREDICTOR
Notes for the evaluator

The Drug Predictor package contains a Kedro application and two complementary apps contained in the folder apps

1. Input data

3 datasets are provided that are used by the kedro application to work correctly.
Since running times are in the range of hours when the aproximately 10000 molecules are analysed,
shorter versions of the datasets are provided (marked as mock_).
You can determine which ones to use in the Kedro catalog.

In the documentation acompanying this work, there are some assays made with variations of the datasets, 
a dataset with fewer labels and a dataset with a more balanced distributions of labels.
This datasets are available in the folder 03_primary and can be run from pipeline 'Molecular_calculations' onwards.

Folder data/09_examples contains small datasets that can be used with the application Drug Predictor High Throughput.
They contain several examples of use with possible exceptions and errors.

2. Complementary apps

They can be found in the apps folder and run through Streamlit
Esta carpeta contiene el trabajo final del Master en Data Science de KSCHOOL
Realizado por Jose Couceiro

El trabajo se compone de una serie de notebooks y una pequeña app.
El objetivo es determinar qué tipo de fármaco es una molécula a partir de sus determinantes estructurales of fingerprints.
Como material de partida contamos con una tabla en la que las moléculas vienen identificadas por códigos y etiquetadas en 12 tipos distintos de fármacos.
El trabajo está pensado para obtener información de forma progresiva siguiendo el siguiente orden:

1 - fp_df_construction.ipynb

Para construir un dataframe que contenga las etiquetas y junto con los fingerprints.
Carga los datos de origen
Obtiene un nuevo dataframe con los fingerprints y lo salva en formato pickle

2 - select_fingerprints.ipynb

Selecciona que fingerprints funcionan mejor en un modelo de redes neuronales convolucionales.
Carga el pickle con los fingerprints
Obtiene la información

3 - tune_model.ipyn

Construye un modelo con la mejor precisión para los datos de entrenamiento elegidos.
Carga el pickle con los fingerprints
Obtiene el modelo y lo salva

4 - lipinski.ipynb

Analiza la relevancia que el modelo le ha dado a las propiedades que tiene en cuenta la regla de Lipinski.
Carga el pickle con los fingerprints
Obtiene información

5 - Correr la app

drug_finder es un script the python que permite introducir el 'CID' (identificador de una molécula en la base de datos PubChem) y te devuelve la predicción sobre la clase de fármaco
Carga el modelo CNN
Devuelve una predicción

NOTAS PARA EL EVALUADOR
He incluido un archivo requirements.txt, creo que debido a mi impericia yo no soy capaz de cargarlo, pero creo que el archivo debería funcionar bien en manos más habiles.
En caso de encontrar algún problema, el código de instalación de las librerías no usadas en clase está incluido en los notebooks.




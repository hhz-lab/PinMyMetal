import pandas as pd
from joblib import dump, load
from keras.models import load_model
import numpy as np
import os, sys
import getopt
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils"))
import db_utils

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'chedh_metalclass_result.csv')
# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()

model_file = os.path.join(script_directory, 'ml_model_scripts', 'chedh_classmodel.joblib')

classmodel = load(model_file)

sql = """
select distinct * from classmodel_feature_neighbor_chedh
"""
cur.execute(sql)
data = cur.fetchall()
cur.close()
conn.close()

raw = pd.DataFrame(data)
raw.columns = [desc[0] for desc in cur.description] # add column name

raw.fillna(0, inplace=True)

X_test = raw.drop(columns=['id', 'sitetype']).values

# Predictions from metalclass models
# returns the probability that each sample belongs to each class
probabilities = classmodel.predict_proba(X_test)
# returns the model's final predicted class for each sample, which is the class with the highest probability.
predicted_classes = classmodel.predict(X_test)

max_probabilities = np.max(probabilities, axis=1)  #The maximum probability for each sample

output_df = raw[['id', 'sitetype']].copy()
output_df['result'] = predicted_classes
output_df['proba'] = max_probabilities
output_df.to_csv(file_path, index=False)



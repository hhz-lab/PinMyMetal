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

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'edhmodel_result.csv')

# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()

edh_model = os.path.join(script_directory, 'ml_model_scripts', 'edh_models.joblib')
edh_FFNN = os.path.join(script_directory, 'ml_model_scripts', 'edh_FFNN.h5')

models = load(edh_model)
LR = models['LR']
DT = models['DT']
MLP = models['MLP']
SVM = models['SVM']
FFNN =load_model(edh_FFNN)


sql = """
select id,site_count,ab_angle1,ab_angle2,ab_angle3,ac_angle1,ac_angle2,ac_angle3,bc_angle1,bc_angle2,bc_angle3,
    ma_angle,mb_angle,mc_angle,mo_angle_a,mo_angle_b,mo_angle_c,
    ca_dist_ab,ca_dist_ac,ca_dist_bc,cb_dist_ab,cb_dist_ac,cb_dist_bc,mca_dist_a,mca_dist_b,mca_dist_c,mcb_dist_a,mcb_dist_b,mcb_dist_c,mo_dist_a,mo_dist_b,mo_dist_c,
    c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,
    a_his,a_glu,a_asp,b_his,b_glu,b_asp,c_his,c_glu,c_asp,resitype_ed,resitype_edh,resitype_hh
    from chedh_features where site_count<=3 and sitetype='edh'
"""


cur.execute(sql)
data = cur.fetchall()
cur.close()
conn.close()

raw = pd.DataFrame(data)
raw.set_index([0], inplace=True) # Set the first column (id) as the index

raw.fillna(0, inplace=True)

X_test=raw.values

# Predictions from individual models
LR_probas = LR.predict_proba(X_test)[:, 1]
DT_probas = DT.predict_proba(X_test)[:, 1]
MLP_probas = MLP.predict_proba(X_test)[:, 1]
SVM_probas = SVM.predict_proba(X_test)[:, 1]

FFNN_probas = FFNN.predict(X_test)[:, 0]

# Ensemble predictions (simple averaging)
ensemble_probas = (LR_probas + DT_probas + MLP_probas + SVM_probas + FFNN_probas) / 5

y_pred_result = [1 if proba > 0.5 else 0 for proba in ensemble_probas]

# Combine predictions and probabilities into a single DataFrame
df_combined = pd.DataFrame({
    'result': y_pred_result,
    'proba': ensemble_probas
}, index=raw.index)

# Output to CSV files
df_combined.to_csv(file_path)


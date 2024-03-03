import psycopg2 as pg
import pandas as pd
from joblib import dump, load
from keras.models import load_model
import numpy as np

import getopt
import os, sys

data_dir = "/zinc_prediction/script"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

models = load('models.joblib')

model1 = models['model1']
model2 = models['model2']
model3 = models['model3']
model4 = models['model4']
fcnn =load_model('fcnn.h5')

conn = pg.connect("dbname="+DBname+" password='' port='5432'")
cur = conn.cursor()

sql = "select id,resname_a,resname_b,a_ne2,a_nd1,a_sg,a_ce1,a_cd2,b_ne2,b_nd1,b_sg,b_ce1,b_cd2,\
        atom_dist,ca_dist,cb_dist,angle,angle1,angle2,angle3,\
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,\
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21 from ml_pre"

cur.execute(sql)
data = cur.fetchall()

raw = pd.DataFrame(data)

raw.set_index([0], inplace=True)

X_test=raw.values

y_pred_model1 = model1.predict_proba(X_test)[:, 1]
y_pred_model2 = model2.predict_proba(X_test)[:, 1]
y_pred_model3 = model3.predict_proba(X_test)[:, 1]
y_pred_model4 = model4.predict_proba(X_test)[:, 1]


X_test = X_test.astype(float)
y_pred_fcnn = fcnn.predict(X_test)[:, 0]

y_pred_proba = (y_pred_fcnn + y_pred_model1 + y_pred_model2 + y_pred_model3 + y_pred_model4) / 5
y_pred_result = [1 if proba > 0.5 else 0 for proba in y_pred_proba]

df = pd.DataFrame(data=y_pred_result,
            index=raw.index)

df.columns = ['result']

df_proba = pd.DataFrame(data=y_pred_proba,
            index=raw.index)
df_proba.columns = ['proba']

result_file = data_dir + '/' +pdbid+ '_nb_result' + '/' + 'ml_result.csv'
df.to_csv(result_file)

proba_file = data_dir + '/' + pdbid+ '_nb_result' + '/' + 'mvc_proba.csv'
df_proba.to_csv(proba_file)

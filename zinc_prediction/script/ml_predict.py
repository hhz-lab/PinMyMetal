#from ml_build_model import *
import psycopg2 as pg
import pickle
import pandas as pd

import getopt
import os, sys

data_dir = "/zinc_prediction/script"
loaded_model = pickle.load(open("/zinc_prediction/script/mvc.dat","rb"))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname='"+DBname+"' user='sg_display' password='' port='5432'")
cur = conn.cursor()

sql = "select id,resname_a,resname_b,a_ne2,a_nd1,a_sg,a_ce1,a_cd2,b_ne2,b_nd1,b_sg,b_ce1,b_cd2,\
        atom_dist,ca_dist,cb_dist,angle,angle1,angle2,angle3,\
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,\
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21 from ml_pre"

cur.execute(sql)
data = cur.fetchall()
raw = pd.DataFrame(data)

raw.set_index([0], inplace=True)
y_pred_result = loaded_model.predict(raw.values)

df = pd.DataFrame(data=y_pred_result,
            index=raw.index)


df.columns = ['result']

result_file = data_dir + '/' +pdbid+ '_nb_result' + '/' + 'ml_result.csv'
df.to_csv(result_file)

y_pred_proba = loaded_model.predict_proba(raw.values)[:,1]
df_proba = pd.DataFrame(data=y_pred_proba,
            index=raw.index)
df_proba.columns = ['proba']

proba_file = data_dir + '/' + pdbid+ '_nb_result' + '/' + 'mvc_proba.csv'
df_proba.to_csv(proba_file)

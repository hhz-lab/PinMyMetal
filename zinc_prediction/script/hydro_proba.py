import psycopg2 as pg
import pickle
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
import getopt
import os, sys

data_dir = "/zinc_prediction/script"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname='"+DBname+"' user='sg_display' password='' port='5432'")
cur = conn.cursor()

sql = "select distinct id,tag,c_value,solv,c_value_exp,solv_exp from hydrophobic_pre where site_count >= 3 order by id,tag"

cur.execute(sql)
data = cur.fetchall()

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'hydro_proba.csv','w')

def pearson_correlation(pre_list,exp_list):
    similarity,_=pearsonr(pre_list, exp_list)
    pvalue = (similarity + 1) / 2
    return pvalue

pre_list=[]; exp_list=[]; lastid=0; pre_list2=[]; exp_list2=[]
for i in data:
    pid = i[0]; tag = i[1]; c_value = i[2]; solv = i[3]; c_value_exp = i[4]; solv_exp = i[5]

    if pid == lastid:
        pre_list.append(c_value)
        exp_list.append(c_value_exp)
        pre_list2.append(solv)
        exp_list2.append(solv_exp)
    else:
        if lastid==0:
            pass
        else:
            pearson_c=pearson_correlation(pre_list,exp_list)
            pearson_s=pearson_correlation(pre_list2,exp_list2)
            print('%s,%.3f,%.3f'%(lastid,pearson_c,pearson_s),file=print_log)
        
        lastid=pid
        pre_list=[]; exp_list=[];pre_list2=[]; exp_list2=[]

pearson_c=pearson_correlation(pre_list,exp_list)
pearson_s=pearson_correlation(pre_list2,exp_list2)
print('%s,%.3f,%.3f'%(lastid,pearson_c,pearson_s),file=print_log)







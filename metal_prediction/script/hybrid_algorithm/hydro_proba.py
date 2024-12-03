import pickle
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
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

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'hydro_proba.csv')
print_log = open(file_path,'w')

# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()

sql = "select distinct id,tag,c_value,solv,c_value_exp,solv_exp,sitetype from hydro_pre_chedh \
        where (sitetype='ch' and site_count >= 3) or (sitetype='edh' and site_count >= 4) order by sitetype,id,tag"


cur.execute(sql)
data = cur.fetchall()
cur.close()
conn.close()

raw = pd.DataFrame(data)

def cosine_similarity(pre_list,exp_list):
    similarity=1 - cosine(pre_list, exp_list)
    pvalue = (similarity + 1) / 2
    return pvalue

def pearson_correlation(pre_list,exp_list):
    similarity,_=pearsonr(pre_list, exp_list)
    pvalue = (similarity + 1) / 2
    return pvalue

lastid=0; last_sitetype=0; pre_list=[]; exp_list=[]; pre_list2=[]; exp_list2=[]
for i in data:
    pid = i[0]; tag = i[1]; c_value = i[2]; solv = i[3]; c_value_exp = i[4]; solv_exp = i[5]; sitetype = i[6]

    if pid == lastid and sitetype == last_sitetype:
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
            p_value = (pearson_c + pearson_s)/2
            print('%s,%s,%.3f,%.3f,%.3f'%(lastid,last_sitetype,pearson_c,pearson_s,p_value),file=print_log)
        
        lastid=pid; last_sitetype = sitetype
        pre_list=[]; exp_list=[];pre_list2=[]; exp_list2=[]

# process the last group after the loop ends
if lastid != 0:
    pearson_c = pearson_correlation(pre_list, exp_list)
    pearson_s = pearson_correlation(pre_list2, exp_list2)
    p_value = (pearson_c + pearson_s)/2
    print('%s,%s,%.3f,%.3f,%.3f' % (lastid, last_sitetype, pearson_c, pearson_s,p_value),file=print_log)

print_log.close()



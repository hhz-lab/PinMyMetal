import psycopg2 as pg
import os, sys
from math import *
import getopt
import numpy as np
import itertools

data_dir = "/zinc_prediction/script/"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432'")
cur = conn.cursor()


sql = """ select distinct pdbid,chainid_a, conc_comma(id||'_'||zinc_x||'_'||zinc_y||'_'||zinc_z) as zincs  from \
    (select distinct pdbid,chainid_a,id,zinc_x,zinc_y,zinc_z from zinc_predict_site234_2\
    where id not in (select distinct pre_id from exp_pre_dist where dist <= 2.5)) a group by pdbid,chainid_a order by pdbid,chainid_a """


cur.execute(sql)
data = cur.fetchall()


def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'pre_pre_dist.csv','w')

for i in data:
    pdbid=i[0]; chainid=i[1]; 
    zincs_data=i[2].split(', ')
    zincs_data2=list(itertools.combinations(zincs_data,2))
    if len(zincs_data2) == 0:
        pid_a=i[2].split('_')[0]
        pid_b=-9
        dist_ab=9
        print("exp_exp,%s,%s,%s,%.2f"%(pdbid,pid_a,pid_b,dist_ab),file=print_log)
    else:
        for line in zincs_data2:
            zinc_a=line[0]
            zinc_b=line[1]

            pid_a = zinc_a.split('_')[0]
            a_coord=[0,0,0]
            a_coord[0]=float(zinc_a.split('_')[1])
            a_coord[1]=float(zinc_a.split('_')[2])
            a_coord[2]=float(zinc_a.split('_')[3])

            pid_b = zinc_b.split('_')[0]
            b_coord=[0,0,0]
            b_coord[0]=float(zinc_b.split('_')[1])
            b_coord[1]=float(zinc_b.split('_')[2])
            b_coord[2]=float(zinc_b.split('_')[3])

            dist_ab=dist(a_coord,b_coord)
            print("exp_exp,%s,%s,%s,%.2f"%(pdbid,pid_a,pid_b,dist_ab),file=print_log)
print_log.close()





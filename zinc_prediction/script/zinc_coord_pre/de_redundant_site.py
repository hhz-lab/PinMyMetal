import psycopg2 as pg
import os, sys
from math import *
import itertools
import getopt

data_dir = "/zinc_prediction/script/"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname='"+DBname+"' user='sg_display' password='' port='5432'")
cur = conn.cursor()

sql = "select distinct pdbid, conc_comma(id||'_'||zinc_x||'_'||zinc_y||'_'||zinc_z) as zincs from zinc_predict_site234_2 where deredundant_site is true group by pdbid order by pdbid "

cur.execute(sql)
data = cur.fetchall()

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def merge_list(L):
    lenth = len(L)
    for i in range(1, lenth):
        for j in range(i):
            if L[i] == {0} or L[j] == {0}:
                continue
            x = L[i].union(L[j])
            y = len(L[i]) + len(L[j])
            if len(x) < y:
                L[i] = x
                L[j] = {0}

    return [i for i in L if i != {0}]

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'redundant_site','w')

re_zinc=[]
for r in data:
    pdbid=r[0]
    zincs=r[1]
    zincs_data=zincs.split(', ')
    zincs_data2=list(itertools.combinations(zincs_data,2))

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
        if dist_ab < 2.5:
            pids_set={int(pid_a),int(pid_b)}
            re_zinc.append(pids_set)
    if re_zinc is not None:
        zinc_id=0
        for z in merge_list(re_zinc):
            zinc_id+=1
            print("%s_%s_%s" %(pdbid,z,zinc_id),file=print_log)
print_log.close()


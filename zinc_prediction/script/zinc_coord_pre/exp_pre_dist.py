import psycopg2 as pg
import os, sys
from math import *
import getopt
import numpy as np
from itertools import product

data_dir = "/zinc_prediction/script/"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname='"+DBname+"' user='sg_display' password='' port='5432'")
cur = conn.cursor()


sql = """ select distinct pdbid, conc_comma(site_type||'#'||zincs) as zincs from \
        (select * from exp_pre_site where redundancy is true order by pdbid,site_type)a \
        group by  pdbid  order by pdbid """

cur.execute(sql)
data = cur.fetchall()


def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'exp_pre_dist.csv','w')

for i in data:
    pdbid=i[0]; 
    zincs_data=i[1].split('#')
    
    a_list = zincs_data[1].strip(', pre').split(", ")
    b_list = zincs_data[2].split(", ")
    ab_list = list(product(a_list,b_list))

    for line in ab_list:
        line=list(line)
        zinc_a=line[0]
        zinc_b=line[1]

        a_id=zinc_a.split('_')[0]
        b_id=zinc_b.split('_')[0]

        a_coord=[0,0,0]
        a_coord[0]=float(zinc_a.split('_')[1])
        a_coord[1]=float(zinc_a.split('_')[2])
        a_coord[2]=float(zinc_a.split('_')[3])

        b_coord=[0,0,0]
        b_coord[0]=float(zinc_b.split('_')[1])
        b_coord[1]=float(zinc_b.split('_')[2])
        b_coord[2]=float(zinc_b.split('_')[3])

        dist_ab=dist(a_coord,b_coord)
        print("%s,%s,%s,%.2f"%(pdbid,a_id,b_id,dist_ab),file=print_log)
print_log.close()



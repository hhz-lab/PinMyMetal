import atomium
from atomium import Model, Atom
import psycopg2 as pg
import os, sys
import copy
import hydrophobicity
from hydrophobicity import *
import getopt


data_dir = "/zinc_prediction/script"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432'")
cur = conn.cursor()

sql = "select distinct id,pdbid,zinc_x,zinc_y,zinc_z,chainid_a from zinc_predict_site234_2 order by id "

cur.execute(sql)
data = cur.fetchall()

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'zinc_hydro.csv','w')
for i in data:
    pid = i[0]
    pdb_code = ''
    pdb_code = i[1]

    mid_dir  = pdb_code[1:3]
    pdbfile  = data_dir + '/' + pdb_code + "_nb_result" +"/pdb"+ pdb_code + ".ent"
    
    pdb1 = atomium.open(pdbfile)
    exp_zinc=(i[2],i[3],i[4])
    x=exp_zinc[0]
    y=exp_zinc[1]
    z=exp_zinc[2]
    chainid=i[5].strip()

    c_value = hydrophobic_contrast_function(pdb1,exp_zinc,chainid)
    
    c_value_dict = c_value[0]
    solv_dict = c_value[1]
    for (key,value) in zip(c_value_dict.items(), solv_dict.items()):
        c_r = (key,value)
        print("%s,%s,%s,%s" %(pid,c_r[0][0], c_r[0][1], c_r[1][1]),file=print_log)

print_log.close()

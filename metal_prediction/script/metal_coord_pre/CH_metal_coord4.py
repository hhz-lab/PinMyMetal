from math import *
from Bio.PDB import *
import numpy as np
import os, sys
import copy
from itertools import product
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

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'CH_metal_coord4.csv')

# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()

sql = """select id,
    cd2_coords_a, ne2_coords_a, ce1_coords_a, nd1_coords_a, sg_coords_a,
    cd2_coords_b, ne2_coords_b, ce1_coords_b, nd1_coords_b, sg_coords_b,
    cd2_coords_c, ne2_coords_c, ce1_coords_c, nd1_coords_c, sg_coords_c,
    cd2_coords_d, ne2_coords_d, ce1_coords_d, nd1_coords_d, sg_coords_d,
    resname_a,resname_b,resname_c,resname_d
    from pre_zncu_coord_sites where site_count = 4 order by id"""


cur.execute(sql)
data = cur.fetchall()
cur.close()
conn.close()

def atom_coord(coords_string):
    if not coords_string:
        return [0,0,0]
    #convert_coords_string_to_list
    coords_list = [float(coord) for coord in coords_string.split(',')]
    return coords_list

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1):
    if a_name == 'CYS':
        list_coord = [a_sg]
        atom_name= ['-SG-']
    else:
        list_coord = [a_cd2,a_ne2,a_ce1,a_nd1]
        atom_name = ['-CD2','-NE2','-CE1','-ND1']
    return list_coord,atom_name

def get_coord_atom(a,b):
    list1=a[0]; list2=b[0]; name1=a[1]; name2=b[1]
    coord_list = list(product(list1,list2))
    name_list = list(product(name1,name2))
    dst_d = {}
    for i in range(len(coord_list)):
        dst_d[i] = dist(coord_list[i][0],coord_list[i][1])
    d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
    min_dist = d_order[0][0]

    for i in range(len(coord_list)):
        a_coord = coord_list[min_dist][0]
        b_coord = coord_list[min_dist][1]
    for r in range(len(name_list)):
        atom_a = name_list[min_dist][0]
        atom_b = name_list[min_dist][1]
    return a_coord,b_coord,atom_a,atom_b

class AtomSelect():
    def __init__(self,a,b,c,d):
        self.a_info = a; self.b_info = b; self.c_info = c; self.d_info=d
    def get_best_atom(self):
        a_info=self.a_info; b_info=self.b_info; c_info=self.c_info; d_info=self.d_info
        ab = get_coord_atom(a_info,b_info)
        ac = get_coord_atom(a_info,c_info)
        ad = get_coord_atom(a_info,d_info)
        bc = get_coord_atom(b_info,c_info)
        bd = get_coord_atom(b_info,d_info)
        cd = get_coord_atom(c_info,d_info)

        a1_coord=ab[0]; b1_coord=ab[1]; a1_name=ab[2]; b1_name=ab[3]
        a2_coord=ac[0]; c1_coord=ac[1]; a2_name=ac[2]; c1_name=ac[3]
        b2_coord=bc[0]; c2_coord=bc[1]; b2_name=bc[2]; c2_name=bc[3]
        a3_coord=ad[0]; d1_coord=ad[1]; a3_name=ad[2]; d1_name=ad[3]
        b3_coord=bd[0]; d2_coord=bd[1]; b3_name=bd[2]; d2_name=bd[3]
        c3_coord=cd[0]; d3_coord=cd[1]; c3_name=cd[2]; d3_name=cd[3]

        coord_abcd = list(product([a1_coord,a2_coord,a3_coord],[b1_coord,b2_coord,b3_coord],[c1_coord,c2_coord,c3_coord],[d1_coord,d2_coord,d3_coord]))
        name_abcd = list(product([a1_name,a2_name,a3_name],[b1_name,b2_name,b3_name],[c1_name,c2_name,c3_name],[d1_name,d2_name,d3_name]))
        dist_abcd = {}
        for i in range(len(coord_abcd)):
            dist1=dist(coord_abcd[i][0],coord_abcd[i][1]); dist2=dist(coord_abcd[i][0],coord_abcd[i][2])
            dist3=dist(coord_abcd[i][1],coord_abcd[i][2]); dist4=dist(coord_abcd[i][0],coord_abcd[i][3])
            dist5=dist(coord_abcd[i][1],coord_abcd[i][3]); dist6=dist(coord_abcd[i][2],coord_abcd[i][3])
            dist_abcd[i] = dist1+dist2+dist3+dist4+dist5+dist6
        d_order=sorted(dist_abcd.items(),key=lambda x:x[1],reverse=False)
        min_abcd = d_order[0][0]
        for i in range(len(coord_abcd)):
            a_coord = coord_abcd[min_abcd][0]
            b_coord = coord_abcd[min_abcd][1]
            c_coord = coord_abcd[min_abcd][2]
            d_coord = coord_abcd[min_abcd][3]
        for r in range(len(name_abcd)):
            atom_a = name_abcd[min_abcd][0]
            atom_b = name_abcd[min_abcd][1]
            atom_c = name_abcd[min_abcd][2]
            atom_d = name_abcd[min_abcd][3]
        return a_coord,b_coord,c_coord,d_coord,atom_a,atom_b,atom_c,atom_d


def metal_coord(f_a,f_b,f_c,f_d):
    f_metal=[0,0,0]
    f_metal[0] = (f_a[0]+f_b[0]+f_c[0]+f_d[0])/4
    f_metal[1] = (f_a[1]+f_b[1]+f_c[1]+f_d[1])/4
    f_metal[2] = (f_a[2]+f_b[2]+f_c[2]+f_d[2])/4
    return f_metal

sql = """select id,
    cd2_coords_a, ne2_coords_a, ce1_coords_a, nd1_coords_a, sg_coords_a,
    cd2_coords_b, ne2_coords_b, ce1_coords_b, nd1_coords_b, sg_coords_b,
    cd2_coords_c, ne2_coords_c, ce1_coords_c, nd1_coords_c, sg_coords_c,
    cd2_doords_d, ne2_doords_d, ce1_doords_d, nd1_doords_d, sg_doords_d,
    resname_a,resname_b,resname_c,resname_d
    from pre_zncu_coord_sites where site_count = 4 order by"""

print_log = open(file_path,'w')

for i in data:
    pid=i[0]

    a_cd2=atom_coord(i[1]);  a_ne2=atom_coord(i[2]);  a_ce1=atom_coord(i[3]);  a_nd1=atom_coord(i[4]); a_sg=atom_coord(i[5])
    b_cd2=atom_coord(i[6]);  b_ne2=atom_coord(i[7]);  b_ce1=atom_coord(i[8]);  b_nd1=atom_coord(i[9]); b_sg=atom_coord(i[10])
    c_cd2=atom_coord(i[11]);  c_ne2=atom_coord(i[12]);  c_ce1=atom_coord(i[13]);  c_nd1=atom_coord(i[14]); c_sg=atom_coord(i[15])
    d_cd2=atom_coord(i[16]);  d_ne2=atom_coord(i[17]);  d_ce1=atom_coord(i[18]);  d_nd1=atom_coord(i[19]); d_sg=atom_coord(i[20])

    a_name=i[21];b_name=i[22];c_name=i[23];d_name=i[24]
    
    a = get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1)
    b = get_atom(b_name,b_sg,b_cd2,b_ne2,b_ce1,b_nd1)
    c = get_atom(c_name,c_sg,c_cd2,c_ne2,c_ce1,c_nd1)
    d = get_atom(d_name,d_sg,d_cd2,d_ne2,d_ce1,d_nd1)
    abcd = AtomSelect(a,b,c,d).get_best_atom()
    a_coord=abcd[0]; b_coord=abcd[1]; c_coord=abcd[2]; d_coord=abcd[3] 
    atom_a=abcd[4]; atom_b=abcd[5]; atom_c=abcd[6]; atom_d=abcd[7] 
    
    f_metal=metal_coord(a_coord,b_coord,c_coord,d_coord)

    dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal); dist_dm = dist(d_coord,f_metal)

    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%s,%8.3f,%8.3f,%8.3f,%8.3f" %(pid,f_metal[0],f_metal[1],f_metal[2],atom_a,atom_b,atom_c,atom_d,dist_am,dist_bm,dist_cm,dist_dm),file=print_log)

print_log.close()




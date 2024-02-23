from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
from itertools import product
import getopt

data_dir = "/zinc_prediction/script/"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)


conn = pg.connect("dbname="+DBname+" password='' port='5432'")
cur = conn.cursor()


sql = """select id,
    cd2_x,cd2_y,cd2_z,ne2_x,ne2_y,ne2_z,ce1_x,ce1_y,ce1_z,nd1_x,nd1_y,nd1_z,
    cd2_x1,cd2_y1,cd2_z1,ne2_x1,ne2_y1,ne2_z1,ce1_x1,ce1_y1,ce1_z1,nd1_x1,nd1_y1,nd1_z1,
    cd2_x2,cd2_y2,cd2_z2,ne2_x2,ne2_y2,ne2_z2,ce1_x2,ce1_y2,ce1_z2,nd1_x2,nd1_y2,nd1_z2,
    cd2_x3,cd2_y3,cd2_z3,ne2_x3,ne2_y3,ne2_z3,ce1_x3,ce1_y3,ce1_z3,nd1_x3,nd1_y3,nd1_z3,
    sg_x,sg_y,sg_z,sg_x1,sg_y1,sg_z1,sg_x2,sg_y2,sg_z2,sg_x3,sg_y3,sg_z3,
    resname_a,resname_b,resname_c,resname_d
    from pre_zinc_coord_site234 where resname_d != 'x' order by id """


cur.execute(sql)
data = cur.fetchall()


def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1):
    if a_name == 'CYS':
        list_coord = [a_sg]
        atom_name= ['a_-SG-']
    else:
        list_coord = [a_cd2,a_ne2,a_ce1,a_nd1]
        atom_name = ['a_-CD2','a_-NE2','a_-CE1','a_-ND1']
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
        atom_a = name_list[min_dist][0].split('_')[1]
        atom_b = name_list[min_dist][1].split('_')[1]
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


def zinc_coord(f_a,f_b,f_c,f_d):
    f_zinc=[0,0,0]
    f_zinc[0] = (f_a[0]+f_b[0]+f_c[0]+f_d[0])/4
    f_zinc[1] = (f_a[1]+f_b[1]+f_c[1]+f_d[1])/4
    f_zinc[2] = (f_a[2]+f_b[2]+f_c[2]+f_d[2])/4
    return f_zinc

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'zinc_coord4.csv','w')

for i in data:
    pid=i[0]

    a_cd2=[0,0,0]; a_ne2=[0,0,0]; a_ce1=[0,0,0]; a_nd1=[0,0,0]

    b_cd2=[0,0,0]; b_ne2=[0,0,0]; b_ce1=[0,0,0]; b_nd1=[0,0,0]

    c_cd2=[0,0,0]; c_ne2=[0,0,0]; c_ce1=[0,0,0]; c_nd1=[0,0,0]

    d_cd2=[0,0,0]; d_ne2=[0,0,0]; d_ce1=[0,0,0]; d_nd1=[0,0,0]

    a_sg=[0,0,0]; b_sg=[0,0,0]; c_sg=[0,0,0]; d_sg=[0,0,0]

    a_cd2[0]=i[1];a_cd2[1]=i[2];a_cd2[2]=i[3];a_ne2[0]=i[4];a_ne2[1]=i[5];a_ne2[2]=i[6];a_ce1[0]=i[7];a_ce1[1]=i[8];a_ce1[2]=i[9];a_nd1[0]=i[10];a_nd1[1]=i[11];a_nd1[2]=i[12]
    b_cd2[0]=i[13];b_cd2[1]=i[14];b_cd2[2]=i[15];b_ne2[0]=i[16];b_ne2[1]=i[17];b_ne2[2]=i[18];b_ce1[0]=i[19];b_ce1[1]=i[20];b_ce1[2]=i[21];b_nd1[0]=i[22];b_nd1[1]=i[23];b_nd1[2]=i[24]
    c_cd2[0]=i[25];c_cd2[1]=i[26];c_cd2[2]=i[27];c_ne2[0]=i[28];c_ne2[1]=i[29];c_ne2[2]=i[30];c_ce1[0]=i[31];c_ce1[1]=i[32];c_ce1[2]=i[33];c_nd1[0]=i[34];c_nd1[1]=i[35];c_nd1[2]=i[36]
    d_cd2[0]=i[37];d_cd2[1]=i[38];d_cd2[2]=i[39];d_ne2[0]=i[40];d_ne2[1]=i[41];d_ne2[2]=i[42];d_ce1[0]=i[43];d_ce1[1]=i[44];d_ce1[2]=i[45];d_nd1[0]=i[46];d_nd1[1]=i[47];d_nd1[2]=i[48]
    a_sg[0]=i[49];a_sg[1]=i[50];a_sg[2]=i[51];b_sg[0]=i[52];b_sg[1]=i[53];b_sg[2]=i[54];c_sg[0]=i[55];c_sg[1]=i[56];c_sg[2]=i[57];d_sg[0]=i[58];d_sg[1]=i[59];d_sg[2]=i[60]
    a_name=i[61];b_name=i[62];c_name=i[63];d_name=i[64]

    
    a = get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1)
    b = get_atom(b_name,b_sg,b_cd2,b_ne2,b_ce1,b_nd1)
    c = get_atom(c_name,c_sg,c_cd2,c_ne2,c_ce1,c_nd1)
    d = get_atom(d_name,d_sg,d_cd2,d_ne2,d_ce1,d_nd1)
    abcd = AtomSelect(a,b,c,d).get_best_atom()
    a_coord=abcd[0]; b_coord=abcd[1]; c_coord=abcd[2]; d_coord=abcd[3] 
    atom_a=abcd[4]; atom_b=abcd[5]; atom_c=abcd[6]; atom_d=abcd[7] 
    
    f_zinc=zinc_coord(a_coord,b_coord,c_coord,d_coord)
    dist_af = dist(f_zinc,a_coord); dist_bf = dist(f_zinc,b_coord); dist_cf = dist(f_zinc,c_coord); dist_df = dist(f_zinc,d_coord)

    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%s,%8.2f,%8.2f,%8.2f,%8.2f" %(pid,f_zinc[0],f_zinc[1],f_zinc[2],atom_a,atom_b,atom_c,atom_d,dist_af,dist_bf,dist_cf,dist_df),file=print_log)

print_log.close()




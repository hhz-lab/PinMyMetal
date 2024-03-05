from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
from itertools import product
import getopt

data_dir = "../"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

sql = """select id,
    cg_x,cg_y,cg_z,cd2_x,cd2_y,cd2_z,ne2_x,ne2_y,ne2_z,ce1_x,ce1_y,ce1_z,nd1_x,nd1_y,nd1_z,
    cg_x1,cg_y1,cg_z1,cd2_x1,cd2_y1,cd2_z1,ne2_x1,ne2_y1,ne2_z1,ce1_x1,ce1_y1,ce1_z1,nd1_x1,nd1_y1,nd1_z1,
    cg_x2,cg_y2,cg_z2,cd2_x2,cd2_y2,cd2_z2,ne2_x2,ne2_y2,ne2_z2,ce1_x2,ce1_y2,ce1_z2,nd1_x2,nd1_y2,nd1_z2,
    site_count
    from pre_zinc_coord_site234 where resi_type= 'H_H' and site_count < 4 order by id """


cur.execute(sql)
data = cur.fetchall()

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def get_atom(a_cd2,a_ne2,a_ce1,a_nd1):
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

    a_coord = coord_list[min_dist][0]
    b_coord = coord_list[min_dist][1]
    atom_a = name_list[min_dist][0].split('_')[1]
    atom_b = name_list[min_dist][1].split('_')[1]
    return a_coord,b_coord,atom_a,atom_b

class AtomSelect():
    def __init__(self,a,b,c):
        self.a_info = a; self.b_info = b; self.c_info = c
    def get_best_atom(self):
        a_info=self.a_info; b_info=self.b_info; c_info=self.c_info
        ab = get_coord_atom(a_info,b_info)
        ac = get_coord_atom(a_info,c_info)
        bc = get_coord_atom(b_info,c_info)
        a1_coord=ab[0]; b1_coord=ab[1]; a1_name=ab[2]; b1_name=ab[3]
        a2_coord=ac[0]; c1_coord=ac[1]; a2_name=ac[2]; c1_name=ac[3]
        b2_coord=bc[0]; c2_coord=bc[1]; b2_name=bc[2]; c2_name=bc[3]

        coord_abc = list(product([a1_coord,a2_coord],[b1_coord,b2_coord],[c1_coord,c2_coord]))
        name_abc = list(product([a1_name,a2_name],[b1_name,b2_name],[c1_name,c2_name]))
        dist_abc = {}
        for i in range(len(coord_abc)):
            dist_abc[i] = dist(coord_abc[i][0],coord_abc[i][1]) + dist(coord_abc[i][0],coord_abc[i][2]) + dist(coord_abc[i][1],coord_abc[i][2])
        d_order=sorted(dist_abc.items(),key=lambda x:x[1],reverse=False)
        min_abc = d_order[0][0]
        for i in range(len(coord_abc)):
            a_coord = coord_abc[min_abc][0]
            b_coord = coord_abc[min_abc][1]
            c_coord = coord_abc[min_abc][2]
        for r in range(len(name_abc)):
            atom_a = name_abc[min_abc][0]
            atom_b = name_abc[min_abc][1]
            atom_c = name_abc[min_abc][2]
        return a_coord,b_coord,c_coord,atom_a,atom_b,atom_c


def f_coord(a,b):
    x1=a[0]; y1=a[1]; z1=a[2];
    x2=b[0]; y2=b[1]; z2=b[2];
    # the angle of CG_NE2_0
    a1=(x2-x1)   #CG-NE2/CE2 vector
    b1=(y2-y1)
    c1=(z2-z1)
    a2=0
    b2=0
    c2=0
    f=[0,0,0]
    if x2 > x1:
        a2=2
        m=(a1*a2+b1*b2+c1*c2)/(sqrt(a1*a1+b1*b1+c1*c1)*sqrt(a2*a2+b2*b2+c2*c2))
        angle=(acos(m))
        f[0]=x2+abs(2.1*cos(angle))
        f[1]=((f[0]-x1)*(y2-y1))/(x2-x1) + y1
        f[2]=((f[0]-x1)*(z2-z1))/(x2-x1) + z1
    elif x2 < x1:
        a2=-2
        m=(a1*a2+b1*b2+c1*c2)/(sqrt(a1*a1+b1*b1+c1*c1)*sqrt(a2*a2+b2*b2+c2*c2))
        angle=(acos(m))
        f[0]=x2-abs(2.1*cos(angle))
        f[1]=((f[0]-x1)*(y2-y1))/(x2-x1) + y1
        f[2]=((f[0]-x1)*(z2-z1))/(x2-x1) + z1
    return f

print_log = open(data_dir +pdbid+ '_nb_result' + '/' + 'His_zinc_coord_23.csv','w')

for i in data:
    a_cg=[0,0,0]
    a_cd2=[0,0,0]
    a_ne2=[0,0,0]
    a_ce1=[0,0,0]
    a_nd1=[0,0,0]

    b_cg=[0,0,0]
    b_cd2=[0,0,0]
    b_ne2=[0,0,0]
    b_ce1=[0,0,0]
    b_nd1=[0,0,0]

    c_cg=[0,0,0]
    c_cd2=[0,0,0]
    c_ne2=[0,0,0]
    c_ce1=[0,0,0]
    c_nd1=[0,0,0]

    pid=i[0]

    a_cg[0]=i[1]
    a_cg[1]=i[2]
    a_cg[2]=i[3]
    a_cd2[0]=i[4]
    a_cd2[1]=i[5]
    a_cd2[2]=i[6]
    a_ne2[0]=i[7]
    a_ne2[1]=i[8]
    a_ne2[2]=i[9]
    a_ce1[0]=i[10]
    a_ce1[1]=i[11]
    a_ce1[2]=i[12]
    a_nd1[0]=i[13]
    a_nd1[1]=i[14]
    a_nd1[2]=i[15]

    b_cg[0]=i[16]
    b_cg[1]=i[17]
    b_cg[2]=i[18]
    b_cd2[0]=i[19]
    b_cd2[1]=i[20]
    b_cd2[2]=i[21]
    b_ne2[0]=i[22]
    b_ne2[1]=i[23]
    b_ne2[2]=i[24]
    b_ce1[0]=i[25]
    b_ce1[1]=i[26]
    b_ce1[2]=i[27]
    b_nd1[0]=i[28]
    b_nd1[1]=i[29]
    b_nd1[2]=i[30]

    c_cg[0]=i[31]
    c_cg[1]=i[32]
    c_cg[2]=i[33]
    c_cd2[0]=i[34]
    c_cd2[1]=i[35]
    c_cd2[2]=i[36]
    c_ne2[0]=i[37]
    c_ne2[1]=i[38]
    c_ne2[2]=i[39]
    c_ce1[0]=i[40]
    c_ce1[1]=i[41]
    c_ce1[2]=i[42]
    c_nd1[0]=i[43]
    c_nd1[1]=i[44]
    c_nd1[2]=i[45]

    site_count=i[46]

    a_CG=[0,0,0]
    b_CG=[0,0,0]
    c_CG=[0,0,0]

    a_CG[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3)
    a_CG[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3)
    a_CG[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)
    
    b_CG[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3)
    b_CG[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3)
    b_CG[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)

    a = get_atom(a_cd2,a_ne2,a_ce1,a_nd1)
    b = get_atom(b_cd2,b_ne2,b_ce1,b_nd1)
    c = get_atom(c_cd2,c_ne2,c_ce1,c_nd1)
    
    f_zinc = [0,0,0]
    if site_count == 2:
        a_coord = get_coord_atom(a,b)[0];b_coord = get_coord_atom(a,b)[1]
        atom_a = get_coord_atom(a,b)[2]; atom_b = get_coord_atom(a,b)[3]; atom_c = 'x'
        fa= f_coord(a_CG,a_coord); fb= f_coord(b_CG,b_coord)
        f_zinc[0] = (fa[0]+fb[0])/2
        f_zinc[1] = (fa[1]+fb[1])/2
        f_zinc[2] = (fa[2]+fb[2])/2
        dist_af = dist(f_zinc,a_coord); dist_bf = dist(f_zinc,b_coord); dist_cf =-9999

    else:
        c_CG[0]=round((c_cg[0]+c_cd2[0]+c_ne2[0]+c_ce1[0]+c_nd1[0])/5,3)
        c_CG[1]=round((c_cg[1]+c_cd2[1]+c_ne2[1]+c_ce1[1]+c_nd1[1])/5,3)
        c_CG[2]=round((c_cg[2]+c_cd2[2]+c_ne2[2]+c_ce1[2]+c_nd1[2])/5,3)
        abc = AtomSelect(a,b,c).get_best_atom()
        a_coord=abc[0]; b_coord=abc[1]; c_coord=abc[2]; atom_a=abc[3]; atom_b=abc[4]; atom_c=abc[5]

        fa= f_coord(a_CG,a_coord); fb= f_coord(b_CG,b_coord);  fc= f_coord(c_CG,c_coord)
        f_zinc[0] = (fa[0]+fb[0]+fc[0])/3
        f_zinc[1] = (fa[1]+fb[1]+fc[1])/3
        f_zinc[2] = (fa[2]+fb[2]+fc[2])/3
        dist_af = dist(f_zinc,a_coord); dist_bf = dist(f_zinc,b_coord); dist_cf = dist(f_zinc,c_coord)
    
    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%8.2f,%8.2f,%8.2f" %(pid,f_zinc[0],f_zinc[1],f_zinc[2],atom_a,atom_b,atom_c,dist_af,dist_bf,dist_cf),file=print_log)

print_log.close()


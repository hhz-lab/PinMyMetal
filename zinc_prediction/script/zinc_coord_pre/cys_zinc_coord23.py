from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
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

sql = """select id,sg_x,sg_y,sg_z,sg_x1,sg_y1,sg_z1,sg_x2,sg_y2,sg_z2,
    cb_x_a,cb_y_a,cb_z_a,cb_x_b,cb_y_b,cb_z_b,cb_x_c,cb_y_c,cb_z_c,
    resname_c
    from pre_zinc_coord_site234 where resi_type= 'C_C' and resname_d = 'x' order by id"""


cur.execute(sql)
data = cur.fetchall()

def dist(x,y):
    return sqrt((x.coord[0]-y.coord[0])**2 + (x.coord[1]-y.coord[1])**2 + (x.coord[2]-y.coord[2])**2)

def dist2(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def get_matrix(a1,b1,e1):

    a2=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C')
    b2=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C')
    e2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')

    a2.coord[0]=0
    a2.coord[1]=0
    a2.coord[2]=dist(a1,b1)/2

    b2.coord[0]=0
    b2.coord[1]=0
    b2.coord[2]=-dist(a1,b1)/2

    fixed  = [a2,b2,e2]
    moving = [a1,b1,e1]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup.rotran

def rotran(c1,matrix) :
    matrix_r = np.mat(matrix[0])
    matrix_t = matrix[1]
    matrix_t = np.ravel(matrix_t).tolist()
    t1=matrix_t[0]
    t2=matrix_t[1]
    t3=matrix_t[2]

    c2 = (c1.coord)*matrix_r
    c2_list = np.concatenate(c2).ravel().tolist()
    c2 = c2_list[0]

    c2[0] =c2[0] + t1
    c2[1] =c2[1] + t2
    c2[2] =c2[2] + t3
    c2 = [round(i,3) for i in c2]
    return c2

def get_r_matrix(a1, b1, e1):
    a2=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C')
    b2=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C')
    e2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
    a2.coord = [0,0,dist(a1,b1)/2]
    b2.coord = [0,0,-dist(a1,b1)/2]
    moving  = [a2,b2,e2]
    fixed = [a1,b1,e1]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup.rotran

def r_rotran(c2,r_matrix) :
    r_matrix_r = np.mat(r_matrix[0])
    r_matrix_t = r_matrix[1]
    c1 = c2*r_matrix_r
    c1_list = np.concatenate(c1).ravel().tolist()
    c1 = c1_list[0]

    c1[0] =c1[0] + r_matrix_t[0]
    c1[1] =c1[1] + r_matrix_t[1]
    c1[2] =c1[2] + r_matrix_t[2]
    c1 = [round(i,3) for i in c1]
    return c1


def angle_AF(A,B,C):
    c = dist2(A,B)
    a = dist2(B,C)
    b = dist2(A,C)
    return (acos((a*a + c*c - b*b)/(2*a*c)))*180/pi

def score1(f1,a1,b1,c1,d1):
    score_c=abs(angle_AF(f1,a1.coord,c1.coord)-109)
    score_d=abs(angle_AF(f1,b1.coord,d1.coord)-109)
    score=score_c+score_d
    return score

def angle_min(l,a1,b1,c1,d1):
    f1_list=l
    score=score1(f1_list[0],a1,b1,c1,d1)
    for i in range(0,360):
        scorei=score1(f1_list[i],a1,b1,c1,d1)
        if scorei <= score:
            score=scorei
            f1_coord=f1_list[i]
    return f1_coord


class ZincCoord:
    def __init__(self,a,b,c,d):
        self.a1 = a
        self.b1 = b
        self.c1 = c
        self.d1 = d
        
    def f1_coord(self):
        a1 = self.a1
        b1 = self.b1
        c1 = self.c1
        d1 = self.d1

        e1=Atom.Atom('CE',[0,0,0], 20, 1, None, 'C1', 1, element='C')
        e1.coord[0] = (a1.coord[0]+b1.coord[0])/2
        e1.coord[1] = (a1.coord[1]+b1.coord[1])/2
        e1.coord[2] = (a1.coord[2]+b1.coord[2])/2
        
        matrix=get_matrix(a1, b1, e1)
        r_matrix = get_r_matrix(a1, b1, e1)
        a2 = rotran(a1,matrix)  
        b2 = rotran(b1,matrix)
        c2 = rotran(c1,matrix)
        d2 = rotran(d1,matrix)
    
        f1_list=[]
        for i in range(1,361):
            B = i*(pi/180)
            f2 = [1.2*sin(B), 1.2*cos(B), 0]

            score_c=abs(angle_AF(f2,a2,c2)-109)  
            score_d=abs(angle_AF(f2,b2,d2)-109)
            score=score_c+score_d

            f1=r_rotran(f2,r_matrix)  
            f1_list.append(f1)
        f1=angle_min(f1_list,a1,b1,c1,d1) 
        return f1

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'cys_zinc_coord23.csv','w')

for i in data:
    a1=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C')
    b1=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C')
    a2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
    b2=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')
    c1=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C')
    c2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')

    
    pid=i[0]

    a1.coord[0]=i[1]
    a1.coord[1]=i[2]
    a1.coord[2]=i[3]
    b1.coord[0]=i[4]
    b1.coord[1]=i[5]
    b1.coord[2]=i[6]
    c1.coord[0]=i[7]
    c1.coord[1]=i[8]
    c1.coord[2]=i[9]

    a_Cb_x=i[10]
    a_Cb_y=i[11]
    a_Cb_z=i[12]
    b_Cb_x=i[13]
    b_Cb_y=i[14]
    b_Cb_z=i[15]
    c_Cb_x=i[16]
    c_Cb_y=i[17]
    c_Cb_z=i[18]

    resname_c=i[19]

    a2.coord[0] = a_Cb_x
    a2.coord[1] = a_Cb_y
    a2.coord[2] = a_Cb_z
    b2.coord[0] = b_Cb_x
    b2.coord[1] = b_Cb_y
    b2.coord[2] = b_Cb_z
    c2.coord[0] = c_Cb_x
    c2.coord[1] = c_Cb_y
    c2.coord[2] = c_Cb_z
   
    
    f_ab = ZincCoord(a1,b1,a2,b2).f1_coord()
    if resname_c == 'x':
        f_zinc = f_ab
        atom_c = 'x'
        dist_cf =-9999
    else:
        f_ac = ZincCoord(a1,c1,a2,c2).f1_coord()
        f_bc = ZincCoord(b1,c1,b2,c2).f1_coord()

        f_zinc = [0,0,0]
        atom_c = '-SG-'; dist_cf = dist2(f_zinc,c1.coord)
        d_f_ab_ac=dist2(f_ab,f_ac)
        d_f_ab_bc=dist2(f_ab,f_bc)
        d_f_ac_bc=dist2(f_ac,f_bc)
        zinc_list = [(f_ab,f_ac),(f_ab,f_bc),(f_ac,f_bc)]
        dst_d = {}
        for i in range(len(zinc_list)):
            dst_d[i] = dist2(zinc_list[i][0],zinc_list[i][1])
        d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
        min_dist = d_order[0][0]
        a = zinc_list[min_dist][0]
        b = zinc_list[min_dist][1]
        f_zinc[0] =(a[0]+b[0])/2
        f_zinc[1] =(a[1]+b[1])/2
        f_zinc[2] =(a[2]+b[2])/2

    atom_a = '-SG-'
    atom_b = '-SG-'
    dist_af = dist2(f_zinc,a1.coord); dist_bf = dist2(f_zinc,b1.coord)

    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%8.2f,%8.2f,%8.2f"%(pid,f_zinc[0],f_zinc[1],f_zinc[2],atom_a,atom_b,atom_c,dist_af,dist_bf,dist_cf),file=print_log)

print_log.close()
    

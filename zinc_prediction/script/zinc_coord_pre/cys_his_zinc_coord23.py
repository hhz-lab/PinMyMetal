from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
import His_3
import cys_3
import cys_his_3
import cys_his_2
from itertools import product
import getopt

data_dir = '/zinc_prediction/script/'

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432'")
cur = conn.cursor()


sql = """select id,
    cg_x,cg_y,cg_z,cd2_x,cd2_y,cd2_z,ne2_x,ne2_y,ne2_z,ce1_x,ce1_y,ce1_z,nd1_x,nd1_y,nd1_z,
    cg_x1,cg_y1,cg_z1,cd2_x1,cd2_y1,cd2_z1,ne2_x1,ne2_y1,ne2_z1,ce1_x1,ce1_y1,ce1_z1,nd1_x1,nd1_y1,nd1_z1,
    cg_x2,cg_y2,cg_z2,cd2_x2,cd2_y2,cd2_z2,ne2_x2,ne2_y2,ne2_z2,ce1_x2,ce1_y2,ce1_z2,nd1_x2,nd1_y2,nd1_z2,
    sg_x,sg_y,sg_z,sg_x1,sg_y1,sg_z1,sg_x2,sg_y2,sg_z2,
    ca_x_a,ca_y_a,ca_z_a,cb_x_a,cb_y_a,cb_z_a,
    ca_x_b,ca_y_b,ca_z_b,cb_x_b,cb_y_b,cb_z_b,
    ca_x_c,ca_y_c,ca_z_c,cb_x_c,cb_y_c,cb_z_c,
    resname_a,resname_b,resname_c,
    ab_type,ac_type,bc_type
    from pre_zinc_coord_site234 where resi_type='C_H' and resname_d='x' order by id """


cur.execute(sql)
data = cur.fetchall()

def get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1):
    if a_name == 'CYS':
        list_coord = [a_sg]
        atom_name= ['a_-SG-']
    else:
        list_coord = [a_cd2,a_ne2,a_ce1,a_nd1]
        atom_name = ['a_-CD2','a_-NE2','a_-CE1','a_-ND1']
    return list_coord,atom_name


def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

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


class ZincCoords():
    def __init__(self,rt,cg,cd,ne,ce,nd,cg2,cd2,ne2,ce2,nd2,cbx,cby,cbz,cbx2,cby2,cbz2,a1,b1,a1_ch,b1_ch,atom,cax,cay,caz,cax2,cay2,caz2,a_coord,atom_a,b_coord,atom_b):
        self.resi_type=rt
        self.a_cg=cg; self.a_cd2=cd; self.a_ne2=ne; self.a_ce1=ce; self.a_nd1=nd
        self.b_cg=cg2; self.b_cd2=cd2; self.b_ne2=ne2; self.b_ce1=ce2; self.b_nd1=nd2
        self.a_Cb_x=cbx; self.a_Cb_y=cby; self.a_Cb_z=cbz; self.b_Cb_x=cbx2; self.b_Cb_y=cby2; self.b_Cb_z=cbz2
        self.a1=a1;self.b1=b1
        self.a1_ch=a1_ch; self.b1_ch=b1_ch; self.resname_b=atom
        self.a_Ca_x=cax;self.a_Ca_y=cay;self.a_Ca_z=caz;self.b_Ca_x=cax2;self.b_Ca_y=cay2;self.b_Ca_z=caz2
        self.a_coord=a_coord; self.atom_a=atom_a; self.b_coord=b_coord; self.atom_b=atom_b

    def zinc_coore(self):
        resi_type=self.resi_type
        a_cg=self.a_cg; a_cd2=self.a_cd2; a_ne2=self.a_ne2; a_ce1=self.a_ce1; a_nd1=self.a_nd1
        b_cg=self.b_cg; b_cd2=self.b_cd2; b_ne2=self.b_ne2; b_ce1=self.b_ce1; b_nd1=self.b_nd1
        a_Cb_x=self.a_Cb_x; a_Cb_y=self.a_Cb_y; a_Cb_z=self.a_Cb_z; b_Cb_x=self.b_Cb_x; b_Cb_y=self.b_Cb_y; b_Cb_z=self.b_Cb_z
        a1=self.a1;b1=self.b1
        a1_ch=self.a1_ch; b1_ch=self.b1_ch; resname_b=self.resname_b
        a_Ca_x=self.a_Ca_x;a_Ca_y=self.a_Ca_y;a_Ca_z=self.a_Ca_z;b_Ca_x=self.b_Ca_x;b_Ca_y=self.b_Ca_y;b_Ca_z=self.b_Ca_z
        a_coord=self.a_coord; atom_a=self.atom_a; b_coord=self.b_coord; atom_b=self.atom_b

        if resi_type == 'H_H':
            
            a_CG[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3); a_CG[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3); a_CG[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)

            b_CG[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3); b_CG[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3); b_CG[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)

            fa= His_3.f_coord(a_CG,a_coord); fb= His_3.f_coord(b_CG,b_coord)
            f_zinc = His_3.f1_zinc(fa,fb,atom_a,atom_b)

        elif resi_type == 'C_C':
            a2.coord[0] = a_Cb_x; a2.coord[1] = a_Cb_y; a2.coord[2] = a_Cb_z; b2.coord[0] = b_Cb_x; b2.coord[1] = b_Cb_y; b2.coord[2] = b_Cb_z;
            f_zinc = cys_3.ZincCoord(a1,b1,a2,b2).f1_coord()
        else:
            f_zinc = cys_his_3.ZincCoord2(a1_ch,b1_ch,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_Ca_x,a_Ca_y,a_Ca_z,b_Ca_x,b_Ca_y,b_Ca_z,a_coord,atom_a,b_coord,atom_b).f_zinc_result()

        return f_zinc

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'cys_his_zinc_coord23.csv','w')

for i in data:
    pid=i[0];
    # H_H
    a_cg=[0,0,0]; a_cd2=[0,0,0]; a_ne2=[0,0,0]; a_ce1=[0,0,0]; a_nd1=[0,0,0]

    b_cg=[0,0,0]; b_cd2=[0,0,0]; b_ne2=[0,0,0]; b_ce1=[0,0,0]; b_nd1=[0,0,0]

    c_cg=[0,0,0]; c_cd2=[0,0,0]; c_ne2=[0,0,0]; c_ce1=[0,0,0]; c_nd1=[0,0,0]

    a_cg[0]=i[1]; a_cg[1]=i[2]; a_cg[2]=i[3]; a_cd2[0]=i[4]; a_cd2[1]=i[5]; a_cd2[2]=i[6]; a_ne2[0]=i[7]; a_ne2[1]=i[8]; a_ne2[2]=i[9]
    a_ce1[0]=i[10]; a_ce1[1]=i[11]; a_ce1[2]=i[12]; a_nd1[0]=i[13]; a_nd1[1]=i[14];a_nd1[2]=i[15]

    b_cg[0]=i[16]; b_cg[1]=i[17]; b_cg[2]=i[18]; b_cd2[0]=i[19]; b_cd2[1]=i[20]; b_cd2[2]=i[21]; b_ne2[0]=i[22]; b_ne2[1]=i[23]; b_ne2[2]=i[24]
    b_ce1[0]=i[25]; b_ce1[1]=i[26]; b_ce1[2]=i[27]; b_nd1[0]=i[28]; b_nd1[1]=i[29];b_nd1[2]=i[30]

    c_cg[0]=i[31]; c_cg[1]=i[32]; c_cg[2]=i[33]; c_cd2[0]=i[34]; c_cd2[1]=i[35]; c_cd2[2]=i[36]; c_ne2[0]=i[37]; c_ne2[1]=i[38]; c_ne2[2]=i[39]
    c_ce1[0]=i[40]; c_ce1[1]=i[41]; c_ce1[2]=i[42]; c_nd1[0]=i[43]; c_nd1[1]=i[44];c_nd1[2]=i[45]

    # C_C
    a1=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C'); a2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
    b1=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C'); b2=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')
    c1=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C1', 1, element='C'); c2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')

    a1.coord[0]=i[46]; a1.coord[1]=i[47]; a1.coord[2]=i[48]; b1.coord[0]=i[49]; b1.coord[1]=i[50]; b1.coord[2]=i[51]; c1.coord[0]=i[52]; c1.coord[1]=i[53]; c1.coord[2]=i[54]

    a_Cb_x=i[58]; a_Cb_y=i[59]; a_Cb_z=i[60]; b_Cb_x=i[64]; b_Cb_y=i[65]; b_Cb_z=i[66]; c_Cb_x=i[70]; c_Cb_y=i[71]; c_Cb_z=i[72];


    # C_H
    a1_ch=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C'); a2_ch=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
    b1_ch=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C'); b2_ch=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')
    c1_ch=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C1', 1, element='C'); c2_ch=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')

    resname_a=i[73].strip();resname_b=i[74].strip();resname_c=i[75].strip()
    
    a1_ch.coord[0]=i[46]; a1_ch.coord[1]=i[47]; a1_ch.coord[2]=i[48]; b1_ch.coord[0]=i[49]; b1_ch.coord[1]=i[50]; b1_ch.coord[2]=i[51]; c1_ch.coord[0]=i[52]; c1_ch.coord[1]=i[53]; c1_ch.coord[2]=i[54]
    
    a_Ca_x=i[55]; a_Ca_y=i[56]; a_Ca_z=i[57]

    b_Ca_x=i[61]; b_Ca_y=i[62]; b_Ca_z=i[63]

    c_Ca_x=i[67]; c_Ca_y=i[68]; c_Ca_z=i[69]
    
    ab_type=i[76]; ac_type=i[77]; bc_type=i[78]

    a_CG=[0,0,0]; b_CG=[0,0,0]; c_CG=[0,0,0]

    f_zinc = [0,0,0]
    if resname_c == 'x':
        f_ab = cys_his_2.ZincCoord2(a1_ch,b1_ch,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_Ca_x,a_Ca_y,a_Ca_z,b_Ca_x,b_Ca_y,b_Ca_z).f_zinc_result()
        f_zinc = f_ab[0]
        a_atom = f_ab[1]
        b_atom = f_ab[2]
        c_atom = 'x'
        atom_coord_a = f_ab[3]; atom_coord_b = f_ab[4]
        dist_af = dist(f_zinc,atom_coord_a); dist_bf = dist(f_zinc,atom_coord_b); dist_cf =-9999
    else:
        a = get_atom(resname_a,a1.coord,a_cd2,a_ne2,a_ce1,a_nd1)
        b = get_atom(resname_b,b1.coord,b_cd2,b_ne2,b_ce1,b_nd1)
        c = get_atom(resname_c,c1.coord,c_cd2,c_ne2,c_ce1,c_nd1)
        abc = AtomSelect(a,b,c).get_best_atom()
        a_coord=abc[0]; b_coord=abc[1]; c_coord=abc[2]; atom_a=abc[3]; atom_b=abc[4]; atom_c=abc[5]
        f_ab = ZincCoords(ab_type,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_Cb_x,a_Cb_y,a_Cb_z,b_Cb_x,b_Cb_y,b_Cb_z,a1,b1,a1_ch,b1_ch,resname_b,a_Ca_x,a_Ca_y,a_Ca_z,b_Ca_x,b_Ca_y,b_Ca_z,a_coord,atom_a,b_coord,atom_b).zinc_coore()
        f_ab_coord = f_ab[0]
        a_atom = f_ab[1]; b_atom=f_ab[2]
        f_ac = ZincCoords(ac_type,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,a_Cb_x,a_Cb_y,a_Cb_z,c_Cb_x,c_Cb_y,c_Cb_z,a1,c1,a1_ch,c1_ch,resname_c,a_Ca_x,a_Ca_y,a_Ca_z,c_Ca_x,c_Ca_y,c_Ca_z,a_coord,atom_a,c_coord,atom_c).zinc_coore()
        f_ac_coord = f_ac[0]
        c_atom = f_ac[2]
        f_bc = ZincCoords(bc_type,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,b_Cb_x,b_Cb_y,b_Cb_z,c_Cb_x,c_Cb_y,c_Cb_z,b1,c1,b1_ch,c1_ch,resname_c,b_Ca_x,b_Ca_y,b_Ca_z,c_Ca_x,c_Ca_y,c_Ca_z,b_coord,atom_b,c_coord,atom_c).zinc_coore()
        f_bc_coord = f_bc[0]

        d_f_ab_ac=dist(f_ab_coord,f_ac_coord)
        d_f_ab_bc=dist(f_ab_coord,f_bc_coord)
        d_f_ac_bc=dist(f_ac_coord,f_bc_coord)
        zinc_list = [(f_ab_coord,f_ac_coord),(f_ab_coord,f_bc_coord),(f_ac_coord,f_bc_coord)]
        dst_d = {}
        for i in range(len(zinc_list)):
            dst_d[i] = dist(zinc_list[i][0],zinc_list[i][1])
        d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
        min_dist = d_order[0][0]
        a = zinc_list[min_dist][0]
        b = zinc_list[min_dist][1]
        f_zinc[0] =(a[0]+b[0])/2
        f_zinc[1] =(a[1]+b[1])/2
        f_zinc[2] =(a[2]+b[2])/2
        dist_af = dist(f_zinc,a_coord); dist_bf = dist(f_zinc,b_coord); dist_cf = dist(f_zinc,c_coord)

    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%8.2f,%8.2f,%8.2f" %(pid,f_zinc[0],f_zinc[1],f_zinc[2],a_atom,b_atom,c_atom,dist_af,dist_bf,dist_cf),file=print_log)

print_log.close()




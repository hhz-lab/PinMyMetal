from math import *
from Bio.PDB import *
import numpy as np
import os, sys
import copy


def dist(x,y):
    return sqrt((x.coord[0]-y.coord[0])**2 + (x.coord[1]-y.coord[1])**2 + (x.coord[2]-y.coord[2])**2)

def dist2(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def angle_cd(x,y,z):
    a = dist(y,z)
    b = dist(x,z)
    c = dist(x,y)
    return acos((c*c+a*a-b*b)/(2*a*c))

def get_matrix(c1, d1, e1):
    c2 = copy.deepcopy(c1)
    d2 = copy.deepcopy(d1)
    e2 = copy.deepcopy(e1)

    c2.coord[0]=0
    c2.coord[1]=-dist(c1,e1)
    c2.coord[2]=0

    d2.coord[0]=-dist(d1,e1)*sin(angle_cd(d1,e1,c1))
    d2.coord[1]=-dist(d1,e1)*cos(angle_cd(d1,e1,c1))
    d2.coord[2]=0

    e2.coord[0]=0
    e2.coord[1]=0
    e2.coord[2]=0

    fixed  = [c2,d2,e2]
    moving = [c1,d1,e1]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup.rotran

def rotran(a1,matrix) :
    matrix_r = np.mat(matrix[0])
    matrix_t = matrix[1]
    matrix_t = np.ravel(matrix_t).tolist()
    t1=matrix_t[0]
    t2=matrix_t[1]
    t3=matrix_t[2]

    a2 = (a1.coord)*matrix_r
    a2_list = np.concatenate(a2).ravel().tolist()
    a2 = a2_list[0]
    
    a2[0] =a2[0] + t1
    a2[1] =a2[1] + t2
    a2[2] =a2[2] + t3
    a2 = [round(i,3) for i in a2]
    return a2


def get_r_matrix(c1, d1, e1,matrix):
    c2=Atom.Atom('CA', [0,0,0], 20, 1, None, 'C1', 1, element='C')
    d2=Atom.Atom('CB', [0,0,0], 20, 1, None, 'C2', 1, element='C')
    e2=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
    c2.coord = rotran(c1,matrix)
    d2.coord = rotran(d1,matrix)
    moving  = [c2,d2,e2]
    fixed = [c1,d1,e1]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup.rotran

# the angle of x_0_a2
def zinc_f2(a2,c2):
    a = abs(a2[0])
    b = abs(a2[1])
    if a==0.0:
        a=0.01
    else:
        a = abs(a2[0])

    m = atan(b/a)
    angle = (90 -(m*180/pi))*(pi/180)
    
    f2=[0,0,0]
    f2_1=[1.2*cos(angle),1.2*sin(angle),0]
    f2_2=[-1.2*cos(angle),-1.2*sin(angle),0]

    if dist2(f2_1,c2) > dist2(f2_2,c2):
        f2[0]=1.2*cos(angle)
        f2[1]=1.2*sin(angle)
    else:
        f2[0]=-1.2*cos(angle)
        f2[1]=-1.2*sin(angle)
    return f2

def r_rotran(a2,r_matrix) :
    r_matrix_r = np.mat(r_matrix[0])
    r_matrix_t = r_matrix[1]
    a1 = a2*r_matrix_r
    a1_list = np.concatenate(a1).ravel().tolist()
    a1 = a1_list[0]

    a1[0] =a1[0] + r_matrix_t[0]
    a1[1] =a1[1] + r_matrix_t[1]
    a1[2] =a1[2] + r_matrix_t[2]
    a1 = [round(i,3) for i in a1]
    return a1


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



def f1_zinc(f1,f_H):

    f_zinc=[0,0,0]
    f_zinc[0] = (f1[0] + f_H[0])/2
    f_zinc[1] = (f1[1] + f_H[1])/2
    f_zinc[2] = (f1[2] + f_H[2])/2
    return f_zinc

class ZincCoord2:
    def __init__(self,a1,b1,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17):
        self.a1=a1; self.b1=b1
        self.s12=s12; self.s13=s13; self.s14=s14; self.s15=s15; self.s16=s16; self.s17=s17
        self.s1 = s1
        self.s2=s2; self.s3=s3; self.s4=s4; self.s5=s5; self.s6=s6
        self.s7=s7; self.s8=s8; self.s9=s9; self.s10=s10; self.s11=s11

    def f_zinc_result(self):
        a1=self.a1; b1=self.b1
        a_Ca_x=self.s12; a_Ca_y=self.s13; a_Ca_z=self.s14; b_Ca_x=self.s15; b_Ca_y=self.s16; b_Ca_z=self.s17
        b_atom =self. s1
        a_cg=self.s2; a_cd2=self.s3; a_ne2=self.s4; a_ce1=self.s5; a_nd1=self.s6
        b_cg=self.s7; b_cd2=self.s8; b_ne2=self.s9; b_ce1=self.s10; b_nd1=self.s11

        c1=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
        d1=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')

        f_H=[0,0,0]

        if b_atom == 'CYS':
            a_CG = [0,0,0]
            a_CG[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3)
            a_CG[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3)
            a_CG[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)
            d_ne2_s = dist2(a_ne2,b1.coord)
            d_ce1_s = dist2(a_ce1,b1.coord)
            d_nd1_s = dist2(a_nd1,b1.coord)
            d_cd2_s = dist2(a_cd2,b1.coord)

            dst_d = {'a_-NE2':d_ne2_s,'a_-CE1':d_ce1_s,'a_-ND1':d_nd1_s,'a_-CD2':d_cd2_s}
            d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
            if d_order[0][0] == 'a_-NE2':
                min_dist_atom = a_ne2
            elif d_order[0][0] == 'a_-CE1':
                min_dist_atom = a_ce1
            elif d_order[0][0] == 'a_-ND1':
                min_dist_atom = a_nd1
            elif d_order[0][0] == 'a_-CD2':
                min_dist_atom = a_cd2
            a1.coord[0] = min_dist_atom[0]
            a1.coord[1] = min_dist_atom[1]
            a1.coord[2] = min_dist_atom[2]

            c1.coord[0] = round((a_Ca_x+a1.coord[0])/2,3)
            c1.coord[1] = round((a_Ca_y+a1.coord[1])/2,3)
            c1.coord[2] = round((a_Ca_z+a1.coord[2])/2,3)
            d1.coord[0] = round((b_Ca_x+b1.coord[0])/2,3)
            d1.coord[1] = round((b_Ca_y+b1.coord[1])/2,3)
            d1.coord[2] = round((b_Ca_z+b1.coord[2])/2,3)

            e2=Atom.Atom('CE2', [0,0,0], 20, 1, None, 'C5', 1, element='C')

            e1=Atom.Atom('CE',[0,0,0], 20, 1, None, 'C1', 1, element='C')
            e1.coord[0] = (a1.coord[0]+b1.coord[0])/2
            e1.coord[1] = (a1.coord[1]+b1.coord[1])/2
            e1.coord[2] = (a1.coord[2]+b1.coord[2])/2
            matrix=get_matrix(c1, d1, e1)

            a2 = rotran(a1,matrix)
            c2 = rotran(c1,matrix)
            r_matrix = get_r_matrix(c1, d1, e1, matrix)
            f2 = zinc_f2(a2,c2)
            f1 = r_rotran(f2,r_matrix)
            
            f_H = f_coord(a_CG,min_dist_atom)
            atom_b = '-SG-'
            atom_a = d_order[0][0].split('_')[1]
            atom_coord_a = a1.coord; atom_coord_b = b1.coord

        else:
            b_CG = [0,0,0]
            b_CG[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3)
            b_CG[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3)
            b_CG[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)
            
            d_ne2_s = dist2(b_ne2,a1.coord)
            d_ce1_s = dist2(b_ce1,a1.coord)
            d_nd1_s = dist2(b_nd1,a1.coord)
            d_cd2_s = dist2(b_cd2,a1.coord)
            dst_d = {'b_-NE2':d_ne2_s,'b_-CE1':d_ce1_s,'b_-ND1':d_nd1_s,'b_-CD2':d_cd2_s}
            d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
            if d_order[0][0] == 'b_-NE2':
                min_dist_atom = b_ne2
            elif d_order[0][0] == 'b_-CE1':
                min_dist_atom = b_ce1
            elif d_order[0][0] == 'b_-ND1':
                min_dist_atom = b_nd1
            elif d_order[0][0] == 'b_-CD2':
                min_dist_atom = b_cd2

            b1.coord[0] = min_dist_atom[0]
            b1.coord[1] = min_dist_atom[1]
            b1.coord[2] = min_dist_atom[2]

            c1.coord[0] = round((a_Ca_x+a1.coord[0])/2,3)
            c1.coord[1] = round((a_Ca_y+a1.coord[1])/2,3)
            c1.coord[2] = round((a_Ca_z+a1.coord[2])/2,3)
            d1.coord[0] = round((b_Ca_x+b1.coord[0])/2,3)
            d1.coord[1] = round((b_Ca_y+b1.coord[1])/2,3)
            d1.coord[2] = round((b_Ca_z+b1.coord[2])/2,3)

            e2=Atom.Atom('CE2', [0,0,0], 20, 1, None, 'C5', 1, element='C')

            e1=Atom.Atom('CE',[0,0,0], 20, 1, None, 'C1', 1, element='C')
            e1.coord[0] = (a1.coord[0]+b1.coord[0])/2
            e1.coord[1] = (a1.coord[1]+b1.coord[1])/2
            e1.coord[2] = (a1.coord[2]+b1.coord[2])/2
            matrix=get_matrix(c1, d1, e1)

            a2 = rotran(a1,matrix)
            c2 = rotran(c1,matrix)
            r_matrix = get_r_matrix(c1, d1, e1, matrix)
            f2 = zinc_f2(a2,c2)
            f1 = r_rotran(f2,r_matrix)

            f_H = f_coord(b_CG,min_dist_atom)
            atom_a = '-SG-'
            atom_b = d_order[0][0].split('_')[1]
            atom_coord_a = a1.coord; atom_coord_b = b1.coord

        f_zinc = f1_zinc(f1,f_H)
        return f_zinc,atom_a,atom_b,atom_coord_a,atom_coord_b
    



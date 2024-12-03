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
def metal_f2(a2,c2):
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



def f1_metal(f1,f_H):

    f_metal=[0,0,0]
    f_metal[0] = (f1[0] + f_H[0])/2
    f_metal[1] = (f1[1] + f_H[1])/2
    f_metal[2] = (f1[2] + f_H[2])/2
    return f_metal

class MetalCoord2:
    def __init__(self,a_atom,b_atom,resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_ca,b_ca,a_cb,b_cb):
        self.a_atom=a_atom; self.b_atom=b_atom
        self.resname_a=resname_a; self.resname_b=resname_b
        self.a_cg=a_cg; self.a_cd2=a_cd2; self.a_ne2=a_ne2; self.a_ce1=a_ce1; self.a_nd1=a_nd1
        self.b_cg=b_cg; self.b_cd2=b_cd2; self.b_ne2=b_ne2; self.b_ce1=b_ce1; self.b_nd1=b_nd1
        self.a_ca=a_ca; self.b_ca=b_ca; self.a_cb=a_cb; self.b_cb=b_cb

    def f_metal_result(self):
        a1 = self.a_atom; b1 = self.b_atom
        resname_a=self.resname_a; resname_b=self.resname_b
        a_cg=self.a_cg; a_cd2=self.a_cd2; a_ne2=self.a_ne2; a_ce1=self.a_ce1; a_nd1=self.a_nd1
        b_cg=self.b_cg; b_cd2=self.b_cd2; b_ne2=self.b_ne2; b_ce1=self.b_ce1; b_nd1=self.b_nd1
        a_ca=self.a_ca; b_ca=self.b_ca; a_cb=self.a_cb; b_cb=self.b_cb

        c1=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
        d1=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')

        f_H=[0,0,0]

        if resname_a == 'HIS':
            a_CG = [0,0,0]
            a_CG[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3)
            a_CG[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3)
            a_CG[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)

            #c and d are the midpoints of the ca atom and the binding metal atom
            c1.coord[0] = round((a_ca[0]+a1.coord[0])/2,3)
            c1.coord[1] = round((a_ca[1]+a1.coord[1])/2,3)
            c1.coord[2] = round((a_ca[2]+a1.coord[2])/2,3)
            d1.coord[0] = round((b_ca[0]+b1.coord[0])/2,3)
            d1.coord[1] = round((b_ca[1]+b1.coord[1])/2,3)
            d1.coord[2] = round((b_ca[2]+b1.coord[2])/2,3)

            e2=Atom.Atom('CE2', [0,0,0], 20, 1, None, 'C5', 1, element='C')
            
            #e is the midpoint of a and b
            e1=Atom.Atom('CE',[0,0,0], 20, 1, None, 'C1', 1, element='C')
            e1.coord[0] = (a1.coord[0]+b1.coord[0])/2
            e1.coord[1] = (a1.coord[1]+b1.coord[1])/2
            e1.coord[2] = (a1.coord[2]+b1.coord[2])/2
            
            # Obtain the matrices for the forward transformation and the inverse transformation in turn.
            matrix=get_matrix(c1, d1, e1)
            a2 = rotran(a1,matrix)
            c2 = rotran(c1,matrix)
            r_matrix = get_r_matrix(c1, d1, e1, matrix)
            # Get the coordinates for predicting the metal.
            f2 = metal_f2(a2,c2)
            f1 = r_rotran(f2,r_matrix)
            f_H = f_coord(a_CG,a1.coord)

        else:
            b_CG = [0,0,0]
            b_CG[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3)
            b_CG[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3)
            b_CG[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)
            
            c1.coord[0] = round((a_ca[0]+a1.coord[0])/2,3)
            c1.coord[1] = round((a_ca[1]+a1.coord[1])/2,3)
            c1.coord[2] = round((a_ca[2]+a1.coord[2])/2,3)
            d1.coord[0] = round((b_ca[0]+b1.coord[0])/2,3)
            d1.coord[1] = round((b_ca[1]+b1.coord[1])/2,3)
            d1.coord[2] = round((b_ca[2]+b1.coord[2])/2,3)

            e2=Atom.Atom('CE2', [0,0,0], 20, 1, None, 'C5', 1, element='C')

            e1=Atom.Atom('CE',[0,0,0], 20, 1, None, 'C1', 1, element='C')
            e1.coord[0] = (a1.coord[0]+b1.coord[0])/2
            e1.coord[1] = (a1.coord[1]+b1.coord[1])/2
            e1.coord[2] = (a1.coord[2]+b1.coord[2])/2
            matrix=get_matrix(c1, d1, e1)

            a2 = rotran(a1,matrix)
            c2 = rotran(c1,matrix)
            r_matrix = get_r_matrix(c1, d1, e1, matrix)
            f2 = metal_f2(a2,c2)
            f1 = r_rotran(f2,r_matrix)

            f_H = f_coord(b_CG,b1.coord)

        f_metal = f1_metal(f1,f_H)
        return f_metal
    



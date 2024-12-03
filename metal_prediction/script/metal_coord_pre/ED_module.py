from math import *
from Bio.PDB import *
import numpy as np
import os, sys
import copy
from itertools import product

def dist(x,y):
    return sqrt((x.coord[0]-y.coord[0])**2 + (x.coord[1]-y.coord[1])**2 + (x.coord[2]-y.coord[2])**2)

def dist2(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def get_C_coord(resname_a,a_cg,a_cb):
    if resname_a=='ASP':
        C_coord = a_cb #=CA of cys
    elif resname_a=='GLU':
        C_coord = a_cg #=CA of cys
    else:
        pass
    return C_coord


def get_metal(a,b):
    # when the ange of aod_metal_b_od is 180
    f_metal=[0,0,0]
    f_metal[0] = (a.coord[0] + b.coord[0])/2
    f_metal[1] = (a.coord[1] + b.coord[1])/2
    f_metal[2] = (a.coord[2] + b.coord[2])/2
    return f_metal

def angle_cd(x,y,z):
    a = dist(y,z)
    b = dist(x,z)
    c = dist(x,y)
    return acos((c*c+a*a-b*b)/(2*a*c))

# Obtain the matrix to move the points c1,d1,e1 to the positions c2,d2,e2 
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
    return sup.rotran  # Use the rotran attribute to obtain the rotation and translation matrix.

# Get coordinates through a matrix
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

# Obtain the inverse matrix to move the points c2,d2,e2 to c1,d1,e1
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
    
    h_dist=1.48
    f2=[0,0,0]
    f2_1=[h_dist*cos(angle),h_dist*sin(angle),0]
    f2_2=[-h_dist*cos(angle),-h_dist*sin(angle),0]

    if dist2(f2_1,c2) > dist2(f2_2,c2):
        f2[0]=h_dist*cos(angle)
        f2[1]=h_dist*sin(angle)
    else:
        f2[0]=-h_dist*cos(angle)
        f2[1]=-h_dist*sin(angle)
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

class MetalCoord:
    # when the ange of aod_metal_b_od is 90
    #a,b: od1/od1 atom; c,d: CB/CG atom
    def __init__(self,a,b,a_cx,b_cx):
        self.a1 = a # a_atom 
        self.b1 = b # b_atom
        # a_cb/cg
        self.a_cx = a_cx
        # b_cb/cg
        self.b_cx = b_cx
        
    def f_metal_result(self):
        a1 = self.a1
        b1 = self.b1
        a_cx = self.a_cx
        b_cx = self.b_cx

        c1=Atom.Atom('CC', [0,0,0], 20, 1, None, 'C3', 1, element='C')
        d1=Atom.Atom('CD', [0,0,0], 20, 1, None, 'C4', 1, element='C')
        #c and d are the midpoints of the cb/cg atom and the binding metal atom

        c1.coord[0] = round((a_cx[0]+a1.coord[0])/2,3)
        c1.coord[1] = round((a_cx[1]+a1.coord[1])/2,3)
        c1.coord[2] = round((a_cx[2]+a1.coord[2])/2,3)

        d1.coord[0] = round((b_cx[0]+b1.coord[0])/2,3)
        d1.coord[1] = round((b_cx[1]+b1.coord[1])/2,3)
        d1.coord[2] = round((b_cx[2]+b1.coord[2])/2,3)

        e2=Atom.Atom('CE2', [0,0,0], 20, 1, None, 'C5', 1, element='C')
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
        return f1

def FinalMetal(a_atom,b_atom,resname_a,acb,acg,resname_b,bcb,bcg):
    if dist(a_atom,b_atom) < 4:
        a_cx = get_C_coord(resname_a,acg,acb)
        b_cx = get_C_coord(resname_b,bcg,bcb)
        f_ab = MetalCoord(a_atom,b_atom,a_cx,b_cx).f_metal_result()
    else:
        f_ab = get_metal(a_atom,b_atom)
    return f_ab


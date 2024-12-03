from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy

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
        c1_v0 = self.c1
        d1_v0 = self.d1

        c1=Atom.Atom('CC', [c1_v0[0],c1_v0[1],c1_v0[2]], 20, 1, None, 'C3', 1, element='C')
        d1=Atom.Atom('CD', [d1_v0[0],d1_v0[1],d1_v0[2]], 20, 1, None, 'C4', 1, element='C')

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


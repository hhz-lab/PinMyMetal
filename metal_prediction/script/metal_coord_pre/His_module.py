from math import *
import numpy as np
import psycopg2 as pg
import os, sys

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)


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
    else:  # x2 == x1
        if y2 != y1:
            f[1] = y2 + 2.1 if y2 > y1 else y2 - 2.1
            f[0] = x1
            f[2] = z1
        elif z2 != z1:
            f[2] = z2 + 2.1 if z2 > z1 else z2 - 2.1
            f[0] = x1
            f[1] = y1
    return f

def f1_metal(fa,fb):
    f_metal=[0,0,0]
    f_metal[0] = (fa[0]+fb[0])/2
    f_metal[1] = (fa[1]+fb[1])/2
    f_metal[2] = (fa[2]+fb[2])/2
    return f_metal


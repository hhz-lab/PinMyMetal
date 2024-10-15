from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
import His_module
import Cys_module
import CH_module
from itertools import product
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'CH_metal_coord23.csv')

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()


sql = """select id,
    cg_coords_a, cd2_coords_a, ne2_coords_a, ce1_coords_a, nd1_coords_a, sg_coords_a, ca_coords_a, cb_coords_a, c_coords_a, o_coords_a, 
    cg_coords_b, cd2_coords_b, ne2_coords_b, ce1_coords_b, nd1_coords_b, sg_coords_b, ca_coords_b, cb_coords_b, c_coords_b, o_coords_b, 
    cg_coords_c, cd2_coords_c, ne2_coords_c, ce1_coords_c, nd1_coords_c, sg_coords_c, ca_coords_c, cb_coords_c, c_coords_c, o_coords_c, 
    resname_a,resname_b,resname_c,site_count
    from pre_zncu_coord_sites where site_count < 4 order by id"""


cur.execute(sql)
data = cur.fetchall()

def atom_coord(coords_string):
    if not coords_string:
        return [0,0,0]
    #convert_coords_string_to_list
    coords_list = [float(coord) for coord in coords_string.split(',')]
    return coords_list


def get_atom(a_name,a_sg,a_cd2,a_ne2,a_ce1,a_nd1):
    if a_name == 'CYS':
        list_coord = [a_sg]
        atom_name= ['-SG-']
    else:
        list_coord = [a_cd2,a_ne2,a_ce1,a_nd1]
        atom_name = ['-CD2','-NE2','-CE1','-ND1']
    return list_coord,atom_name


def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

class AtomSelector:
    def __init__(self):
        pass

    def dist(self, coord1, coord2):
        # calculate Euclidean distance between two 3D points
        return sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2 + (coord1[2]-coord2[2])**2)

    def get_min_distance_combination(self, *molecules):
        # Generate all possible combinations of one coordinate from each molecule
        all_coords = [m[0] for m in molecules]  # List of all coordinate lists
        all_names = [m[1] for m in molecules]  # List of all name lists
        all_combinations = list(product(*all_coords))
        all_name_combinations = list(product(*all_names))

        min_distance_sum = float('inf')
        best_combination = None
        best_names = None

        # Iterate over all combinations
        for coord_combo, name_combo in zip(all_combinations, all_name_combinations):
            # Compute the sum of all pairwise distances within the combination
            distance_sum = sum(self.dist(coord_combo[i], coord_combo[j])
                               for i in range(len(coord_combo))
                               for j in range(i + 1, len(coord_combo)))
            if distance_sum < min_distance_sum:
                min_distance_sum = distance_sum
                best_combination = coord_combo
                best_names = name_combo

        return best_combination, best_names


class MetalCoords():
    def __init__(self,resname_a,resname_b,cg,cd,ne,ce,nd,cg2,cd2,ne2,ce2,nd2,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord):
        self.resname_a=resname_a; self.resname_b=resname_b
        self.a_cg=cg; self.a_cd2=cd; self.a_ne2=ne; self.a_ce1=ce; self.a_nd1=nd
        self.b_cg=cg2; self.b_cd2=cd2; self.b_ne2=ne2; self.b_ce1=ce2; self.b_nd1=nd2
        self.a_cb=a_cb; self.b_cb=b_cb;
        self.a_ca=a_ca; self.b_ca=b_ca
        self.a_coord=a_coord; self.b_coord=b_coord

    def metal_coord(self):
        resname_a=self.resname_a; resname_b=self.resname_b
        a_cg=self.a_cg; a_cd2=self.a_cd2; a_ne2=self.a_ne2; a_ce1=self.a_ce1; a_nd1=self.a_nd1
        b_cg=self.b_cg; b_cd2=self.b_cd2; b_ne2=self.b_ne2; b_ce1=self.b_ce1; b_nd1=self.b_nd1
        a_cb=self.a_cb; b_cb=self.b_cb
        a_ca=self.a_ca; b_ca=self.b_ca
        a_coord=self.a_coord; b_coord=self.b_coord

        a_atom=Atom.Atom('OD', [0,0,0], 20, 1, None, 'C1', 1, element='C')
        b_atom=Atom.Atom('OD', [0,0,0], 20, 1, None, 'C1', 1, element='C')
        a_atom.coord=a_coord; b_atom.coord=b_coord


        if resname_a=='HIS' and resname_b=='HIS':
            a_GC=[0,0,0]; b_GC=[0,0,0] # the centroid of the imidazole ring
            a_GC[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3); a_GC[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3); a_GC[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)

            b_GC[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3); b_GC[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3); b_GC[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)

            fa= His_module.f_coord(a_GC,a_coord); fb= His_module.f_coord(b_GC,b_coord)
            f_metal = His_module.f1_metal(fa,fb)

        elif resname_a=='CYS' and resname_b=='CYS':
            f_metal = Cys_module.ZincCoord(a_atom,b_atom,a_cb,b_cb).f1_coord()

        else:
            f_metal = CH_module.MetalCoord2(a_atom,b_atom,resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_ca,b_ca,a_cb,b_cb).f_metal_result()

        return f_metal

# Calculate the angle between vectors AB and CD
def calculate_angle(A, B, C, D):
    # Convert lists or tuples to NumPy arrays
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    D = np.array(D)

    # Calculate vectors AB and CD
    AB = B - A
    CD = D - C

    # Calculate the dot product of vectors AB and CD
    dot_product = np.dot(AB, CD)

    # Calculate the magnitude (norm) of vectors AB and CD
    norm_AB = np.linalg.norm(AB)
    norm_CD = np.linalg.norm(CD)

    # Calculate the angle in radians between the vectors AB and CD
    angle = np.arccos(dot_product / (norm_AB * norm_CD))

    # Convert the angle from radians to degrees
    angle_degrees = np.degrees(angle)

    return angle_degrees

print_log = open(file_path,'w')

for i in data:
    pid=i[0];
    # H_H cg,cd2,ne2,ce1,nd1
    a_cg=atom_coord(i[1]);  a_cd2=atom_coord(i[2]);  a_ne2=atom_coord(i[3]);  a_ce1=atom_coord(i[4]);  a_nd1=atom_coord(i[5])
    a_sg=atom_coord(i[6]);  a_ca=atom_coord(i[7]);   a_cb=atom_coord(i[8]);   a_c=atom_coord(i[9]);   a_o=atom_coord(i[10])

    b_cg=atom_coord(i[11]); b_cd2=atom_coord(i[12]); b_ne2=atom_coord(i[13]); b_ce1=atom_coord(i[14]); b_nd1=atom_coord(i[15])
    b_sg=atom_coord(i[16]); b_ca=atom_coord(i[17]);  b_cb=atom_coord(i[18]);   b_c=atom_coord(i[19]);   b_o=atom_coord(i[20])

    c_cg=atom_coord(i[21]); c_cd2=atom_coord(i[22]); c_ne2=atom_coord(i[23]); c_ce1=atom_coord(i[24]); c_nd1=atom_coord(i[25])
    c_sg=atom_coord(i[26]); c_ca=atom_coord(i[27]);  c_cb=atom_coord(i[28]);   c_c=atom_coord(i[29]);   c_o=atom_coord(i[30])
    
    resname_a=i[31].strip();resname_b=i[32].strip();resname_c=i[33].strip(); site_count=i[34]

    f_metal = [0,0,0]

    a_info = get_atom(resname_a,a_sg,a_cd2,a_ne2,a_ce1,a_nd1)
    b_info = get_atom(resname_b,b_sg,b_cd2,b_ne2,b_ce1,b_nd1)

    a_name = 'x'; b_name = 'x'; c_name = 'x'
    dist_am = 0; dist_bm = 0; dist_cm = 0
    selector = AtomSelector()
    if site_count == 2:
        ab_coords, ab_names = selector.get_min_distance_combination(a_info, b_info)
        a_coord=ab_coords[0]; b_coord=ab_coords[1]; a_name=ab_names[0]; b_name=ab_names[1]
        f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord).metal_coord()

        f_metal = f_ab
        dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal)
    elif site_count ==3:
        c_info = get_atom(resname_c,c_sg,c_cd2,c_ne2,c_ce1,c_nd1)
        abc_coords,abc_names = selector.get_min_distance_combination(a_info, b_info, c_info) # min_abc
        a_coord=abc_coords[0]; b_coord=abc_coords[1]; c_coord=abc_coords[2]
        a_name=abc_names[0]; b_name=abc_names[1]; c_name=abc_names[2]
        f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord).metal_coord()
        f_ac = MetalCoords(resname_a,resname_c,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,a_cb,c_cb,a_ca,c_ca,a_coord,c_coord).metal_coord()
        f_bc = MetalCoords(resname_b,resname_c,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,b_cb,c_cb,b_ca,c_ca,b_coord,c_coord).metal_coord()

        if resname_a=='HIS' and resname_b=='HIS' and resname_c=='HIS':
            f_metal[0] =(f_ab[0]+f_ac[0]+f_bc[0])/3
            f_metal[1] =(f_ab[1]+f_ac[1]+f_bc[1])/3
            f_metal[2] =(f_ab[2]+f_ac[2]+f_bc[2])/3
            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal)

        else:
            metal_list = [(f_ab,f_ac),(f_ab,f_bc),(f_ac,f_bc)]
            dst_d = {}
            for i in range(len(metal_list)):
                dst_d[i] = dist(metal_list[i][0],metal_list[i][1])
            d_order=sorted(dst_d.items(),key=lambda x:x[1],reverse=False)
            min_dist = d_order[0][0]
            a = metal_list[min_dist][0]
            b = metal_list[min_dist][1]
            f_metal[0] =(a[0]+b[0])/2
            f_metal[1] =(a[1]+b[1])/2
            f_metal[2] =(a[2]+b[2])/2

            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal)

    print("%s\t%8.3f\t%8.3f\t%8.3f\t%s\t%s\t%s\t%8.3f\t%8.3f\t%8.3f" %(pid,f_metal[0],f_metal[1],f_metal[2],a_name,b_name,c_name,dist_am,dist_bm,dist_cm),file=print_log)

print_log.close()





from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
import His_module
import ED_module
import EDH_module
from itertools import product, combinations
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'EDH_metal_coord2345.csv')
conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()


sql = """select id,
    cg_coords_a, cd2_coords_a, ne2_coords_a, ce1_coords_a, nd1_coords_a, od1_coords_a, od2_coords_a, cd_coords_a, cb_coords_a, ca_coords_a,
    cg_coords_b, cd2_coords_b, ne2_coords_b, ce1_coords_b, nd1_coords_b, od1_coords_b, od2_coords_b, cd_coords_b, cb_coords_b, ca_coords_b,
    cg_coords_c, cd2_coords_c, ne2_coords_c, ce1_coords_c, nd1_coords_c, od1_coords_c, od2_coords_c, cd_coords_c, cb_coords_c, ca_coords_c,
    cg_coords_d, cd2_coords_d, ne2_coords_d, ce1_coords_d, nd1_coords_d, od1_coords_d, od2_coords_d, cd_coords_d, cb_coords_d, ca_coords_d,
    cg_coords_e, cd2_coords_e, ne2_coords_e, ce1_coords_e, nd1_coords_e, od1_coords_e, od2_coords_e, cd_coords_e, cb_coords_e, ca_coords_e,
    resname_a,resname_b,resname_c,resname_d,resname_e,site_count,
    c_coords_a, o_coords_a, c_coords_b, o_coords_b, c_coords_c, o_coords_c
    from pre_metal_coord_sites where site_count<6 order by id"""


cur.execute(sql)
data = cur.fetchall()

def atom_coord(coords_string):
    if not coords_string:
        return [0,0,0]
    #convert_coords_string_to_list
    coords_list = [float(coord) for coord in coords_string.split(',')]
    return coords_list


def get_atom(resname_a,a_od1,a_od2,a_cd2,a_ne2,a_ce1,a_nd1):
    if resname_a == 'ASP':
        list_coord = [a_od1,a_od2]
        atom_name= ['-OD1','-OD2']
    elif resname_a == 'GLU':
        list_coord = [a_od1,a_od2]
        atom_name= ['-OE1','-OE2']
    else:
        list_coord = [a_cd2,a_ne2,a_ce1,a_nd1]
        atom_name = ['-CD2','-NE2','-CE1','-ND1']
    return list_coord,atom_name

def get_C_coord(resname_a,acg,acd):
    if resname_a=='ASP':
        C_coord = acg
    else:
        C_coord = acd
    return C_coord

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
    def __init__(self,resname_a,resname_b,cg,cd,ne,ce,nd,cg2,cd2,ne2,ce2,nd2,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord,a_cd,b_cd):
        self.resname_a=resname_a; self.resname_b=resname_b
        self.a_cg=cg; self.a_cd2=cd; self.a_ne2=ne; self.a_ce1=ce; self.a_nd1=nd
        self.b_cg=cg2; self.b_cd2=cd2; self.b_ne2=ne2; self.b_ce1=ce2; self.b_nd1=nd2
        self.a_cb=a_cb; self.b_cb=b_cb;
        self.a_ca=a_ca; self.b_ca=b_ca
        self.a_coord=a_coord; self.b_coord=b_coord
        self.a_cd=a_cd; self.b_cd=b_cd
    
    def metal_coord(self):
        resname_a=self.resname_a; resname_b=self.resname_b
        a_cg=self.a_cg; a_cd2=self.a_cd2; a_ne2=self.a_ne2; a_ce1=self.a_ce1; a_nd1=self.a_nd1
        b_cg=self.b_cg; b_cd2=self.b_cd2; b_ne2=self.b_ne2; b_ce1=self.b_ce1; b_nd1=self.b_nd1
        a_cb=self.a_cb; b_cb=self.b_cb
        a_ca=self.a_ca; b_ca=self.b_ca
        a_coord=self.a_coord; b_coord=self.b_coord
        a_cd=self.a_cd; b_cd=self.b_cd

        a_atom=Atom.Atom('OD', [0,0,0], 20, 1, None, 'C1', 1, element='C')
        b_atom=Atom.Atom('OD', [0,0,0], 20, 1, None, 'C1', 1, element='C')
        a_atom.coord=a_coord; b_atom.coord=b_coord


        if resname_a=='HIS' and resname_b=='HIS':
            a_GC=[0,0,0]; b_GC=[0,0,0] # the centroid of the imidazole ring
            a_GC[0]=round((a_cg[0]+a_cd2[0]+a_ne2[0]+a_ce1[0]+a_nd1[0])/5,3); a_GC[1]=round((a_cg[1]+a_cd2[1]+a_ne2[1]+a_ce1[1]+a_nd1[1])/5,3); a_GC[2]=round((a_cg[2]+a_cd2[2]+a_ne2[2]+a_ce1[2]+a_nd1[2])/5,3)

            b_GC[0]=round((b_cg[0]+b_cd2[0]+b_ne2[0]+b_ce1[0]+b_nd1[0])/5,3); b_GC[1]=round((b_cg[1]+b_cd2[1]+b_ne2[1]+b_ce1[1]+b_nd1[1])/5,3); b_GC[2]=round((b_cg[2]+b_cd2[2]+b_ne2[2]+b_ce1[2]+b_nd1[2])/5,3)

            fa= His_module.f_coord(a_GC,a_coord); fb= His_module.f_coord(b_GC,b_coord)
            f_metal = His_module.f1_metal(fa,fb)

        elif resname_a in ('ASP','GLU') and resname_b in ('ASP','GLU'):
            f_metal = ED_module.FinalMetal(a_atom,b_atom,resname_a,a_cb,a_cg,resname_b,b_cb,b_cg)

        else:
            f_metal = EDH_module.MetalCoord2(a_atom,b_atom,resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_ca,b_ca,a_cb,b_cb).f_metal_result()

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
    pid=i[0]
    # atom coord 
    a_cg=atom_coord(i[1]); a_cd2=atom_coord(i[2]); a_ne2=atom_coord(i[3]); a_ce1=atom_coord(i[4]); a_nd1=atom_coord(i[5])
    a_od1=atom_coord(i[6]); a_od2=atom_coord(i[7]); a_cd=atom_coord(i[8]); a_cb=atom_coord(i[9]); a_ca=atom_coord(i[10])

    b_cg=atom_coord(i[11]); b_cd2=atom_coord(i[12]); b_ne2=atom_coord(i[13]); b_ce1=atom_coord(i[14]); b_nd1=atom_coord(i[15])
    b_od1=atom_coord(i[16]); b_od2=atom_coord(i[17]); b_cd=atom_coord(i[18]); b_cb=atom_coord(i[19]); b_ca=atom_coord(i[20])

    c_cg=atom_coord(i[21]); c_cd2=atom_coord(i[22]); c_ne2=atom_coord(i[23]); c_ce1=atom_coord(i[24]); c_nd1=atom_coord(i[25])
    c_od1=atom_coord(i[26]); c_od2=atom_coord(i[27]); c_cd=atom_coord(i[28]); c_cb=atom_coord(i[29]); c_ca=atom_coord(i[30])

    d_cg=atom_coord(i[31]); d_cd2=atom_coord(i[32]); d_ne2=atom_coord(i[33]); d_ce1=atom_coord(i[34]); d_nd1=atom_coord(i[35])
    d_od1=atom_coord(i[36]); d_od2=atom_coord(i[37]); d_cd=atom_coord(i[38]); d_cb=atom_coord(i[39]); d_ca=atom_coord(i[40])

    e_cg=atom_coord(i[41]); e_cd2=atom_coord(i[42]); e_ne2=atom_coord(i[43]); e_ce1=atom_coord(i[44]); e_nd1=atom_coord(i[45])
    e_od1=atom_coord(i[46]); e_od2=atom_coord(i[47]); e_cd=atom_coord(i[48]); e_cb=atom_coord(i[49]); e_ca=atom_coord(i[50])

    # AA name
    resname_a=i[51]; resname_b=i[52]; resname_c=i[53]; resname_d=i[54]; resname_e=i[55];
    site_count=i[56]

    a_c=atom_coord(i[57]); a_o=atom_coord(i[58]); b_c=atom_coord(i[59]); b_o=atom_coord(i[60]); c_c=atom_coord(i[61]); c_o=atom_coord(i[62])

    f_metal = [0,0,0]

    try:
        a_info = get_atom(resname_a,a_od1,a_od2,a_cd2,a_ne2,a_ce1,a_nd1)
        b_info = get_atom(resname_b,b_od1,b_od2,b_cd2,b_ne2,b_ce1,b_nd1)
        a_name = 'x'; b_name = 'x'; c_name = 'x'; d_name = 'x'; e_name = 'x'
        dist_am = 0; dist_bm = 0; dist_cm = 0; dist_dm = 0; dist_em = 0
        selector = AtomSelector()
        if site_count == 2:
            ab_coords, ab_names = selector.get_min_distance_combination(a_info, b_info)
            a_coord=ab_coords[0]; b_coord=ab_coords[1]; a_name=ab_names[0]; b_name=ab_names[1]
            f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord,a_cd,b_cd).metal_coord()
    
            f_metal = f_ab
            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal)
        elif site_count ==3:
            c_info = get_atom(resname_c,c_od1,c_od2,c_cd2,c_ne2,c_ce1,c_nd1)
            abc_coords,abc_names = selector.get_min_distance_combination(a_info, b_info, c_info) # min_abc
            a_coord=abc_coords[0]; b_coord=abc_coords[1]; c_coord=abc_coords[2]
            a_name=abc_names[0]; b_name=abc_names[1]; c_name=abc_names[2]
            f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord,a_cd,b_cd).metal_coord()
            f_ac = MetalCoords(resname_a,resname_c,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,a_cb,c_cb,a_ca,c_ca,a_coord,c_coord,a_cd,c_cd).metal_coord()
            f_bc = MetalCoords(resname_b,resname_c,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,b_cb,c_cb,b_ca,c_ca,b_coord,c_coord,b_cd,c_cd).metal_coord()
    
            f_metal[0] =(f_ab[0]+f_ac[0]+f_bc[0])/3
            f_metal[1] =(f_ab[1]+f_ac[1]+f_bc[1])/3
            f_metal[2] =(f_ab[2]+f_ac[2]+f_bc[2])/3
            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal)
            
        elif site_count ==4:
            c_info = get_atom(resname_c,c_od1,c_od2,c_cd2,c_ne2,c_ce1,c_nd1)
            d_info = get_atom(resname_d,d_od1,d_od2,d_cd2,d_ne2,d_ce1,d_nd1)
            abcd_coords,abcd_names = selector.get_min_distance_combination(a_info, b_info, c_info, d_info) # min_abcd
            a_coord=abcd_coords[0]; b_coord=abcd_coords[1]; c_coord=abcd_coords[2]; d_coord=abcd_coords[3]
            a_name=abcd_names[0]; b_name=abcd_names[1]; c_name=abcd_names[2]; d_name=abcd_names[3]
            f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,a_cb,b_cb,a_ca,b_ca,a_coord,b_coord,a_cd,b_cd).metal_coord()
            f_ac = MetalCoords(resname_a,resname_c,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,a_cb,c_cb,a_ca,c_ca,a_coord,c_coord,a_cd,c_cd).metal_coord()
            f_bc = MetalCoords(resname_b,resname_c,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,b_cb,c_cb,b_ca,c_ca,b_coord,c_coord,b_cd,c_cd).metal_coord()
            f_ad = MetalCoords(resname_a,resname_d,a_cg,a_cd2,a_ne2,a_ce1,a_nd1,d_cg,d_cd2,d_ne2,d_ce1,d_nd1,a_cb,d_cb,a_ca,d_ca,a_coord,d_coord,a_cd,d_cd).metal_coord()
            f_bd = MetalCoords(resname_b,resname_d,b_cg,b_cd2,b_ne2,b_ce1,b_nd1,d_cg,d_cd2,d_ne2,d_ce1,d_nd1,b_cb,d_cb,b_ca,d_ca,b_coord,d_coord,b_cd,d_cd).metal_coord()
            f_cd = MetalCoords(resname_c,resname_d,c_cg,c_cd2,c_ne2,c_ce1,c_nd1,d_cg,d_cd2,d_ne2,d_ce1,d_nd1,c_cb,d_cb,c_ca,d_ca,c_coord,d_coord,c_cd,d_cd).metal_coord()
    
            f_metal[0] =(f_ab[0]+f_ac[0]+f_bc[0]+f_ad[0]+f_bd[0]+f_cd[0])/6
            f_metal[1] =(f_ab[1]+f_ac[1]+f_bc[1]+f_ad[1]+f_bd[1]+f_cd[1])/6
            f_metal[2] =(f_ab[2]+f_ac[2]+f_bc[2]+f_ad[2]+f_bd[2]+f_cd[2])/6
            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal); dist_dm = dist(d_coord,f_metal)
        else:
            c_info = get_atom(resname_c,c_od1,c_od2,c_cd2,c_ne2,c_ce1,c_nd1)
            d_info = get_atom(resname_d,d_od1,d_od2,d_cd2,d_ne2,d_ce1,d_nd1)
            e_info = get_atom(resname_e,e_od1,e_od2,e_cd2,e_ne2,e_ce1,e_nd1)
            abcde_coords,abcde_names = selector.get_min_distance_combination(a_info, b_info, c_info, d_info, e_info) # min_abcd
            a_coord=abcde_coords[0]; b_coord=abcde_coords[1]; c_coord=abcde_coords[2]; d_coord=abcde_coords[3]; e_coord=abcde_coords[4]
            a_name=abcde_names[0]; b_name=abcde_names[1]; c_name=abcde_names[2]; d_name=abcde_names[3]; e_name=abcde_names[4]
    
            f_ab = MetalCoords(resname_a,resname_b,a_cg,a_cd2,a_ne2,a_ce1,a_nd1, b_cg,b_cd2,b_ne2,b_ce1,b_nd1, a_cb,b_cb, a_ca,b_ca, a_coord,b_coord, a_cd,b_cd).metal_coord()
            f_ac = MetalCoords(resname_a,resname_c,a_cg,a_cd2,a_ne2,a_ce1,a_nd1, c_cg,c_cd2,c_ne2,c_ce1,c_nd1, a_cb,c_cb, a_ca,c_ca, a_coord,c_coord, a_cd,c_cd).metal_coord()
            f_bc = MetalCoords(resname_b,resname_c,b_cg,b_cd2,b_ne2,b_ce1,b_nd1, c_cg,c_cd2,c_ne2,c_ce1,c_nd1, b_cb,c_cb, b_ca,c_ca, b_coord,c_coord, b_cd,c_cd).metal_coord()
            f_ad = MetalCoords(resname_a,resname_d,a_cg,a_cd2,a_ne2,a_ce1,a_nd1, d_cg,d_cd2,d_ne2,d_ce1,d_nd1, a_cb,d_cb, a_ca,d_ca, a_coord,d_coord, a_cd,d_cd).metal_coord()
            f_bd = MetalCoords(resname_b,resname_d,b_cg,b_cd2,b_ne2,b_ce1,b_nd1, d_cg,d_cd2,d_ne2,d_ce1,d_nd1, b_cb,d_cb, b_ca,d_ca, b_coord,d_coord, b_cd,d_cd).metal_coord()
            f_cd = MetalCoords(resname_c,resname_d,c_cg,c_cd2,c_ne2,c_ce1,c_nd1, d_cg,d_cd2,d_ne2,d_ce1,d_nd1, c_cb,d_cb, c_ca,d_ca, c_coord,d_coord, c_cd,d_cd).metal_coord()
            f_ae = MetalCoords(resname_a,resname_e,a_cg,a_cd2,a_ne2,a_ce1,a_nd1, e_cg,e_cd2,e_ne2,e_ce1,e_nd1, a_cb,e_cb, a_ca,e_ca, a_coord,e_coord, a_cd,e_cd).metal_coord()
            f_be = MetalCoords(resname_b,resname_e,b_cg,b_cd2,b_ne2,b_ce1,b_nd1, e_cg,e_cd2,e_ne2,e_ce1,e_nd1, b_cb,e_cb, b_ca,e_ca, b_coord,e_coord, b_cd,e_cd).metal_coord()
            f_ce = MetalCoords(resname_c,resname_e,c_cg,c_cd2,c_ne2,c_ce1,c_nd1, e_cg,e_cd2,e_ne2,e_ce1,e_nd1, c_cb,e_cb, c_ca,e_ca, c_coord,e_coord, c_cd,e_cd).metal_coord()
            f_de = MetalCoords(resname_d,resname_e,d_cg,d_cd2,d_ne2,d_ce1,d_nd1, e_cg,e_cd2,e_ne2,e_ce1,e_nd1, d_cb,e_cb, d_ca,e_ca, d_coord,e_coord, d_cd,e_cd).metal_coord()
    
            f_metal[0] =(f_ab[0]+f_ac[0]+f_bc[0]+f_ad[0]+f_bd[0]+f_cd[0]+f_ae[0]+f_be[0]+f_ce[0]+f_de[0])/10
            f_metal[1] =(f_ab[1]+f_ac[1]+f_bc[1]+f_ad[1]+f_bd[1]+f_cd[1]+f_ae[1]+f_be[1]+f_ce[1]+f_de[1])/10
            f_metal[2] =(f_ab[2]+f_ac[2]+f_bc[2]+f_ad[2]+f_bd[2]+f_cd[2]+f_ae[2]+f_be[2]+f_ce[2]+f_de[2])/10
            dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal); dist_dm = dist(d_coord,f_metal); dist_em = dist(e_coord,f_metal)

    
        print("%s\t%8.3f\t%8.3f\t%8.3f\t%s\t%s\t%s\t%s\t%s\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f" %(pid,f_metal[0],f_metal[1],f_metal[2],a_name,b_name,c_name,d_name,e_name,dist_am,dist_bm,dist_cm,dist_dm,dist_em),file=print_log)
    except:
        pass

print_log.close()




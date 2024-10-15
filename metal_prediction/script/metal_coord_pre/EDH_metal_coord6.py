from math import *
from Bio.PDB import *
import numpy as np
import psycopg2 as pg
import os, sys
import copy
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

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'EDH_metal_coord6.csv')
conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

sql = """select id,
    cd2_coords_a, ne2_coords_a, ce1_coords_a, nd1_coords_a, od1_coords_a, od2_coords_a,
    cd2_coords_b, ne2_coords_b, ce1_coords_b, nd1_coords_b, od1_coords_b, od2_coords_b,
    cd2_coords_c, ne2_coords_c, ce1_coords_c, nd1_coords_c, od1_coords_c, od2_coords_c,
    cd2_coords_d, ne2_coords_d, ce1_coords_d, nd1_coords_d, od1_coords_d, od2_coords_d,
    cd2_coords_e, ne2_coords_e, ce1_coords_e, nd1_coords_e, od1_coords_e, od2_coords_e,
    cd2_coords_f, ne2_coords_f, ce1_coords_f, nd1_coords_f, od1_coords_f, od2_coords_f,
    resname_a,resname_b,resname_c,resname_d,resname_e,resname_f,site_count
    from pre_metal_coord_sites where site_count = 6 order by id"""


cur.execute(sql)
data = cur.fetchall()

def atom_coord(coords_string):
    if not coords_string:
        return [0,0,0]
    #convert_coords_string_to_list
    coords_list = [float(coord) for coord in coords_string.split(',')]
    return coords_list

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

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

print_log = open(file_path,'w')

for i in data:
    pid=i[0]
    # H_H cg,cd2,ne2,ce1,nd1
    a_cd2=atom_coord(i[1]); a_ne2=atom_coord(i[2]); a_ce1=atom_coord(i[3]); a_nd1=atom_coord(i[4])
    # E_D OD1/OD2
    a_od1=atom_coord(i[5]); a_od2=atom_coord(i[6])

    b_cd2=atom_coord(i[7]); b_ne2=atom_coord(i[8]); b_ce1=atom_coord(i[9]); b_nd1=atom_coord(i[10])
    b_od1=atom_coord(i[11]); b_od2=atom_coord(i[12])
    c_cd2=atom_coord(i[13]); c_ne2=atom_coord(i[14]); c_ce1=atom_coord(i[15]); c_nd1=atom_coord(i[16])
    c_od1=atom_coord(i[17]); c_od2=atom_coord(i[18])
    d_cd2=atom_coord(i[19]); d_ne2=atom_coord(i[20]); d_ce1=atom_coord(i[21]); d_nd1=atom_coord(i[22])
    d_od1=atom_coord(i[23]); d_od2=atom_coord(i[24])
    e_cd2=atom_coord(i[25]); e_ne2=atom_coord(i[26]); e_ce1=atom_coord(i[27]); e_nd1=atom_coord(i[28])
    e_od1=atom_coord(i[29]); e_od2=atom_coord(i[30])
    f_cd2=atom_coord(i[31]); f_ne2=atom_coord(i[32]); f_ce1=atom_coord(i[33]); f_nd1=atom_coord(i[34])
    f_od1=atom_coord(i[35]); f_od2=atom_coord(i[36])

    # AA name
    resname_a=i[37]; resname_b=i[38]; resname_c=i[39]; resname_d=i[40]; resname_e=i[41]; resname_f=i[42]
    site_count=i[43]

    a_info = get_atom(resname_a,a_od1,a_od2,a_cd2,a_ne2,a_ce1,a_nd1)
    b_info = get_atom(resname_b,b_od1,b_od2,b_cd2,b_ne2,b_ce1,b_nd1)
    c_info = get_atom(resname_c,c_od1,c_od2,c_cd2,c_ne2,c_ce1,c_nd1)
    d_info = get_atom(resname_d,d_od1,d_od2,d_cd2,d_ne2,d_ce1,d_nd1)
    e_info = get_atom(resname_e,e_od1,e_od2,e_cd2,e_ne2,e_ce1,e_nd1)
    f_info = get_atom(resname_f,f_od1,f_od2,f_cd2,f_ne2,f_ce1,f_nd1)
    selector = AtomSelector()
    abcdef_coords,abcdef_names = selector.get_min_distance_combination(a_info, b_info, c_info, d_info, e_info, f_info)
    
    a_coord=abcdef_coords[0]; b_coord=abcdef_coords[1]; c_coord=abcdef_coords[2]; d_coord=abcdef_coords[3]; e_coord=abcdef_coords[4]; f_coord=abcdef_coords[5]
    a_name=abcdef_names[0]; b_name=abcdef_names[1]; c_name=abcdef_names[2]; d_name=abcdef_names[3]; e_name=abcdef_names[4]; f_name=abcdef_names[5]

    f_metal = [0,0,0]

    f_metal[0]=(a_coord[0]+b_coord[0]+c_coord[0]+d_coord[0]+e_coord[0]+f_coord[0])/6
    f_metal[1]=(a_coord[1]+b_coord[1]+c_coord[1]+d_coord[1]+e_coord[1]+f_coord[1])/6
    f_metal[2]=(a_coord[2]+b_coord[2]+c_coord[2]+d_coord[2]+e_coord[2]+f_coord[2])/6
    dist_am = dist(a_coord,f_metal); dist_bm = dist(b_coord,f_metal); dist_cm = dist(c_coord,f_metal); dist_dm = dist(d_coord,f_metal); dist_em = dist(e_coord,f_metal); dist_fm = dist(f_coord,f_metal)

    print("%s,%8.3f,%8.3f,%8.3f,%s,%s,%s,%s,%s,%s,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f" %(pid,f_metal[0],f_metal[1],f_metal[2],a_name,b_name,c_name,d_name,e_name,f_name,dist_am,dist_bm,dist_cm,dist_dm,dist_em,dist_fm), file=print_log)

print_log.close()



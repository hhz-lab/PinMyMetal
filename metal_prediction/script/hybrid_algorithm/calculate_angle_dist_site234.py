from math import *
from Bio.PDB import *
import numpy as np
import os, sys
import copy
from itertools import product, combinations
import getopt
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils"))
import db_utils

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'calculate_angle_dist_site234.csv')
print_log = open(file_path,'w')

# Get a database connection
conn = db_utils.create_connection(DBname)
cur = conn.cursor()


sql = """select distinct a.id,metal_x,metal_y,metal_z,
    CA_coords_a,CB_coords_a,O_coords_a,C_coords_a,
    CA_coords_b,CB_coords_b,O_coords_b,C_coords_b,
    CA_coords_c,CB_coords_c,O_coords_c,C_coords_c,
    CA_coords_d,CB_coords_d,O_coords_d,C_coords_d,site_count, 'ch' as sitetype from
    (select distinct id,metal_x,metal_y,metal_z from zncu_predict_sites_2) a left join pre_zncu_coord_sites b
    on a.id=b.id union
    select distinct c.id,metal_x,metal_y,metal_z,
    CA_coords_a,CB_coords_a,O_coords_a,C_coords_a,
    CA_coords_b,CB_coords_b,O_coords_b,C_coords_b,
    CA_coords_c,CB_coords_c,O_coords_c,C_coords_c,
    CA_coords_d,CB_coords_d,O_coords_d,C_coords_d,site_count, 'edh' as sitetype from
    (select distinct id,metal_x,metal_y,metal_z from metal_predict_sites_2) c left join pre_metal_coord_sites d
    on c.id=d.id"""


cur.execute(sql)
data = cur.fetchall()
cur.close()
conn.close()

def atom_coord(coords_string):
    #convert_coords_string_to_list
    if not coords_string:
        return [0,0,0]
    coords_list = [float(coord) for coord in coords_string.split(',')]
    return coords_list

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

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
    
    if norm_AB == 0 or norm_CD == 0:
        # Handle the case where the norm is zero
        return np.nan

    # Calculate the angle in radians between the vectors AB and CD
    angle = np.arccos(dot_product / (norm_AB * norm_CD))
    
    # Convert the angle from radians to degrees
    angle_degrees = np.degrees(angle)
    
    return angle_degrees

for i in data:
    pid=i[0];
    # H_H cg,cd2,ne2,ce1,nd1
    f_metal=[i[1],i[2],i[3]] 
    a_ca=atom_coord(i[4]); a_cb=atom_coord(i[5]); a_o=atom_coord(i[6]); a_c=atom_coord(i[7])
    b_ca=atom_coord(i[8]); b_cb=atom_coord(i[9]); b_o=atom_coord(i[10]); b_c=atom_coord(i[11])
    
    site_count=i[20]; sitetype=i[21]

    as_ca=a_ca
    ab_angle1 = calculate_angle(a_ca,a_cb,as_ca,b_ca) # the angle1 of CA_a to CB_a and CA_a to CA_b
    bs_ca=b_ca
    ab_angle2 = calculate_angle(b_ca,b_cb,bs_ca,a_ca) # the angle2 of CA_B to CB_B and CA_B to CA_a
    ab_angle3 = calculate_angle(a_ca,a_cb,b_ca,b_cb) # the angle3 of CA_a to CB_a and CA_b to CB_b


    ca_dist_ab = dist(a_ca,b_ca) # the distance between their CA atoms
    cb_dist_ab = dist(a_cb,b_cb) # 5the distance between their CB atoms
    
    ma_angle = calculate_angle(a_ca,f_metal,as_ca,a_cb) # the angle between the metal, CA and CB 
    mb_angle = calculate_angle(b_ca,f_metal,bs_ca,b_cb) 

    mCA_dist_a = dist(f_metal,a_ca); mCA_dist_b = dist(f_metal,b_ca) # the dist from metal to CA
    mCB_dist_a = dist(f_metal,a_cb); mCB_dist_b = dist(f_metal,b_cb) # 4the dist from metal to CB

    mo_dist_a = dist(f_metal,a_o); mo_dist_b = dist(f_metal,b_o) # the dist from metal to the backbone oxygen
    as_o=a_o
    mo_angle_a = calculate_angle(a_o,f_metal,as_o,a_c) # the angle between the metal, backbone oxygen and backbone carbon 
    bs_o=b_o
    mo_angle_b = calculate_angle(b_o,f_metal,bs_o,b_c) 

    ### 3
    ac_angle1 = 0; ac_angle2 = 0; ac_angle3= 0
    mc_angle = 0; ca_dist_ac = 0; ca_dist_bc = 0
    mCA_dist_c =0; mo_dist_c=0; mo_angle_c=0;
    bc_angle1=0; bc_angle2=0; bc_angle3=0;
    cb_dist_ac = 0; cb_dist_bc = 0; mCB_dist_c =0 # 3
    ### 4
    ad_angle1 = 0; ad_angle2 = 0; ad_angle3= 0
    md_angle = 0; ca_dist_ad = 0; ca_dist_bd = 0; ca_dist_cd = 0
    mCA_dist_d =0; mo_dist_d=0; mo_angle_d=0;
    bd_angle1=0; bd_angle2=0; bd_angle3=0;
    cb_dist_ad = 0; cb_dist_bd = 0; cb_dist_cd = 0; mCB_dist_d =0
    cd_angle1=0; cd_angle2=0; cd_angle3=0;

    if site_count== 2:
        pass
    else:
        c_ca=atom_coord(i[12]); c_cb=atom_coord(i[13]); c_o=atom_coord(i[14]); c_c=atom_coord(i[15])
        # ac
        ac_angle1 = calculate_angle(a_ca,a_cb,as_ca,c_ca) # the angle1 of CA_a to CB_a and CA_a to CA_c
        cs_ca=c_ca
        ac_angle2 = calculate_angle(c_ca,c_cb,cs_ca,a_ca) # the angle2 of CA_B to CB_B and CA_B to CA_a
        ac_angle3 = calculate_angle(a_ca,a_cb,c_ca,c_cb) # the angle3 of CA_a to CB_a and CA_b to CB_b

        mc_angle = calculate_angle(c_ca,f_metal,cs_ca,c_cb)

        ca_dist_ac = dist(a_ca,c_ca); ca_dist_bc = dist(b_ca,c_ca)
        cb_dist_ac = dist(a_cb,c_cb); cb_dist_bc = dist(b_cb,c_cb) # 2

        # new
        mCA_dist_c = dist(f_metal,c_ca) #
        mCB_dist_c = dist(f_metal,c_cb) # 1
        mo_dist_c = dist(f_metal,c_o) #
        cs_o=c_o
        mo_angle_c = calculate_angle(c_o,f_metal,cs_o,c_c) #
        
        # bc
        bc_angle1 = calculate_angle(b_ca,b_cb,bs_ca,c_ca) 
        bc_angle2 = calculate_angle(c_ca,c_cb,cs_ca,b_ca) 
        bc_angle3 = calculate_angle(b_ca,b_cb,c_ca,c_cb) 
        if site_count==3:
            pass
        else:
            d_ca=atom_coord(i[16]); d_cb=atom_coord(i[17]); d_o=atom_coord(i[18]); d_c=atom_coord(i[19])
    
        # ad
            ad_angle1 = calculate_angle(a_ca,a_cb,as_ca,d_ca) 
            ds_ca=d_ca
            ad_angle2 = calculate_angle(d_ca,d_cb,ds_ca,a_ca) 
            ad_angle3 = calculate_angle(a_ca,a_cb,d_ca,d_cb) 
    
            md_angle = calculate_angle(d_ca,f_metal,ds_ca,d_cb)
    
            ca_dist_ad = dist(a_ca,d_ca)
            cb_dist_ad = dist(a_cb,d_cb)
    
            # new
            mCA_dist_d = dist(f_metal,d_ca) #
            mCB_dist_d = dist(f_metal,d_cb) # 1
            mo_dist_d = dist(f_metal,d_o) #
            ds_o=d_o
            mo_angle_d = calculate_angle(d_o,f_metal,ds_o,d_c) #
    
            # bd
            bd_angle1 = calculate_angle(b_ca,b_cb,bs_ca,d_ca)
            bd_angle2 = calculate_angle(d_ca,d_cb,ds_ca,b_ca)
            bd_angle3 = calculate_angle(b_ca,b_cb,d_ca,d_cb)
            ca_dist_bd = dist(b_ca,d_ca)
            cb_dist_bd = dist(b_cb,d_cb)

            # cd
            cd_angle1 = calculate_angle(c_ca,c_cb,cs_ca,d_ca)
            cd_angle2 = calculate_angle(d_ca,d_cb,ds_ca,c_ca)
            cd_angle3 = calculate_angle(c_ca,c_cb,d_ca,d_cb)
            ca_dist_cd = dist(c_ca,d_ca)
            cb_dist_cd = dist(c_cb,d_cb)
    

    values = [
            pid, sitetype, ab_angle1, ab_angle2, ab_angle3, ac_angle1, ac_angle2, ac_angle3, ad_angle1, ad_angle2, ad_angle3,
            bc_angle1, bc_angle2, bc_angle3, bd_angle1, bd_angle2, bd_angle3, cd_angle1, cd_angle2, cd_angle3,
            ma_angle, mb_angle, mc_angle, md_angle,
            mo_angle_a, mo_angle_b,mo_angle_c,mo_angle_d,
            ca_dist_ab, ca_dist_ac, ca_dist_bc,ca_dist_ad, ca_dist_bd, ca_dist_cd,
            cb_dist_ab,cb_dist_ac, cb_dist_bc, cb_dist_ad, cb_dist_bd, cb_dist_cd,
            mCA_dist_a, mCA_dist_b, mCA_dist_c, mCA_dist_d,
            mCB_dist_a, mCB_dist_b, mCB_dist_c, mCB_dist_d,
            mo_dist_a, mo_dist_b, mo_dist_c, mo_dist_d
            ]

    print("\t".join(map(str, values)), file=print_log)




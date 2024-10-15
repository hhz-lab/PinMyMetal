#!/usr/bin/python3
import os, sys
import getopt
import gzip
import psycopg2 as pg
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

data_dir = os.path.join(script_directory, '../', f"{pdbid}_nb_result")

file_path = os.path.join(data_dir, 'pre_metals_info.csv')
print_log = open(file_path,'w')

pdbfile = os.path.join(data_dir, f"pdb{pdbid}.ent")
output_dir = os.path.join(data_dir, 'premetal')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#get ion_bindingsites features for pre_site (CH+EDH) 

sql = """
select distinct pdbid,metal_x,metal_y,metal_z, chainid_ion, pre_id, sitetype from pre_sites
"""

cur.execute(sql)
data = cur.fetchall()

def metal_pdbfile(pdbid):
    list_r1=[]
    list_r2=[]
    resi_max=0
    resq_max=0

    outputfile = output_dir + '/' +"pdb" +pdbid+'.ent'
    if os.path.isfile(outputfile):
        pass
    else:
        with open(pdbfile, 'r') as infile, open(outputfile, 'w') as outfile:
            for line in infile:
                if line.startswith('ENDMDL'):
                    break
                if line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END'):
                    continue
                outfile.write(line)

    with open(outputfile,"r") as fd:
        for line in fd:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                residueid = line[6:12]
                chainid = line[21:22]
                resseq  = line[22:26]
                list_r1.append(residueid)
                resi_max= max(list_r1)
                list_r2.append(resseq)
                resq_max= max(list_r2)
        resi_zinc= int(resi_max)
        resq_zinc= int(resq_max)
    return resi_zinc,resq_zinc,outputfile


last_pdbid=''; last_resi=0; last_resq=0; resi_zinc =0; resq_zinc=0
for i in data:
    pdb_code = i[0]
    zinc_coord=(i[1],i[2],i[3]); x, y, z = zinc_coord 
    chainid_a = i[4].strip(); pid = i[5]; sitetype=i[6]

    if last_pdbid==pdb_code:
        resi_zinc = last_resi + 1
        resq_zinc = last_resq + 1
        if resq_zinc > 9999:
            resq_zinc = 1
    else:
        zinc_info = metal_pdbfile(pdb_code)
        outputfile=zinc_info[2]
        resi_zinc = zinc_info[0] + 1
        resq_zinc = zinc_info[1] + 1

    zinc_line = "HETATM%5d ZN    ZN%2s%4d    %8.3f%8.3f%8.3f  1.00 46.74          ZN\n" % (resi_zinc,chainid_a,resq_zinc,x,y,z)
    with open(outputfile, "a+") as f:
        f.write(zinc_line)
    print('%s,%s,%s,%s,%s,%s'%(pid,sitetype,pdb_code,chainid_a,resq_zinc,resi_zinc),file=print_log)
    last_pdbid=pdb_code; last_resi=resi_zinc; last_resq=resq_zinc


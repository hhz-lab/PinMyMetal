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
select distinct pdbid,metal_x,metal_y,metal_z, chainid_ion, pre_id, sitetype,proba_ismetal from pre_sites
"""

cur.execute(sql)
data = cur.fetchall()

def format_metal_line(resi_metal, metal_label, chainid_ion, resseq_metal, metal_coord, proba_ismetal):
    return "HETATM{:>5}{:>5}@{:>3}{:>2}{:>4}@   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f} 00.00{:>12}\n".format(
        resi_metal, metal_label, metal_label, chainid_ion.strip(), resseq_metal,
        metal_coord[0], metal_coord[1], metal_coord[2], proba_ismetal, metal_label
    )


list_r1 = []; list_r2 = []; resi_zinc = 0; resq_zinc = 0
outputfile = output_dir + '/' + "pdb" + pdbid + '.ent'
with open(pdbfile, 'r') as infile, open(outputfile, 'w') as outfile:
    for line in infile:
        if line.startswith('ENDMDL'):
            break
        if line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END'):
            continue

        if line.startswith('ATOM') or line.startswith('HETATM'):
            residueid = line[6:12]
            chainid = line[21:22]
            resseq = line[22:26]

            list_r1.append(residueid)
            list_r2.append(resseq)
            resi_zinc = int(max(list_r1))
            resq_zinc = int(max(list_r2))
        outfile.write(line)


atom_l = []
for i in data:
    pdb_code = i[0]
    zinc_coord=(i[1],i[2],i[3]) 
    chainid_a = i[4].strip(); pid = i[5]; sitetype=i[6]; proba_ismetal=i[7]

    resi_zinc += 1
    resq_zinc += 1
    if resq_zinc > 9999:
        resq_zinc = 1

    metal_name='ZN'
    zinc_line = format_metal_line(resi_zinc, metal_name, chainid_a, resq_zinc, zinc_coord, proba_ismetal)
    atom_l.append(zinc_line)
    print('%s,%s,%s,%s,%s,%s'%(pid,sitetype,pdb_code,chainid_a,resq_zinc,resi_zinc),file=print_log)
    
with open(outputfile, "a+") as f:
    f.write(''.join(atom_l))


#!/usr/bin/python3
import os, sys
import getopt
import gzip
from math import *
import psycopg2 as pg

script_directory = os.path.dirname(os.path.abspath(__file__))

output_dir = os.path.join(script_directory, 'excute/output_data')

output  = 'default'

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['intput=',
                                                                 ])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
    elif opt in ('-o', '--output'):
        output = arg
        output_clstr = arg+'.clstr'

dict_seqlen = {}

pdbid=inputid
DBname="metal"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

sql="""
select distinct pdbid,metal_x,metal_y,metal_z, chainid_a,chainid_b,chainid_c,chainid_d,chainid_e,chainid_f,
resseq_a,resseq_b,resseq_c,resseq_d,resseq_e,resseq_f, 
trim(atom_a,'-'),trim(atom_b,'-'),trim(atom_c,'-'),trim(atom_d,'-'),trim(atom_e,'-'),trim(atom_f,'-'),
resname_a,resname_b,resname_c,resname_d,resname_e,resname_f,
dist_am,dist_bm,dist_cm,dist_dm,dist_em,dist_fm, 
exp_tag,bench,new_preid,proba_ismetal,metal_label,proba_class,chainid_ion,site_count,metal_pdb from pre_sites order by new_preid
"""

cur.execute(sql)
data = cur.fetchall()

sql2="""
select distinct new_preid,bench,chainid,resseq,chainids_lig,resseqs_lig,resnames_lig,metal_label from exp_sites order by bench,chainid,resseq
"""

cur.execute(sql2)
data2 = cur.fetchall()

inputfile = os.path.join(script_directory, f"{pdbid}_nb_result", f"pdb{pdbid}.ent")
outputfile = os.path.join(output_dir, f"{pdbid}_metal.pdb")

lines=[]; list_r1=[]; list_r2=[]; resi_max=0; resq_max=0
f = open(inputfile,'r')
for r in f.readlines():
    if r[0:6] == "CONECT" or r[0:6] == "MASTER" or r[0:6] == "ENDMDL":
        pass
    else:
        lines.append(r)
        if r[0:4] == "ATOM" or r[0:6] == "HETATM":
            residueid = r[6:12]
            chainid = r[21:22]
            resseq  = r[22:26]
            list_r1.append(residueid)
            resi_max= max(list_r1)
            list_r2.append(resseq)
            resq_max= max(list_r2)
    resi_metal=int(resi_max)
    resseq_metal=int(resq_max)
f.close()

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

def  get_exp(a,b,c):
    exp = b+','+a+','+c
    return exp

def format_link(atomname_a, resname_a, chainid_a, resseq_a, metal_label, chainid_ion, resseq_metal, dist_am):
    return "LINK       {:>4}{:>4}{:>2}{:>4}{:>18}{:>4}{:>2}{:>4}     1555   1555  {:>6.2f}\n".format(
        atomname_a, resname_a, chainid_a, resseq_a,
        metal_label, metal_label, chainid_ion.strip(), resseq_metal, dist_am
    )

def format_metal_line(resi_metal, metal_label, chainid_ion, resseq_metal, metal_coord, proba_ismetal):
    return "HETATM{:>5}{:>5}{:>4}{:>2}{:>4}@   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f} 00.00{:>12}\n".format(
        resi_metal, metal_label, metal_label, chainid_ion.strip(), resseq_metal,
        metal_coord[0], metal_coord[1], metal_coord[2], proba_ismetal, metal_label
    )

file_path = os.path.join(output_dir, f"{pdbid}_output.csv")
print_log = open(file_path,'w')

for d in data2:
    pid=d[0]; bench=d[1]
    chainid=d[2]; resseq=d[3]; chainids_lig=d[4]; resseqs_lig=d[5];resnames_lig=d[6]; metal_label=d[7]
    aa = list(map(get_exp,chainids_lig,resseqs_lig,resnames_lig))
    aa = list(set(aa))
    exp_ligands = ';'.join(aa)
    if bench is True:
        print('{"name":"exp","bench":"true","status":"U","pid":"0","chain":"%s","resseq":"%s","ligand":"%s","metal_name":"%s"}'%(chainid.strip(),resseq,exp_ligands,metal_label),file=print_log)
    else:
        if pid==-1:
            print('{"name":"exp","bench":"false","status":"N","pid":"-1","chain":"%s","resseq":"%s","ligand":"%s","metal_name":"%s"}'%(chainid.strip(),resseq,exp_ligands,metal_label),file=print_log)
        else:
            print('{"name":"exp","bench":"false","status":"Y","pid":"%s","chain":"%s","resseq":"%s","ligand":"%s","metal_name":"%s"}'%(pid,chainid.strip(),resseq,exp_ligands,metal_label),file=print_log)

link_l=[]; atom_l=[]
for i in data:
    pdb_code = i[0]; metal_coord = (i[1], i[2], i[3])
    chainid_a = i[4].strip(); chainid_b = i[5].strip(); chainid_c = i[6].strip(); chainid_d = i[7].strip(); chainid_e = i[8].strip(); chainid_f = i[9].strip()
    resseq_a = i[10]; resseq_b = i[11]; resseq_c = i[12]; resseq_d = i[13]; resseq_e = i[14]; resseq_f = i[15]
    atomname_a = i[16]; atomname_b = i[17]; atomname_c = i[18]; atomname_d = i[19]; atomname_e = i[20]; atomname_f = i[21]
    resname_a = i[22]; resname_b = i[23]; resname_c = i[24]; resname_d = i[25]; resname_e = i[26]; resname_f = i[27]
    dist_am = i[28]; dist_bm = i[29]; dist_cm = i[30]; dist_dm = i[31]; dist_em = i[32]; dist_fm = i[33]
    exp_tag = i[34]; bench = i[35]; pid = i[36]; proba_ismetal = i[37]; metal_label = i[38]; proba_class = i[39]; chainid_ion = i[40]
    site_count = i[41]; metal_pdb = i[42]

    ligands = []
    resi_metal += 1
    resseq_metal += 1

    link1 = format_link(atomname_a, resname_a, chainid_a, resseq_a, metal_pdb, chainid_ion, resseq_metal, dist_am)
    link_l.append(link1)
    link2 = format_link(atomname_b, resname_b, chainid_b, resseq_b, metal_pdb, chainid_ion, resseq_metal, dist_bm)
    link_l.append(link2)

    ligands.append("%s,%s,%s" % (resseq_a, chainid_a, resname_a))
    ligands.append("%s,%s,%s" % (resseq_b, chainid_b, resname_b))

    if site_count >= 3:
        ligands.append("%s,%s,%s" % (resseq_c, chainid_c, resname_c))
        link3 = format_link(atomname_c, resname_c, chainid_c, resseq_c, metal_pdb, chainid_ion, resseq_metal, dist_cm)
        link_l.append(link3)

    if site_count >= 4:
        ligands.append("%s,%s,%s" % (resseq_d, chainid_d, resname_d))
        link4 = format_link(atomname_d, resname_d, chainid_d, resseq_d, metal_pdb, chainid_ion, resseq_metal, dist_dm)
        link_l.append(link4)

    if site_count >= 5:
        ligands.append("%s,%s,%s" % (resseq_e, chainid_e, resname_e))
        link5 = format_link(atomname_e, resname_e, chainid_e, resseq_e, metal_pdb, chainid_ion, resseq_metal, dist_em)
        link_l.append(link5)

    if site_count >= 6:
        ligands.append("%s,%s,%s" % (resseq_f, chainid_f, resname_f))
        link6 = format_link(atomname_f, resname_f, chainid_f, resseq_f, metal_pdb, chainid_ion, resseq_metal, dist_fm)
        link_l.append(link6)
        
    # output_file: metal info
    ligand_str = ";".join(ligands)
    print('{"name":"pre","bench":"%s","status":"U","pid":"%s","resseq":"%s","chain":"%s","ligand":"%s","metal_name":"%s","Metal Probability":"%.2f"}' % (bench, pid, resseq_metal, chainid_ion, ligand_str, metal_label,proba_ismetal), file=print_log)

    metal_line = format_metal_line(resi_metal, metal_pdb, chainid_ion, resseq_metal, metal_coord, proba_ismetal)
    atom_l.append(metal_line)

res = [idx for idx in lines if idx.startswith('CRYST1')]
location=lines.index(res[0])
lines[location:location]=link_l
lines[-1:1]=atom_l
s=''.join(lines)
with open(outputfile,'w') as f:
    f.write(s)


#!/usr/bin/python3
import os, sys
import getopt
import atomium
import gzip
from math import *
from atomium import Model, Atom
import psycopg2 as pg

import hydrophobicity
from hydrophobicity import *

data_dir = "/zinc_prediction/script"
output_dir = "/zinc_prediction/script/excute/output_data"

output  = 'default'


options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['intput=',
                                                                'output=',
                                                                 ])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
    elif opt in ('-o', '--output'):
        output = arg
        output_clstr = arg+'.clstr'


dict_seqlen = {}

pid=inputid
DBname="zinc"+str(pid)

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

sql=" select distinct pdbid, zinc_x,zinc_y,zinc_z,\
        chainid_a, chainid_b,chainid_c,chainid_d,\
        resseq_a, resseq_b,resseq_c,resseq_d,\
        residueid_a, residueid_b, \
        trim(atomname_a,'-'),trim(atomname_b,'-'),trim(atomname_c,'-'),trim(atomname_d,'-'),\
        ml_result,exp_tag,resname_a,resname_b,resname_c,resname_d,\
        dist_af,dist_bf,dist_cf,dist_df,bench,id,p_value\
        from zinc_predict_site234_2 where selected_site is true order by id"

cur.execute(sql)
data = cur.fetchall()

sql2="select a.id,a.bench,chainid,resseq,chainids_lig,resseqs_lig,resnames_lig from exp_site234 a left join protein_metal b on a.pdbid=b.pdbid and a.residueid_ion=b.residueid_ion order by bench,id"
cur.execute(sql2)
data2 = cur.fetchall()


inputfile = data_dir + '/' + pid + "_nb_result" +"/pdb"+ pid + ".ent"
outputfile = output_dir + '/' + pid +'_zinc.pdb'


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
    resi_zinc=int(resi_max)
    resq_zinc=int(resq_max)
f.close()

print_log = open(output_dir + '/'+pid+'_output.csv','w')

def dist(x,y):
    return sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)


def  get_exp(a,b,c):
    exp = b+','+a+','+c

    return exp

for d in data2:
    pid=d[0]; bench=d[1]
    chainid=d[2]; resseq=d[3]; chainids_lig=d[4]; resseqs_lig=d[5];resnames_lig=d[6]
    aa = list(map(get_exp,chainids_lig,resseqs_lig,resnames_lig))
    aa = list(set(aa))
    exp = ';'.join(aa)
    if bench is True:
        print('{"name":"exp","bench":"true","status":"U","pid":"0","chain":"%s","resseq":"%s","ligand":"%s"}'%(chainid.strip(),resseq,exp),file=print_log)
    else:
        if pid==-1:
            print('{"name":"exp","bench":"false","status":"N","pid":"-1","chain":"%s","resseq":"%s","ligand":"%s"}'%(chainid.strip(),resseq,exp),file=print_log)
        else:
            print('{"name":"exp","bench":"false","status":"Y","pid":"%s","chain":"%s","resseq":"%s","ligand":"%s"}'%(pid,chainid.strip(),resseq,exp),file=print_log)


link_l=[]; atom_l=[]
for i in data:
    pdb_code = ''
    pdb_code = i[0]

    pdbfile  = inputfile

    pdb1 = atomium.open(pdbfile)
    zinc_coord=(i[1],i[2],i[3])
    x=zinc_coord[0]
    y=zinc_coord[1]
    z=zinc_coord[2]
    chainid_a = i[4].strip(); chainid_b = i[5].strip(); chainid_c = i[6].strip(); chainid_d = i[7].strip()
    resseq_a = i[8]; resseq_b = i[9]; resseq_c = i[10]; resseq_d = i[11]
    residueid_a = i[12]; residueid_b = i[13]
    atomname_a = i[14]; atomname_b = i[15]; atomname_c = i[16]; atomname_d = i[17]

    ml_result = i[18]; exp_site = i[19]
    resname_a=i[20];resname_b=i[21];resname_c=i[22];resname_d=i[23]
    dist_af=i[24];dist_bf=i[25];dist_cf=i[26];dist_df=i[27]
    bench=i[28];pid=i[29]
    p_value=i[30]

    if exp_site == 1 and bench =='true':
        pass
    else:
        resi_zinc += 1
        resq_zinc += 1
        if ml_result == 1:
            link1 = "LINK%12s%4s%2s%4d                ZN    ZN%2s%4d     1555   1555%6.2f\n"%(atomname_a,resname_a,chainid_a,resseq_a,chainid_a,resq_zinc,dist_af)
            link_l.append(link1)
            link2 = "LINK%12s%4s%2s%4d                ZN    ZN%2s%4d     1555   1555%6.2f\n"%(atomname_b,resname_b,chainid_b,resseq_b,chainid_a,resq_zinc,dist_bf)
            link_l.append(link2)
            if atomname_c == 'x':
                print('{"name":"pre","bench":"%s","status":"U","pid":"%s","resseq":"%s","chain":"%s","ligand":"%s,%s,%s;%s,%s,%s"}'%(bench,pid,resq_zinc,chainid_a,resseq_a,chainid_a,resname_a,resseq_b,chainid_b,resname_b),file=print_log)
            else:
                link3 = "LINK%12s%4s%2s%4d                ZN    ZN%2s%4d     1555   1555%6.2f\n"%(atomname_c,resname_c,chainid_c,resseq_c,chainid_a,resq_zinc,dist_cf)
                link_l.append(link3)
                if atomname_d == 'x':
                    print('{"name":"pre","bench":"%s","status":"U","pid":"%s","resseq":"%s","chain":"%s","ligand":"%s,%s,%s;%s,%s,%s;%s,%s,%s"}'%(bench,pid,resq_zinc,chainid_a,resseq_a,chainid_a,resname_a,resseq_b,chainid_b,resname_b,resseq_c,chainid_c,resname_c),file=print_log)
                else:
                    link4 = "LINK%12s%4s%2s%4d                ZN    ZN%2s%4d     1555   1555%6.2f\n"%(atomname_d,resname_d,chainid_d,resseq_d,chainid_a,resq_zinc,dist_df)
                    link_l.append(link4)
                    print('{"name":"pre","bench":"%s","status":"U","pid":"%s","resseq":"%s","chain":"%s","ligand":"%s,%s,%s;%s,%s,%s;%s,%s,%s;%s,%s,%s"}'%(bench,pid,resq_zinc,chainid_a,resseq_a,chainid_a,resname_a,resseq_b,chainid_b,resname_b,resseq_c,chainid_c,resname_c,resseq_d,chainid_d,resname_d),file=print_log)
            zinc_line = "HETATM%5d ZN    ZN%2s%4d@   %8.3f%8.3f%8.3f%6.2f 00.00          ZN\n" % (resi_zinc,chainid_a,resq_zinc,x,y,z,p_value)
            atom_l.append(zinc_line)

res = [idx for idx in lines if idx.startswith('CRYST1')]
location=lines.index(res[0])
lines[location:location]=link_l
lines[-1:1]=atom_l
s=''.join(lines)
with open(outputfile,'w') as f:
    f.write(s)


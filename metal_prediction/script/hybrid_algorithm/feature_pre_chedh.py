import atomium
from atomium import Model, Atom
import psycopg2 as pg
import os, sys
import copy
import hydrophobicity
from hydrophobicity import *
import gzip
import tempfile
import warnings
import numpy as np
import csv
from chemfeatures import *
from Bio.PDB import *    
import warnings
from Bio import BiopythonWarning
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

data_dir = os.path.join(script_directory, '../', f"{pdbid}_nb_result")
conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

pdbfile = os.path.join(data_dir, f"pdb{pdbid}.ent")

warnings.filterwarnings("ignore", category=BiopythonWarning)

#EDH+CH

sql = "select distinct pdbid,id,metal_x,metal_y,metal_z,site_count,conc_comma(chainid_a||'_'||chainid_b||'_'||chainid_c||'_'||chainid_d) as chainids, 'ch' as sitetype from zncu_predict_sites_2\
        group by pdbid,id,metal_x,metal_y,metal_z,site_count union\
        select distinct pdbid,id,metal_x,metal_y,metal_z,site_count,conc_comma(chainid_a||'_'||chainid_b||'_'||chainid_c||'_'||chainid_d) as chainids, 'edh' as sitetype from metal_predict_sites_2\
        group by pdbid,id,metal_x,metal_y,metal_z,site_count order by sitetype,id"

hydro_file = os.path.join(data_dir, 'hydro_pre_chedh.csv')
chem_file =  os.path.join(data_dir, 'chem_pre_chedh.csv')

cur.execute(sql)
data = cur.fetchall()


class ChainSelect(Select):
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_letters

def filter_pdb_chains(input_file_path, chain_ids):
    with open(input_file_path, 'rt') as file:
        parser = PDBParser()
        structure = parser.get_structure("PDB", file)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_file:
        tmp_file_path = tmp_file.name
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_file_path, ChainSelect(chain_ids))
    return tmp_file_path

hydro_log = open(hydro_file,'w')
chem_log = open(chem_file,'w')

for i in data:
    pdb_code = ''
    pdb_code = i[0]; pid = i[1]

    #pdbfile = os.path.join(data_dir, f"pdb{pdb_code}.ent")

    metal_coord=(i[2],i[3],i[4])
    site_count=i[5]; chainids=i[6]; sitetype=i[7]

    unique_chains = set()
    parts = chainids.split('_')
    unique_chains.update([part for part in parts if part != 'NaN' and part != ''])
    unique_chains = list(unique_chains)

    output_pdb_file = filter_pdb_chains(pdbfile, unique_chains)

    pdb1 = atomium.open(output_pdb_file)

    c_value = hydrophobic_contrast_function(pdb1,metal_coord)
     
    c_value_dict = c_value[0]
    solv_dict = c_value[1]
    for (key,value) in zip(c_value_dict.items(), solv_dict.items()):
        c_r = (key,value)
        print(f"{pid},{c_r[0][0]},{c_r[0][1]},{c_r[1][1]},{sitetype}",flush=True,file=hydro_log)
    
    # Atoms within 5A distance of metal ions
    neighbor_atoms = c_value[2]
    neighbor_atoms = [atom for atom in neighbor_atoms]
    # Get the chemical properties of these atoms
    for atom in neighbor_atoms:
        if atom.het.name == 'HOH':
            pass
        else:
            chem_types, distances = get_chemical_property(atom, metal_coord)
            for chem_type, distance in zip(chem_types, distances):
                chemical_property = chem_type
                distance = distance
                print('%s,%s,%s,%s'%(pid,chemical_property,distance,sitetype),file=chem_log)
    os.remove(output_pdb_file)

hydro_log.close()
chem_log.close()

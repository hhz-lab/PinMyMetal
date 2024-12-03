import atomium
import numpy as np

# Define a function to get the chemical property of an atom
def get_chemical_property(atom, metal_coord):
    chem_types = []
    distances = []

    if atom.het.name in {'HIS', 'TRP', 'TYR', 'PHE'}:
        if atom.name in {'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'CD1', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'CZ', 'OH'}:
            chem_type = 1 # aromatic
            distance = np.linalg.norm(np.array(atom.location) - metal_coord)
            chem_types.append(chem_type)
            distances.append(distance)

    if atom.element in {'C'}:
        chem_type = 2 # hydrophobic
        distance = np.linalg.norm(np.array(atom.location) - metal_coord)
        chem_types.append(chem_type)
        distances.append(distance)

    if atom.name in {'N'}:
        chem_type = 3 # hydrogen-bonding donor
        distance = np.linalg.norm(np.array(atom.location) - metal_coord)
        chem_types.append(chem_type)
        distances.append(distance)

    if atom.het.name in {'ASN', 'GLN', 'TRP', 'MSE', 'SER', 'THR', 'MET', 'CYS'}:
        if atom.name in {'ND2', 'NE2', 'NE1', 'SG', 'SE', 'OG', 'OG1'}:
            chem_type = 3 # hydrogen-bonding donor
            distance = np.linalg.norm(np.array(atom.location) - metal_coord)
            chem_types.append(chem_type)
            distances.append(distance)

    if atom.name in {'O'}:
        chem_type = 4 # hydrogen-bonding acceptor
        distance = np.linalg.norm(np.array(atom.location) - metal_coord)
        chem_types.append(chem_type)
        distances.append(distance)

    if atom.het.name in {'ASP', 'GLU', 'HIS', 'SER', 'THR', 'MSE', 'CYS', 'MET'}:
        if atom.name in {'ND2', 'NE2', 'OE1', 'OE2', 'OD1', 'OD2', 'OG', 'OG1', 'SE', 'SG'}:
            chem_type = 4 # hydrogen-bonding acceptor
            distance = np.linalg.norm(np.array(atom.location) - metal_coord)
            chem_types.append(chem_type)
            distances.append(distance)

    if atom.het.name in {'LYS', 'ARG', 'HIS'}:
        if atom.name in {'NZ', 'NH1', 'NH2', 'ND1', 'NE2', 'NE'}:
            chem_type = 5 # positive charge
            distance = np.linalg.norm(np.array(atom.location) - metal_coord)
            chem_types.append(chem_type)
            distances.append(distance)

    if atom.het.name in {'ASP', 'GLU'}:
        if atom.name in {'OD1', 'OD2', 'OE1', 'OE2'}:
            chem_type = 6 # negative charge
            distance = np.linalg.norm(np.array(atom.location) - metal_coord)
            chem_types.append(chem_type)
            distances.append(distance)

    if atom.name in {'ND1', 'NE2', 'SG', 'OE1', 'OE2', 'OD2'}:
        chem_type = 7 # metalbinding
        distance = np.linalg.norm(np.array(atom.location) - metal_coord)
        chem_types.append(chem_type)
        distances.append(distance)

    if not chem_types:  # If no chem_types were assigned, default to 0
        chem_type = 0
        distance = np.linalg.norm(np.array(atom.location) - metal_coord)
        chem_types.append(chem_type)
        distances.append(distance)

    return chem_types, distances



'''
import pandas as pd
import re
import gzip
from Bio.PDB import PDBParser, Selection
import Bio

data_dir = "/data/databases/pdbrepo/rsync/pub/pdb/data/structures/divided/pdb"


def get_occupancy(pdbid, atomid):
    mid_dir  = pdbid[1:3]
    pdbfile  = data_dir + '/' + mid_dir + "/pdb" + pdbid + ".ent.gz"

    with gzip.open(pdbfile, 'rt') as f:
        p = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = p.get_structure(pdbid, f)

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_serial_number() == atomid:
                        return atom.get_occupancy()
    return None
'''

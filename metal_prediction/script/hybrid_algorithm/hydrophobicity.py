import atomium
from atomium import Model, Atom

aa_resname = {'NCY', 'SAM', 'NCD', 'GLN', 'MLY', 'MAA', 'NEP', 'MHS', 'DTR', 'HYP', 'CHG', 'DTY', '5PN', 'CAF', 'SAR', 'SNC', 'CYG', 'SVA', 'DPR', 'DTH', 'HIP', 'GLY', 'ALO', 'DPN', 'CSX', 'PTR', 'GGL', 'AYA', 'TRP', 'ALA', 'NVA', 'DNE', 'SAC', 'DLE', 'PRO', '4FB', 'ASX', 'SEB', 'LYS', 'OMT', 'VAL', 'CSO', 'HIC', 'TYS', 'DLY', 'DCY', 'DAR', 'R1A', 'ASN', 'BMT', 'THR', 'SMC', 'SCY', 'A6Q', 'DGL', 'GLZ', 'CMH', 'SER', 'A6I', 'HIS', '407', 'KCX', '409', 'IIL', 'DSN', 'ARG', 'BAL', 'DAL', 'DGN', 'PCA', 'A64', 'MLZ', 'GLU', 'TYR', 'ASP', 'SEP', 'LEU', 'ALY', 'MVA', 'NLE', 'OCS', 'CAS', 'DAH', 'TPO', 'FME', 'GLX', 'A6N', 'HAR', 'DP9', 'PHI', 'GL3', 'AIB', 'CME', 'M3L', 'MGN', 'BET', 'MET', 'MLE', 'AGM', 'HMR', 'UNK', 'PHE', 'A63', 'MSE', 'CSS', 'LLP', 'ILE', 'IYR', 'DHA', 'DHI', 'SET', 'DSG', 'CCS', 'PFF', 'CYS', 'DAS', 'BPN', 'A6L', 'PAL', '403', 'CGU', 'DVA', 'CSD', 'BCS'}


#def hydrophobic_contrast_function(pdb1,metal_coord,chainid)
def hydrophobic_contrast_function(pdb1,metal_coord):
    location = metal_coord
    model = pdb1.model
    c = hydrophobic_contrast(model, *location, metal=False)
    return c

def hydrophobic_contrast(model, x, y, z, metal=True):
    
    radius_list=[2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7]
    
    sphere = model.atoms_in_sphere((x, y, z), 7, is_metal=metal)
    c_value_dict={}
    solv_dict={}
    
    for index,a in enumerate(sphere):
        if len(sphere) == 0: return 0
        sum_, r2 = 0, 0
        sum_solv=0
        atom_count = 0
        
        a_radius = radius_list[index]
        atom_list = a
        
        for atom in atom_list:
            if atom.het.name in aa_resname:
                distance = atom.distance_to((x, y, z))
                solv = atom_solvation(atom)
                sum_solv += solv
                sum_ += solv * (distance ** 2)
                r2 += (distance ** 2)
                atom_count += 1
            else:
                pass
        if atom_count == 0:
            r2=0
            average_solvation=0
        else:
            r2 /= atom_count
            average_solvation= sum_solv / atom_count
        c_value = sum_ - (atom_count * average_solvation * r2)
        c_value_dict[a_radius]=c_value
        solv_dict[a_radius]=average_solvation
        
        if index == 12:
            a_value_at_index_12 = a

    return c_value_dict, solv_dict, a_value_at_index_12


def atom_solvation(atom):
    """Returns the atomic solvation parameter of an atomium atom. The atomic
    solvation parameters are taken from Yamashita et al (1990).
    :param Atom atom: an atomium atom object.
    :rtype: ``float``"""

    specials = {
     "O": {"GLU": ["OE1", "OE2"], "ASP": ["OD1", "OD2"]},
     "N": {"HIS": ["ND1", "NE2"], "ARG": ["NH1", "NH2"]}
    }
    if atom.element == "C": return 18
    if atom.element == "S": return -5
    if atom.element in specials:
        if atom.charge != 0:
            return -37 if atom.element == "O" else -38
        if atom.het and atom.het.name in specials[atom.element]:
            if atom.name in specials[atom.element][atom.het.name]:
                return -23 if atom.element == "O" else -23.5
        return -9
    return 0


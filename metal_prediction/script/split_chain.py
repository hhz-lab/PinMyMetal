from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

output_path = /var/www/django/zinc/static/split_chain/
filename = "pdb2z9k.ent"
chains = ['A','B']

p = PDBParser(PERMISSIVE=1)
s = p.get_structure('',filename)

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0


outputfile = output_path+filename
with open(outputfile, "a") as t:
    for r in open(filename,'r').readlines():
        if r[0:6]=="HEADER" or r[0:5]=="CRYST" or r[0:5]=="ORIGX" or r[0:5]=="SCALE":
            t.write(r)
    for chain in chains:
        pdb_chain_file = filename.split('.')[0]+'_{}.pdb'.format(chain)
        io = PDBIO()
        io.set_structure(s)
        io.save('{}'.format(pdb_chain_file), ChainSelect(chain))
        for i in open(pdb_chain_file,'r').readlines():
            if i[0:3]=="TER" or i[0:3]=="END":
                pass
            else:
                t.write(i)
    t.write("END")


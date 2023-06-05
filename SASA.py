""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#SASA of NZ atom of lysines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#Import necessary packages/libraries
import pandas as pd
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley

#Get PDB structure
parser = MMCIFParser()
structure = parser.get_structure("2kut", "2kut.cif")

#Removing Hetero Atoms
residue_to_remove = []
chain_to_remove = []

for chain in structure:
    for residue in chain:
        if residue.id[0] != '':
            residue_to_remove.append((chain.id, residue.id))
    if len(chain) == 0:
        chain_to_remove.append(chain.id)

for residue in residue_to_remove:
    structure[residue[0]].detach_child(residue[1])

for chain in chain_to_remove:
    structure.detach_child(chain)

#Compute SASA at level of atom
sr = ShrakeRupley()
sr.compute(structure, level="A")

#Get list of all LYS NZ atoms in structure
nuc_atoms = []

for residue in structure.get_residues():
    if residue.get_resname() == "LYS":
        for atom in residue:
            if atom.id == "NZ":
                nuc_atoms.append(atom)

#Calculating depth of each NZ atom
sasa_atoms = []
for atom in structure.get_atoms():
    if atom.sasa > 0:
        sasa_atoms.append(atom)

depths_nz = []
for atom in nuc_atoms:
    try:
        nuc_sasa = atom.sasa
    except(AttributeError):
        nuc_sasa = pd.NA

    depth = 0

    if pd.isnull(nuc_sasa):
        depth = pd.NA
    if nuc_sasa > 0:
        depth = 0
    if nuc_sasa == 0:
        depth = 9999
        for sasa_atom in sasa_atoms:
            new_depth = atom - sasa_atom
            if new_depth < depth:
                depth = new_depth
    depths_nz.append(depth)

import os
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)

with open('xlms.dat', 'r') as file:
    interactions = file.readlines()

results = []

for interaction in interactions:
    fields = interaction.strip().split()
    protein_A, protein_B = fields[0], fields[1]
    residue_A, residue_B = int(fields[2]), int(fields[3])

    pdb_filename = f"{protein_A}_{protein_B}_1.pdb"

    if os.path.exists(pdb_filename):
        structure = parser.get_structure('PDB', pdb_filename)
        model = structure[0]

        try:
            ca_atom_A = model['A'][residue_A]['CA']
            b_factor_A = ca_atom_A.get_bfactor()
        except KeyError:
            b_factor_A = 'NA'  # print NA if no corresponding residues found in PDB file

        try:
            ca_atom_B = model['B'][residue_B]['CA']
            b_factor_B = ca_atom_B.get_bfactor()
        except KeyError:
            b_factor_B = 'NA'  # print NA if no corresponding residues found in PDB file
    else:
        b_factor_A = 'NA'  # print NA if no corresponding PDB files
        b_factor_B = 'NA'

    results.append(f"{interaction.strip()} {b_factor_A} {b_factor_B}")

with open('output-plddt.txt', 'w') as out_file:
    for result in results:
        out_file.write(result + '\n')


import os
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser()

with open('xlms.dat', 'r') as file:
    interactions = file.readlines()

results = []

for interaction in interactions:
    fields = interaction.strip().split()
    protein_A, protein_B = fields[0], fields[1]
    residue_A, residue_B = int(fields[2]), int(fields[3])

    distances = []

    for x in range(1, 6):
        pdb_filename = f"{protein_A}_{protein_B}_{x}.pdb"

        if os.path.exists(pdb_filename):
            structure = parser.get_structure('PDB', pdb_filename)
            model = structure[0]

            try:
                ca_atom_A = model['A'][residue_A]['CA']
                ca_atom_B = model['B'][residue_B]['CA']
                distance = ca_atom_A - ca_atom_B
            except KeyError:
                distance = 0  # print 0 if no corresponding residues found in PDB file
        else:
            distance = 0  # print 0 if no corresponding PDB file

        distances.append(f"{distance:.2f}")

    result = f"{protein_A} {protein_B} {residue_A} {residue_B} " + " ".join(distances)
    results.append(result)

with open('output-distance.txt', 'w') as out_file:
    for result in results:
        out_file.write(result + '\n')


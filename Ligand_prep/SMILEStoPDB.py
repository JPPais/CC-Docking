#! /usr/bin/python3

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

file_path = '/compounds.txt'  # Replace with the actual path to your text file
try:
    with open(file_path, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()
except FileNotFoundError:
    print(f"File not found: {file_path}")
    exit(1)

# Check if there are entries to convert
if len(lines) <= 1:
    print("No entries to convert.")
    exit(0)

# Extract all entries except the first one (index 0)
smiles_list = [line.split()[2].strip() for line in lines[1:]]
names_list = [line.split()[0].strip() for line in lines[1:]]

# generate mollist
mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

for i in range(len(mol_list)):

# add hydrogens
    mol_h = Chem.AddHs(mol_list[i])
# generate a conformation
    AllChem.EmbedMolecule(mol_h)
# optimize
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters = 1000, nonBondedThresh = 100)

    filename_pdb = "%s.pdb" % names_list[i]
    filename_sdf = "%s.sdf" % names_list[i]
    filename_pdbqt = "%s.pdbqt" % names_list[i]

# write pdb
    Chem.rdmolfiles.MolToPDBFile(mol_h, filename_pdb)

#allow_amide_torsions

# write SDF 
# We strongly advice you against using PDB format for preparing small
# molecules, since it does not contain information about bond connections. Please
# dont forget to always check the protonation state of your molecules before
# docking. Your success can sometimes hang by just an hydrogen atom. ;-)
# https://autodock-vina.readthedocs.io/en/latest/docking_basic.html

    sdf_mol = Chem.SDWriter(filename_sdf)
    sdf_mol.write(mol_h)
    sdf_mol.close()

# write pdbqt using Meeko
# https://github.com/forlilab/Meeko
# https://autodock-vina.readthedocs.io/en/latest/docking_hydrated.html

    #preparator = MoleculePreparation() 
    #preparator.prepare(mol_h)
    #preparator.show_setup()
    #preparator.write_pdbqt_file(filename_pdbqt) 
    #pdbqt_string = PDBQTWriterLegacy.write_string(mol_h) 
    #print(pdbqt_string, end="")

    preparator = MoleculePreparation() 
    preparator.prepare(mol_h)
    #preparator.show_setup()
    preparator.write_pdbqt_file(filename_pdbqt)

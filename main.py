from rdkit import Chem
from rdkit.Chem import Draw

# Example molecule: Aspirin
aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
molecule = Chem.MolFromSmiles(aspirin_smiles)

# Visualize the molecule
Draw.MolToImage(molecule)

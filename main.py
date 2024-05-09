from VF2 import VF2Matcher, molecule_to_graph, print_matched_subgraph, drawCommonMolecule, benchmark_molecules

ethanol = Chem.MolFromSmiles("CCO")
methanol = Chem.MolFromSmiles("CO")
ethanol_graph = molecule_to_graph(ethanol)
methanol_graph = molecule_to_graph(methanol)

matcher = VF2Matcher(ethanol_graph, methanol_graph)

# TODO: Maybe refactor to choose the larger of the two molecules and grab that name for the file?
drawCommonMolecule(matcher, ethanol, "ethanol")
drawCommonMolecule(matcher, methanol, "methanol")

# Time the ones below
# Assuming you have defined ethanol and methanol using RDKit
ibuprofen = Chem.MolFromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
ibuprofen_graph = molecule_to_graph(ibuprofen)
naproxen = Chem.MolFromSmiles("CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O")
naproxen_graph = molecule_to_graph(naproxen)

matcher = VF2Matcher(ibuprofen_graph, naproxen_graph)
match = matcher.is_match()

print("Is there a match:", match)

matcher.print_largest_match()  # Print the largest common subgraph
print_matched_subgraph(matcher)

matched_atoms_ibuprofen_naproxen = list(matcher.largest_match.keys())

# Highlight these atoms in the ethanol molecule
img = Draw.MolToImage(ibuprofen, highlightAtoms=matched_atoms_ibuprofen_naproxen, size=(300, 300))
drawCommonMolecule(matcher, ibuprofen, "ibuprofen")
drawCommonMolecule(matcher, naproxen, "naproxen")

benzene = Chem.MolFromSmiles("C1=CC=CC=C1")
benzene_graph = molecule_to_graph(benzene)
napthalene = Chem.MolFromSmiles("C1=CC=C2C=CC=CC2=C1")
napthalene_graph = molecule_to_graph(napthalene)

indole = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CN2")
indole_graph = molecule_to_graph(indole)
tryptophan = Chem.MolFromSmiles("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
tryptophan_graph = molecule_to_graph(tryptophan)

acetic_acid = Chem.MolFromSmiles("CC(=O)O")
acetic_acid_graph = molecule_to_graph(acetic_acid)
stearic_acid = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCC(=O)O")
stearic_acid_graph = molecule_to_graph(stearic_acid)

# (benzene_graph, "benzene", napthalene_graph, "napthalene"), (indole_graph, "indole", tryptophan_graph, "tryptophan"), (acetic_acid_graph, "acetic_acid", tryptophan_graph, "tryptophan")


molecule_pairs = [(methanol, methanol_graph, "methanol", ethanol, ethanol_graph, "ethanol"), (ibuprofen, ibuprofen_graph, "ibuprofen", naproxen,naproxen_graph, "naproxen"), (benzene, benzene_graph, "benzene", napthalene, napthalene_graph, "napthalene"), (indole, indole_graph, "indole", tryptophan, tryptophan_graph, "tryptophan"), (acetic_acid, acetic_acid_graph, "acetic_acid", stearic_acid, stearic_acid_graph, "stearic_acid")]
benchmark_molecules(molecule_pairs)
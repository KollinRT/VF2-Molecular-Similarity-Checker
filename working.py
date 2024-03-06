from rdkit import Chem
from rdkit.Chem import Draw
from util import timer
import time
import csv


class VF2Matcher:
    def __init__(self, G1, G2):
        self.G1 = G1  # First graph
        self.G2 = G2  # Second graph
        self.core_1 = {}  # Mapping from G1 to G2
        self.core_2 = {}  # Mapping from G2 to G1
        self.largest_match = {}  # Store the largest match found
        # self.name1 = name1  # Name for the first molecule
        # self.name2 = name2  # Name for the second molecule

    def is_match(self):
        # start timing
        tic = time.perf_counter()  # Start timing
        self.match()  # Perform the matching process; no need to get result...
        toc = time.perf_counter()  # End timing
        elapsed_time = toc - tic
        print(f"Total matching time: {elapsed_time:0.4f} seconds") # Print total time
        return len(self.largest_match) > 0, elapsed_time

    def match(self, n=None):
        # Check and update the largest match found so far
        if len(self.core_1) > len(self.largest_match):
            self.largest_match = self.core_1.copy()

        # Checking if M(s) covers all nodes of G2
        # IF M(s) covers all the nodes of G2 THEN OUTPUT M(s)
        if len(self.core_1) == len(self.G1):
            return True

        if n is None:
            try:
                n = next(iter(set(self.G1.keys()) - set(self.core_1.keys())))
            except StopIteration:
                return True
        # FOREACH (n, m)∈P(s)
        for m in self.G2:
            # IF F(s, n, m) THEN
            if m not in self.core_2 and self.semantic_feasibility(n, m):
                # Compute the state s’ obtained by adding (n, m) to M(s) (core_1/2)
                self.core_1[n] = m
                self.core_2[m] = n
                # TODO: look into this print(self.core_1[n], self.core_2[m])
                unmatched_nodes = list(set(self.G1.keys()) - set(self.core_1.keys()))
                # CALL Match(s')
                if not unmatched_nodes or self.match(unmatched_nodes[0]):
                    return True
                # Backtracking portion
                del self.core_1[n]
                del self.core_2[m]
        return False

# TODO: experimenting with other semantic feasibility
    def semantic_feasibility(self, n, m):
        # Check if the number of neighbors is the same
        if len(self.G1[n]['neighbors']) != len(self.G2[m]['neighbors']):
            return False

        # Check if the atoms are of the same type
        if self.G1[n]['atom'] != self.G2[m]['atom']:
            return False

        # If additional checks are needed, like checking bond types or other specific conditions such as isotopes or stereochemistry, add them here

        # If all checks pass, the nodes are semantically feasible
        return True

    def print_largest_match(self):
        # Gather matched atoms and their symbols
        matched_atoms = {self.G1[n]['atom']: self.G2[self.largest_match[n]]['atom'] for n in self.largest_match}
        print("Largest Common Subgraph Matches:")
        for n in self.largest_match:
            g1_atom_symbol = self.G1[n]['atom']
            g2_atom_symbol = self.G2[self.largest_match[n]]['atom']
            print(f"{g1_atom_symbol} (G1) -> {g2_atom_symbol} (G2)")
        print(f"Total Atoms Matched: {len(self.largest_match)}")


def molecule_to_graph(molecule):
    graph = {}
    for atom in molecule.GetAtoms():
        graph[atom.GetIdx()] = {'atom': atom.GetSymbol(), 'neighbors': set()}

    for bond in molecule.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[start]['neighbors'].add(end)
        graph[end]['neighbors'].add(start)

    return graph

def print_matched_subgraph(matcher):
    """
    Prints matched subgraph
    :param matcher:
    :return:
    """
    print("Matched Subgraph:")
    for n, m in matcher.largest_match.items():
        atom_g1 = matcher.G1[n]['atom']
        atom_g2 = matcher.G2[m]['atom']
        # Assuming 'neighbors' are indices; we convert them to atom symbols
        neighbors_g1 = [matcher.G1[neighbor]['atom'] for neighbor in matcher.G1[n]['neighbors'] if neighbor in matcher.largest_match]
        neighbors_g2 = [matcher.G2[neighbor]['atom'] for neighbor in matcher.G2[m]['neighbors'] if neighbor in matcher.largest_match]
        print(f"{atom_g1} (G1) matched with {atom_g2} (G2), Neighbors in G1: {neighbors_g1}, Neighbors in G2: {neighbors_g2}")

# TODO: add the graph imaging portion

# def benchmark_molecules(G1, G2, G3, G4):
#     """
#     Benchmark molecule by adding 1 atom from G3 to G1 and one additional one to G4.
#     Benchmark two molcules +0/+0, +1/+0, +1/+1.
#     :param G1:
#     :param G2:
#     :param G3:
#     :param G4:
#     :return:
#     """
# def benchmark_molecules(molecule_pairs, output_csv="benchmark_results.csv", output_latex="benchmark_results.tex"):
#     """
#     Benchmark multiple pairs of molecule graphs and output the results to a CSV file and a LaTeX table.
#     :param molecule_pairs: List of tuples (G1, name1, G2, name2) representing pairs of molecule graphs to compare and their strings of names.
#     :param output_csv: Filename for the CSV output.
#     :param output_latex: Filename for the LaTeX output.
#     :return: None
#     """
#     benchmark_results = [("Molecule 1 Size", "Molecule 2 Size", "Match Size", "Execution Time (s)")]
#
#     for (G1, name1, G2, name2) in molecule_pairs:
#         matcher = VF2Matcher(G1, G2)
#         _, elapsed_time = matcher.is_match()
#         match_size = len(matcher.largest_match)
#         benchmark_results.append((len(G1), len(G2), match_size, f"{elapsed_time:0.4f}"))
#
#
#     # Output results to CSV
#     with open(output_csv, mode='w', newline='') as csv_file:
#         writer = csv.writer(csv_file)
#         writer.writerows(benchmark_results)
#
#     # Output results to a LaTeX table format
#     with open(output_latex, mode='w') as latex_file:
#         latex_file.write("\\begin{tabular}{cccc}\n")
#         latex_file.write("\\hline\n")
#         latex_file.write("Molecule 1 Size & Molecule 2 Size & Match Size & Execution Time (s) \\\\\n")
#         latex_file.write("\\hline\n")
#         for row in benchmark_results[1:]:  # Skip the header
#             latex_file.write(f"{row[0]} & {row[1]} & {row[2]} & {row[3]} \\\\\n")
#         latex_file.write("\\hline\n")
#         latex_file.write("\\end{tabular}")

# def benchmark_molecules(molecule_pairs, output_csv="benchmark_results.csv", output_latex="benchmark_results.tex"):
#     """
#     Benchmark multiple pairs of molecule graphs and output the results to a CSV file and a LaTeX table.
#     Formats molecule names with their sizes in parentheses.
#     :param molecule_pairs: List of tuples (G1, name1, G2, name2) representing pairs of molecule graphs to compare and their strings of names.
#     :param output_csv: Filename for the CSV output.
#     :param output_latex: Filename for the LaTeX output.
#     :return: None
#     """
#     benchmark_results = [("Molecule 1 (Size)", "Molecule 2 (Size)", "Match Size", "Execution Time (s)")]
#
#     for (mol1, G1, name1, mol2, G2, name2) in molecule_pairs:
#         matcher = VF2Matcher(G1, G2)
#         _, elapsed_time = matcher.is_match()
#         match_size = len(matcher.largest_match)
#         molecule_1_info = f"{name1} ({len(G1)})"
#         molecule_2_info = f"{name2} ({len(G2)})"
#         benchmark_results.append((molecule_1_info, molecule_2_info, match_size, f"{elapsed_time:0.4f}"))
#
#     # Output results to CSV
#     with open(output_csv, mode='w', newline='') as csv_file:
#         writer = csv.writer(csv_file)
#         writer.writerows(benchmark_results)
#
#     # Output results to a LaTeX table format
#     with open(output_latex, mode='w') as latex_file:
#         latex_file.write("\\begin{tabular}{cccc}\n")
#         latex_file.write("\\hline\n")
#         latex_file.write("Molecule 1 (Size) & Molecule 2 (Size) & Match Size & Execution Time (s) \\\\\n")
#         latex_file.write("\\hline\n")
#         for row in benchmark_results[1:]:  # Skip the header
#             latex_file.write(f"{row[0]} & {row[1]} & {row[2]} & {row[3]} \\\\\n")
#         latex_file.write("\\hline\n")
#         latex_file.write("\\end{tabular}")

from rdkit import Chem
from rdkit.Chem import Draw

def drawCommonSubgraphOnLargerMolecule(matcher, mol1, mol2, molecule_name_prefix):
    # Determine the larger molecule based on the number of atoms
    if len(mol1.GetAtoms()) >= len(mol2.GetAtoms()):
        larger_molecule = mol1
        matched_atoms = [mol1.GetAtomWithIdx(i).GetIdx() for i in matcher.largest_match.keys()]
    else:
        larger_molecule = mol2
        # Map matched atoms in G1 to their counterparts in G2 for highlighting
        matched_atoms_in_g2 = [matcher.largest_match[i] for i in matcher.largest_match.keys()]
        matched_atoms = [mol2.GetAtomWithIdx(i).GetIdx() for i in matched_atoms_in_g2]

    molecule_name = f"{molecule_name_prefix}_highlighted.png"
    Draw.MolToFile(larger_molecule, molecule_name, highlightAtoms=matched_atoms, size=(300, 300))

def benchmark_molecules(molecule_pairs, output_csv="benchmark_results.csv", output_latex="benchmark_results.tex"):
    """
    Benchmark multiple pairs of molecule graphs and output the results to a CSV file and a LaTeX table.
    Formats molecule names with their sizes in parentheses. Additionally, draws the common subgraph on the larger molecule for each pair.
    :param molecule_pairs: List of tuples (mol1, G1, name1, mol2, G2, name2) representing RDKit Mol objects and corresponding molecule graphs to compare, along with their names.
    :param output_csv: Filename for the CSV output.
    :param output_latex: Filename for the LaTeX output.
    :return: None
    """
    benchmark_results = [("Molecule 1 (Size)", "Molecule 2 (Size)", "Match Size", "Execution Time (s)")]

    for (mol1, G1, name1, mol2, G2, name2) in molecule_pairs:
        matcher = VF2Matcher(G1, G2)
        _, elapsed_time = matcher.is_match()
        match_size = len(matcher.largest_match)
        molecule_1_info = f"{name1} ({len(G1)})"
        molecule_2_info = f"{name2} ({len(G2)})"
        benchmark_results.append((molecule_1_info, molecule_2_info, match_size, f"{elapsed_time:0.4f}"))

        # Draw the common subgraph on the larger molecule
        image_name_prefix = f"{name1}_vs_{name2}"
        drawCommonSubgraphOnLargerMolecule(matcher, mol1, mol2, image_name_prefix)

    # Output results to CSV
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(benchmark_results)

    # Output results to a LaTeX table format
    with open(output_latex, mode='w') as latex_file:
        latex_file.write("\\begin{tabular}{cccc}\n")
        latex_file.write("\\hline\n")
        latex_file.write("Molecule 1 (Size) & Molecule 2 (Size) & Match Size & Execution Time (s) \\\\\n")
        latex_file.write("\\hline\n")
        for row in benchmark_results[1:]:  # Skip the header
            latex_file.write(f"{row[0]} & {row[1]} & {row[2]} & {row[3]} \\\\\n")
        latex_file.write("\\hline\n")
        latex_file.write("\\end{tabular}")

def highlight_match_in_molecule(molecule, match_indices):
    from rdkit.Chem.Draw import MolsToGridImage
    return MolsToGridImage([molecule], highlightAtomLists=[match_indices])


# def highlight_pharmacophore(mol, match_indices):
#     """
#     Highlight the pharmacophore in a molecule based on matched atom indices.
#
#     Parameters:
#     - mol: RDKit molecule object.
#     - match_indices: List of atom indices in 'mol' that are part of the pharmacophore.
#     """
#     from rdkit.Chem.Draw import rdMolDraw2D
#     from IPython.display import SVG
#
#     # Create a drawer object and set the drawing options
#     drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
#     drawer.drawOptions().addAtomIndices = True  # Optional: Show atom indices
#     rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=match_indices)
#     drawer.FinishDrawing()
#
#     # Display the drawing in Jupyter Notebook
#     svg = drawer.GetDrawingText().replace('svg:', '')
#     display(SVG(svg))

def drawCommonMolecule(matcher, molecule, molecule_name):
    from rdkit import Chem
    from rdkit.Chem import Draw

    # Pass the matcher object for the VF2Matcher
    matched_atoms = list(matcher.largest_match.keys())
    print(matcher.G1)
    print(matcher.G2)
    print(matcher.largest_match)
    print("Matching atoms:", matched_atoms)
    # Choose one of the two molecules for the overlay
    Draw.MolToFile(molecule, f'{molecule_name}.png', highlightAtoms=matched_atoms, size=(300, 300))


ethanol = Chem.MolFromSmiles("CCO")
methanol = Chem.MolFromSmiles("CO")
#%%
ethanol_graph = molecule_to_graph(ethanol)
methanol_graph = molecule_to_graph(methanol)

matcher = VF2Matcher(ethanol_graph, methanol_graph)

# TODO: Maybe refactor to choose the larger of the two molecules and grab that name for the file?
drawCommonMolecule(matcher, ethanol, "ethanol")
drawCommonMolecule(matcher, methanol, "methanol")

# def check_chirality(self, n, m):
#     # Assuming 'chirality' is a property stored in the node attributes
#     chirality_n = self.G1.nodes[n].get('chirality')
#     chirality_m = self.G2.nodes[m].get('chirality')
#     return chirality_n == chirality_m or not chirality_n or not chirality_m
#
# def check_bond_stereochemistry(self, n, m):
#     # This requires checking the spatial arrangement around double bonds
#     # This example assumes you have a way to determine the configuration (e.g., 'cis' or 'trans') of a bond
#     for neighbor in self.G1[n]:
#         if any(self.G1[n][neighbor]['stereo'] != self.G2[m][m_neighbor]['stereo']
#                for m_neighbor in self.G2[m] if self.core_1.get(neighbor) == m_neighbor):
#             return False
#     return True
#
# def check_isotope(self, n, m):
#     isotope_n = self.G1.nodes[n].get('isotope')
#     isotope_m = self.G2.nodes[m].get('isotope')
#     return isotope_n == isotope_m or not isotope_n or not isotope_m
def read_molecules_from_smiles(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    molecules = [Chem.MolFromSmiles(line.split()[0]) for line in lines]
    return molecules



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
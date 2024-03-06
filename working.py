from rdkit import Chem


class VF2Matcher:
    def __init__(self, G1, G2):
        self.G1 = G1  # First graph
        self.G2 = G2  # Second graph
        self.core_1 = {}  # Mapping from G1 to G2
        self.core_2 = {}  # Mapping from G2 to G1
        self.largest_match = {}  # Store the largest match found

    def is_match(self):
        self.match()
        return len(self.largest_match) > 0

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

    def semantic_feasibility(self, n, m):
        # Simplified check: only consider node count match for demonstration
        return len(self.G1[n]['neighbors']) == len(self.G2[m]['neighbors'])  # TODO: DO we need the same neighbors here? I think yeah but the neighbors do not need to be the same numbers just the same amount....

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

def benchmark_molecules(G1, G2, G3, G4):
    """
    Benchmark molecule by adding 1 atom from G3 to G1 and one additional one to G4.
    Benchmark two molcules +0/+0, +1/+0, +1/+1.
    :param G1:
    :param G2:
    :param G3:
    :param G4:
    :return:
    """

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

def drawCommonMolecule(matcher, molecule):
    from rdkit import Chem
    from rdkit.Chem import Draw

    # Pass the matcher object for the VF2Matcher
    matched_atoms = list(matcher.largest_match.keys())

    # Choose one of the two molecules for the overlay
    Draw.MolToFile(molecule, f'{molecule}_common.png', highlightAtoms=matched_atoms, size=(300, 300))


ethanol = Chem.MolFromSmiles("CCO")
methanol = Chem.MolFromSmiles("CO")
#%%
ethanol_graph = molecule_to_graph(ethanol)
methanol_graph = molecule_to_graph(methanol)

matcher = VF2Matcher(ethanol_graph, methanol_graph)

# TODO: Maybe refactor to choose the larger of the two molecules and grab that name for the file?
drawCommonMolecule(matcher, ethanol)
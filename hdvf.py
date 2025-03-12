import random
import json
import copy
from itertools import chain

from complex import Complex, ChainOptimization

class HDVF:
    """Implementation of a HDVF (a labelling of the cells of a complex) and its reduction"""
    def __init__(self, complex: Complex):
        self.cells = complex.simplices
        self.partition = ['c' for c in self.cells] # map : cells -> label ('p', 's', 'c')
        self.dimension = 3
        self.Hc = [{} for _ in range(self.dimension + 1)] # columns of the matrix H
        self.Fr = [{} for _ in range(self.dimension + 1)] # rows    of the matrix F
        self.Gc = [{} for _ in range(self.dimension + 1)] # columns of the matrix G
        self.Dc = [{} for _ in range(self.dimension + 1)] # columns of the matrix D
        for cell in self.cells:
            q = cell.dimension()
            self.Fr[q][cell.index] = set()
            self.Gc[q][cell.index] = set()
            self.Dc[q][cell.index] = set(cell.boundary)

    
    def dim(self, index: int):
        return self.cells[index].dimension()
    
    def nb_cells(self):
        return len(self.cells)
    
    def is_inside(self, index: int, inside=True) -> bool:
        """Tell if the cell (with index `index`) is inside the complex."""
        assert 0 <= index < len(self.cells), "The index is not valid: " + str(index)
        return self.cells[index].inside == inside

    def faces(self, index: int, inside=True) -> bool:
        indices = []
        for t in self.cells[index].boundary:
            if self.is_inside(t, inside):
                indices.append(t)
        return indices

    def cofaces(self, index: int, inside=True):
        indices = []
        for t in self.cells[index].coboundary:
            if self.is_inside(t, inside):
                indices.append(t)
        return indices
    
    def indices(self):
        pass

    def get_type_cell(self, index: int) -> str:
        return self.partition[index]
    
    def set_type_cell(self, index: int, type: str):
        assert type in ['p', 's', 'c']
        self.partition[index] = type
    
    def get_h(self, index: int):
        q = self.dim(index)
        return copy.copy(self.Hc[q][index])
    
    def get_h_star(self, index: int):
        q = self.dim(index)
        chain = set()
        if q-1 >= 0:
            for key, value in self.Hc[q - 1].items():
                if index in value:
                    chain.add(key)
        return chain

    def get_f(self, index: int):
        q = self.dim(index)
        chain = set()
        for key, value in self.Fr[q].items():
            if index in value:
                chain.add(key)
        return chain

    def get_f_star(self, index: int):
        q = self.dim(index)
        chain = copy.copy(self.Fr[q][index])
        chain.add(index)
        return chain

    def get_g(self, index: int):
        q = self.dim(index)
        chain = copy.copy(self.Gc[q][index])
        chain.add(index)
        return chain

    def get_g_star(self, index: int):
        q = self.dim(index)
        chain = set()
        for key, value in self.Gc[q].items():
            if index in value:
                chain.add(key)
        return chain

    def get_d(self, index: int):
        q = self.dim(index)
        return copy.copy(self.Dc[q].get(index))

    def get_d_star(self, index: int):
        q = self.dim(index)
        chain = set()
        if q + 1 <= self.dimension:
            for key, value in self.Dc[q + 1].items():
                if index in value:
                    chain.add(key)
        return chain
    
    def set_perfect_hdvf(self):
        stability = False
        while not stability:
            stability = True
            pair = self.find_pair_to_add()
            if pair is not None:
                print("Adding pair ", pair)
                self.hdvf_a(*pair)
                stability = False

    def find_pair_to_add(self):
        for q in range(1, self.dimension+1):
            for t, chain in self.Dc[q].items():
                if chain:
                    s = random.choice(list(chain))
                    return (s, t)
        return None

    def clear_hdvf(self):
        pass

    def add_cell_to_chain(self, chain, index):
        """
        Add a cell to a chain.
        Since we use Z/2Z coefficients, if the cell is already in the chain, we remove it.
        """
        if index in chain:
            chain.remove(index)
        else:
            chain.add(index)

    def add_chain_to_chain(self, chain, added_chain):
        """Add a chain to a chain."""
        for item in added_chain:
            self.add_cell_to_chain(chain, item)

    def hdvf_a(self, s, t):
        """
        Add a pair of critical cells to the HDVF.
        Return true if it was possible
        """
        # Checking the validity of the operation
        if self.get_type_cell(s) != 'c':
            print('First cell must be critical')
            return False
        if self.get_type_cell(t) != 'c':
            print('Second cell must be critical')
            return False
        if self.dim(s) + 1 != self.dim(t):
            print('Cells must have consecutive dimensions')
            return False
        if s not in self.get_d(t):
            print(f"<d({t}), {s}> == 0: you cannot do this")
            return False
        # print('Adding cells ' + str(s) + ' and ' + str(t) + ' to the HDVF')
        self.set_type_cell(s, 'p')
        self.set_type_cell(t, 's')

        # Get the necessary submatrices
        q = self.dim(s)
        F21 = self.Fr[q][s]
        G12 = self.Gc[q + 1][t]
        D21 = self.get_d_star(s); D21.remove(t)
        D12 = self.get_d(t); D12.remove(s)

        # (1) Update D
        # (1.a) remove d(s:C)
        del self.Dc[q][s]
        # (1.b) update d(C:C)
        for i in D12:
            for j in D21:
                self.add_cell_to_chain(self.Dc[q + 1][j], i)
        # (1.c) remove d(t:C)
        del self.Dc[q + 1][t]
        # (1.d) remove d(C:s)
        for item in D21:
            self.Dc[q + 1][item].remove(s)
        # (1.e) remove d(C:t)
        if q + 2 <= self.dimension:
            for value in self.Dc[q + 2].values():
                value.discard(t)

        # (2) Update H
        # (2.a) update h(P:S)
        for i in G12:
            for j in F21:
                self.add_cell_to_chain(self.Hc[q][j], i)
        # (2.b) add h(s:S)
        self.Hc[q][s] = G12.copy()
        # (2.c) add h(P:t)
        for item in F21:
            self.Hc[q][item].add(t)
        # (2.d) add h(s:t) = 1
        self.Hc[q][s].add(t)

        # (3) Update F
        # (3.a) update f(P:C)
        for i in D12:
            for j in F21:
                self.add_cell_to_chain(self.Fr[q][i], j)
        # (3.b) add f(s:C)
        for item in D12:
            self.Fr[q][item].add(s)
        # (3.c) remove f(P:s)
        del self.Fr[q][s]
        # (3.d) remove f(P:t)
        del self.Fr[q + 1][t]

        # (4) Update G
        # (4.a) delete g(s:C)
        del self.Gc[q][s]
        # (4.b) update g(C:S)
        for i in G12:
            for j in D21:
                self.add_cell_to_chain(self.Gc[q + 1][j], i)
        # (4.c) add g(C:t)
        for item in D21:
            self.Gc[q + 1][item].add(t)
        # (4.d) delete g(t:C)
        del self.Gc[q + 1][t]
        return True
    
    def number_critical_cells(self):
        count = [0 for _ in range(self.dimension + 1)]
        for index, value in enumerate(self.partition):
            if value == 'c':
                count[self.dim(index)] += 1
        print("Number of critical cells: ", count)

    def check_valid(self):
        """
        Check the validity of the data structure.
        I wrote this because I had errors.
        """
        for i, cell in enumerate(self.cells):
            assert cell.index == i
            q = self.dim(i)
            assert self.partition[i] in ('p', 's', 'c')
            if self.partition[i] == 'p':
                assert i in self.Hc[q]
            elif self.partition[i] == 'c':
                assert i in self.Dc[q]
                assert i in self.Gc[q]
                assert i in self.Fr[q]
                assert self.partition[i] == 'c'
        for q, F in enumerate(self.Fr):
            for i, chain in F.items():
                assert self.dim(i) == q
                assert self.partition[i] == 'c'
                for j in chain:
                    assert self.partition[j] == 'p'
        for q, G in enumerate(self.Gc):
            for i, chain in G.items():
                assert self.dim(i) == q
                assert self.partition[i] == 'c'
                for j in chain:
                    assert self.partition[j] == 's'
        for q, H in enumerate(self.Hc):
            for i, chain in H.items():
                assert self.dim(i) == q
                assert self.partition[i] == 'p'
                for j in chain:
                    assert self.partition[j] == 's'





class HomologicalQuartet(HDVF):
    """
    Implementation of the homological quartets.
    We need to compute two perfect relative HDVF and a pairing between the remaining critical cells.
    It is better if the matrix (F D G) is diagonal
    """
    def __init__(self, complex: Complex):
        super().__init__(complex)
        self.pairs = []

    def random_perfect(self):
        """Set random perfect relative HDVF in K and L-K."""
        self._perfect_hdvf(True)
        self._perfect_hdvf(False)

    def _perfect_hdvf(self, inside):
        """Compute a perfect HDVF in K and in L-K."""
        indices = list(range(self.nb_cells()))
        random.shuffle(indices)
        print(f"\nPerfect HDVF inside={inside}...")
        n_cells = sum(1 for t in indices if self.is_inside(t, inside))
        n_paired_cells = 0
        for t in indices:
            if self.partition[t] == 'c' and self.is_inside(t, inside):
                candidates = [s for s in self.get_d(t) if self.is_inside(s, inside)]
                if candidates:
                    s = random.choice(candidates)
                    self.hdvf_a(s, t)
                    n_paired_cells += 2
                    if n_paired_cells % 100 == 0:
                        print(f"Progress: {100*n_paired_cells/n_cells:.2f} % ({n_paired_cells} / {n_cells})", end="\r")
    

    def distance_hdvf(self, vertices, surface_vertices):
        """Make two perfect relative HDVF using the signed distance transform"""
        filtration = {}
        for cell in self.cells:
            if cell.dimension() == 3:
                filtration[cell.index] = (self.filtration_value(cell, vertices, surface_vertices), -3)
        for q in range(2, -1, -1):
            for cell in self.cells:
                if cell.dimension() == q:
                    value = max([filtration[i][0] for i in cell.coboundary])
                    filtration[cell.index] = (value, -q)
        filter = sorted(filtration.keys(), key=lambda k: filtration[k], reverse=False)
        assert len(filter) == len(self.cells)
        n_paired_cells = 0
        print("Filtration computed, building HDVF...")
        for t in filter:
            cell = self.cells[t]
            candidates = self.get_d_star(t)
            if candidates:
                s = max(candidates, key=lambda i: filtration[i])
                if self.cells[t].inside == self.cells[s].inside:
                    self.hdvf_a(t, s)
                    n_paired_cells += 2
                    if n_paired_cells % 100 == 0:
                        print(f"Progress: {100*n_paired_cells/len(filter):.2f} % ({n_paired_cells} / {len(filter)})", end="\r")
                        


    def filtration_value(self, simplex, vertices, surface_vertices):
        assert simplex.dimension() == 3
        c = self.circumcenter([vertices[i] for i in simplex.vertices])
        dist = min([self.distance(c, vertices[i]) for i in surface_vertices])
        if not simplex.inside:
            dist *= -1
        return dist

    def circumcenter(self, vertices):
        """Actually I compute the center."""
        x = sum(v[0] for v in vertices) / 4
        y = sum(v[1] for v in vertices) / 4
        z = sum(v[2] for v in vertices) / 4
        return [x, y, z]
    
    def distance(self, v1, v2):
        return sum((v2[i] - v1[i])**2 for i in range(3))
    
    def hdvf_with_quartets(self, fn='quartets.json'):
        """Build a HDVF having a given homology and cohomology basis."""
        with open(fn) as f:
            quartets = json.load(f)
        # Fix critical cells
        criticals = []
        for quartet in quartets:
            for side in ('inside', 'outside'):
                # assert len(set(quartet[side][0]) & set(quartet[side][1])) == 1
                critical = list(set(quartet[side][0]) & set(quartet[side][1]))[0]
                criticals.append(critical)
        # Make HDVF inside
        self._perfect_hdvf_quartets(quartets, criticals, True)
        self._perfect_hdvf_quartets(quartets, criticals, False)
        self.number_critical_cells()
        
    
    def _perfect_hdvf_quartets(self, quartets, criticals, inside):
        """Compute a perfect HDVF in K or in L-K given a fixed basis of homology and cohomology."""
        for quartet in quartets:
            pair = quartet['inside'] if inside else quartet['outside']
            # Make secondary
            for t in pair[0]:
                if t not in criticals and self.partition[t] == 'c': # @note I added the second condition but I don't undferstand why
                    # print(t, self.partition[t], self.get_d(t))
                    candidates = [s for s in self.get_d(t) if self.is_inside(s, inside) and s not in criticals]
                    if candidates:
                        s = candidates[0]
                        self.hdvf_a(s, t)
            # Make primary inside
            for s in pair[1]:
                if s not in criticals and self.partition[s] == 'c':
                    candidates = [t for t in self.get_d_star(s) if self.is_inside(t, inside) and t not in criticals]
                    if candidates:
                        t = candidates[0]
                        self.hdvf_a(s, t)
        # Fill the rest of the HDVF
        indices = list(range(self.nb_cells()))
        # random.shuffle(indices)
        print(f"\nPerfect HDVF inside={inside}...")
        n_cells = sum(1 for t in indices if self.is_inside(t, inside))
        n_paired_cells = 0
        for t in indices:
            if self.partition[t] == 'c' and self.is_inside(t, inside) and t not in criticals:
                candidates = [s for s in self.get_d(t) if self.is_inside(s, inside) and s not in criticals]
                if candidates:
                    s = random.choice(candidates)
                    self.hdvf_a(s, t)
                    n_paired_cells += 2
                    if n_paired_cells % 100 == 0:
                        print(f"Progress: {100*n_paired_cells/n_cells:.2f} % ({n_paired_cells} / {n_cells})", end="\r")
    

    def hdvf_with_bases(self, bases: dict):
        """Build a HDVF having a given homology and cohomology basis."""
        # Fix critical cells
        cells = {k: set(chain.from_iterable(v)) for k, v in bases.items()}
        criticals = list(cells['incycles'] & cells['incocycles']) + list(cells['outcycles'] & cells['outcocycles'])
        print("Number of critical cells to leave: ", len(criticals), criticals)
        # Make HDVF inside
        self._perfect_hdvf_bases(bases['incycles'], bases['incocycles'], criticals, True)
        self._perfect_hdvf_bases(bases['outcycles'], bases['outcocycles'], criticals, False)
        self.number_critical_cells()
        
    
    def _perfect_hdvf_bases(self, cycles, cocycles, criticals, inside):
        """Compute a perfect HDVF in K or in L-K given a fixed basis of homology and cohomology."""
        # Make secondary
        for cycle in cycles:
            for t in cycle:
                if t not in criticals and self.partition[t] == 'c': # @note I added the second condition but I don't understand why
                    candidates = [s for s in self.get_d(t) if self.is_inside(s, inside) and s not in criticals]
                    if candidates:
                        s = candidates[0]
                        self.hdvf_a(s, t)
        # Make primary
        for cocycle in cocycles:
            for s in cocycle:
                if s not in criticals and self.partition[s] == 'c': # @note I added the second condition but I don't understand why
                    candidates = [t for t in self.get_d_star(s) if self.is_inside(t, inside) and t not in criticals]
                    if candidates:
                        t = candidates[0]
                        self.hdvf_a(s, t)
        # Fill the rest of the HDVF
        indices = list(range(self.nb_cells()))
        # random.shuffle(indices)
        print(f"\nPerfect HDVF inside={inside}...")
        n_cells = sum(1 for t in indices if self.is_inside(t, inside))
        n_paired_cells = 0
        for t in indices:
            if self.partition[t] == 'c' and self.is_inside(t, inside) and t not in criticals:
                candidates = [s for s in self.get_d(t) if self.is_inside(s, inside) and s not in criticals]
                if candidates:
                    s = random.choice(candidates)
                    self.hdvf_a(s, t)
                    n_paired_cells += 2
                    if n_paired_cells % 100 == 0:
                        print(f"Progress: {100*n_paired_cells/n_cells:.2f} % ({n_paired_cells} / {n_cells})", end="\r")


    def make_pairs(self):
        """Make the stitching pairs."""
        self.number_critical_cells()
        pairs = {}
        print('Isomorphism:')
        for cell in self.cells:
            if self.partition[cell.index] == 'c':
                if not cell.inside:
                    print(f"d_{cell.dimension()}({cell.index}) = {self.get_d(cell.index)}")
                pairs[cell.index] = self.pair_inside(cell.index, cell.inside)
        indices = list(range(self.nb_cells()))
        # random.shuffle(indices)
        for t in indices:
            if self.partition[t] == 'c' and self.dim(t) == 2:
                candidates = [s for s in self.get_d(t)]
                if candidates:
                    # s = random.choice(candidates)
                    s = candidates[0]
                    self.pairs.append({
                        'inside': pairs[s],
                        'outside': pairs[t]
                    })
                    self.hdvf_a(s, t)
    
    def pair_inside(self, cell: int, inside: bool):
        """Compute the pair cycle-cocycle of a critical cell"""
        assert self.partition[cell] == 'c'
        cycle = [x for x in self.get_g(cell) if self.is_inside(x, inside)]
        cocycle = [x for x in self.get_f_star(cell) if self.is_inside(x, inside)]
        return (cycle, cocycle)
    
    def stats(self, complex: Complex):
        """Compute statistics of the homological quartets."""
        optim = ChainOptimization(complex)
        stats = {'incycle': [], 'incocycle': [], 'outcycle': [], 'outcocycle': []}
        for pair in self.pairs:
            stats['incycle'].append(optim.evaluate(pair['inside'][0], cycle=True))
            stats['incocycle'].append(optim.evaluate(pair['inside'][1], cycle=False))
            stats['outcycle'].append(optim.evaluate(pair['outside'][0], cycle=True))
            stats['outcocycle'].append(optim.evaluate(pair['outside'][1], cycle=False))
        print('Statistiques of current quartets:')
        print(f"- incycle:    {sum(stats['incycle']):.2f}")
        print(f"- incocycle:  {sum(stats['incocycle']):.2f}")
        print(f"- outcycle:   {sum(stats['outcycle']):.2f}")
        print(f"- outcocycle: {sum(stats['outcocycle']):.2f}")
    

    def improve(self, complex: Complex):
        optim = ChainOptimization(complex)
        n_steps = 10000
        for pair in self.pairs:
            pair['inside'] = (optim.simulated_annealing(pair['inside'][0], True, True, n_steps), 
                              optim.simulated_annealing(pair['inside'][1], False, True, n_steps))
            pair['outside'] = (optim.simulated_annealing(pair['outside'][0], True, False, n_steps), 
                               optim.simulated_annealing(pair['outside'][1], False, False, n_steps))
        with open('quartets.json', 'w') as outfile:
            outfile.write(json.dumps(self.pairs))
    

    def improve_quartets(self, complex: Complex):
        for pair in self.pairs:
            x = complex.boundary(pair['outside'][0]) # candidate incycle
            x = [c for c in x if complex.simplices[c].inside]
            y = complex.coboundary(pair['inside'][1]) # candidate outcocycle
            y = [c for c in y if not complex.simplices[c].inside]
            pair['inside'] = (x, pair['inside'][1]) 
            pair['outside'] = (pair['outside'][0], y)
        # with open('quartets.json', 'w') as outfile:
        #     outfile.write(json.dumps(self.pairs))

    def cofaces_cocycles(self, cocycle, inside: bool):
        """
        Get the list of the cofaces of a cocycle.
        Note that the coboundary of a cocycle is zero (by definition)
        """
        cofaces = []
        for simplex_i in cocycle:
            cofaces.extend(self.cofaces(simplex_i, inside))
        return list(set(cofaces))
    
    def non_standard_outside_cocycle(self, cocycle):
        segments = []
        for simplex_i in cocycle:
            cofaces = self.cofaces(simplex_i, False)
            if len(cofaces) == 2:
                v0 = set(self.cells[cofaces[0]].vertices)
                v1 = set(self.cells[cofaces[1]].vertices)
                vertices = list(v0 & v1) + list(v0.symmetric_difference(v1))
                segments.append(vertices)
        return segments
    
    
    def write(self, complex: Complex, fn='output.json'):
        """Write the output file for the viewer."""
        data = {
            'vertices': complex.vertices, 
            'triangles': complex.surface(),
            'quartets': []
        }
        for pair in self.pairs:
            data['quartets'].append({
                'in_cycle': [self.cells[c].vertices for c in pair['inside'][0]],
                'in_cocycle': [self.cells[c].vertices for c in pair['inside'][1]],
                'in_cocycle_nonstandard': [self.cells[c].vertices for c in self.cofaces_cocycles(pair['inside'][1], inside=True)],
                'out_cycle': [self.cells[c].vertices for c in pair['outside'][0]],
                'out_cocycle': [self.cells[c].vertices for c in pair['outside'][1]],
                'out_cocycle_nonstandard': self.non_standard_outside_cocycle(pair['outside'][1])
            })
        with open("output.json", 'w') as outfile:
            outfile.write(json.dumps(data))

    
def main(name="socket"):
    complex = Complex("data/" + name + ".1.ele", "data/" + name + ".1.node")
    complex.number_simplices()

    hq = HomologicalQuartet(complex)
    # hq.random_perfect()
    hq.distance_hdvf(complex.vertices, complex.surface_vertices())
    hq.make_pairs()
    hq.write(complex, 'output0.json')

    hq.improve(complex)
    hq.write(complex, 'output.json')

def testing(name="socket"):
    complex = Complex("data/" + name + ".1.ele", "data/" + name + ".1.node")
    complex.number_simplices()
    hq = HomologicalQuartet(complex)
    hq.hdvf_with_basis()
    hq.make_pairs()
    # hq.improve(complex)
    hq.stats(complex)
    hq.improve_quartets(complex)
    hq.stats(complex)
    hq.improve(complex)
    hq.stats(complex)
    hq.write(complex, 'output.json')

def testing2(name="socket"):
    """
    I obtain the cycles and cocycles, reduce their sizes, and try to make quartets with them
    """
    complex = Complex("data/" + name + ".1.ele", "data/" + name + ".1.node")
    complex.number_simplices()
    
    hq = HomologicalQuartet(complex)
    hq.random_perfect()
    # hq.distance_hdvf(complex.vertices, complex.surface_vertices())
    hq.make_pairs()
    hq.stats(complex)
    hq.write(complex, 'output0.json')

    optim = ChainOptimization(complex)
    bases = optim.improve(hq.pairs)
    hq.hdvf_with_bases(bases)
    hq.make_pairs()
    hq.stats(complex)
    hq.write(complex, 'output.json')

if __name__ == "__main__":
    # main()
    # testing()
    testing2()
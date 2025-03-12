import json
import random
import copy
import math
import networkx as nx
import numpy as np

from itertools import combinations
from collections import deque

class Simplex:
    """
    Class implementing a simplex.
    A simplex has:
    - A unique index in its complex
    - A list `vertices` of (indices of) its vertices
    - A boolean `inside` that is true if the simplex belong to K (see class Complex)
    - A list `boundary` of (the indices) its faces
    - A list `coboundary` of (the indices) its cofaces
    """
    def __init__(self, index: int, vertices_str: str, inside: bool):
        self.index = index
        self.vertices = [int(v) for v in vertices_str.split(':')]
        self.inside = inside
        self.boundary = []
        self.coboundary = []

    def dimension(self):
        """Return the dimension of the simplex."""
        return len(self.vertices) - 1
    
    def vertices_str(self):
        """Return a string encoding the vertices of the simplex."""
        return ':'.join(f"{v}" for v in self.vertices)

    def faces(self):
        """Return a list of the faces, described by the list of vertices."""
        if self.dimension() == 0:
            return []
        return list(combinations(self.vertices, self.dimension()))


def find_shortest_path(A: int, C: int, vertices: list, edges: list) -> list:
    """
    Construit un graphe à partir de 'edges' (liste de tuples (v1, v2))
    en utilisant les coordonnées dans 'vertices', puis retourne le plus court
    chemin (liste d'indices de sommets) entre A et C en pondérant chaque arête
    par la distance euclidienne.
    """
    G = nx.Graph()
    for (u, v) in edges:
        p_u = np.array(vertices[u])
        p_v = np.array(vertices[v])
        weight = np.linalg.norm(p_u - p_v)
        G.add_edge(u, v, weight=weight)
    try:
        path = nx.shortest_path(G, source=A, target=C, weight='weight')
    except nx.NetworkXNoPath:
        path = []
    return path


class Complex:
    """
    A pair of simplicial complexes K \subset L.
    The simplices in K are labelled as "inside", while the simplices in L-K are "outside".
    A complex has:
    - A list of simplices. Their index is their position in this list.
    - A list of vertices. Each vertex is a list with its three coordinates. 
      The index of a vertex is its position in this list.
    """
    def __init__(self, fn_ele, fn_node):
        self.simplices = []
        self.vertices = []
        self.make_complex(fn_ele)
        self.make_vertices(fn_node)

    def make_complex(self, fn):
        """Make the complex from the list of tetrahedra."""
        # Read the file
        with open(fn, "r") as file:
            n_lines = [int(num) for num in file.readline().strip().split()][0]
            print(f"Number of tetrahedra: {n_lines}")
            lines = [[int(num) for num in file.readline().strip().split()[1:]] for _ in range(n_lines)] 
        # Index the simplices and tell if they belong to K
        # We use a dictionnary `index` that maps each simplex (uniquely represented by
        # the list of the indices of its vertices, sorted and stringified) to a pair (index, inside).
        index = {}
        for line in lines:
            self.add_faces(index, line)
        print(f"Number of simplices: {len(index)}")
        # Make the simplices
        for k, v in index.items():
            simplex = Simplex(v['index'], k, v['inside'])
            self.simplices.append(simplex)
        # Compute the boundary and coboundary
        for simplex in self.simplices:
            for face in simplex.faces():
                face_str = ':'.join(f"{v}" for v in face)
                # print(simplex.vertices_str(), face_str)
                assert face_str in index
                subsimplex = self.simplices[index[face_str]['index']]
                simplex.boundary.append(subsimplex.index)
                subsimplex.coboundary.append(simplex.index)


    def add_faces(self, index, line):
        """
        Add a tetrahedron and its faces to the dictionnary `index`.
        The list `line` contains the indices of the vertices of the tetrahedron.
        Note that a tetrahedron is marked as "outside" if its 4th vertex is -1.
        """
        tet = sorted(line[0:4])
        inside = (line[4] != -1)
        for q in range(4):
            for simplex in combinations(tet, q+1):
                simplex_str = ':'.join(f"{v}" for v in simplex)
                if simplex_str not in index:
                    index[simplex_str] = {'index': len(index), 'inside': inside}
                elif inside:
                    index[simplex_str]['inside'] = True


    def make_vertices(self, fn):
        """Read the list of vertices from the *.node file."""
        with open(fn, "r") as file:
            n_lines = [int(num) for num in file.readline().strip().split()][0]
            print(f"Number of vertices: {n_lines}")
            self.vertices = [[float(num) for num in file.readline().strip().split()[1:]] for _ in range(n_lines)] 
    
    def write(self, fn='complex.json'):
        """Write the complex into a JSON file."""
        data = {
            'complex': [simplex.__dict__ for simplex in self.simplices],
            'vertices': self.vertices
        }
        with open(fn, 'w') as outfile:
            outfile.write(json.dumps(data, indent=2))
    
    def number_simplices(self):
        """Print the number of simplices of each dimension."""
        count = [0, 0, 0, 0]
        for simplex in self.simplices:
            count[simplex.dimension()] += 1
        print("Number of simplices: ", count)

    def surface(self):
        """List of 2-simplices in the boundary of the inside complex K, well oriented."""
        triangles = []
        for simplex in self.simplices:
            if simplex.dimension() == 3 and simplex.inside:
                oriented = self.orientation(simplex)
                for f in simplex.boundary:
                    if len([i for i in self.simplices[f].coboundary if self.simplices[i].inside]) == 1:
                        triangle = copy.copy(self.simplices[f].vertices)
                        for j in range(4):
                            if simplex.vertices[j] not in self.simplices[f].vertices:
                                if (j % 2 == 0) != oriented:
                                    triangle.reverse()
                        triangles.append(triangle)
        return triangles
    
    def surface_vertices(self):
        """Return the list of vertices in the surface of the complex."""
        vertices = []
        for triangle in self.surface():
            vertices += triangle
        return list(set(vertices))

    
    def orientation(self, tetrahedron: Simplex):
        """Retur true if the vertices of the tetrahedron are oriented positively."""
        p0 = self.vertices[tetrahedron.vertices[0]]
        p1 = self.vertices[tetrahedron.vertices[1]]
        p2 = self.vertices[tetrahedron.vertices[2]]
        p3 = self.vertices[tetrahedron.vertices[3]]
        m = [
            [p3[i] - p0[i] for i in range(3)],
            [p3[i] - p1[i] for i in range(3)],
            [p3[i] - p2[i] for i in range(3)]
        ]
        det = m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] 
        return det >= 0
    
    def write_surface(self):
        """Write the surface mesh in a JSON file."""
        data = {
            'vertices': self.vertices, 
            'triangles': self.surface()
        }
        with open("mesh.json", 'w') as outfile:
            outfile.write(json.dumps(data))
        return data
    
    def boundary(self, chain: list) -> list:
        """Compute the boundary of a chain of simplices."""
        x = set()
        for c in chain:
            for s in self.simplices[c].boundary:
                if s in x:
                    x.remove(s)
                else:
                    x.add(s)
        return list(x)
    
    def coboundary(self, chain: list) -> list:
        """Compute the coboundary of a chain of simplices."""
        x = set()
        for c in chain:
            for s in self.simplices[c].coboundary:
                if s in x:
                    x.remove(s)
                else:
                    x.add(s)
        return list(x)
    
    def isCycle(self, chain: list) -> bool:
        """
        Détermine si une chaîne de simplexes (supposés être des arêtes) forme un cycle.
        """
        boundary_counts = {}
        
        # Parcours de chaque arête dans la chaîne
        for c in chain:
            edge = self.simplices[c]
            # On vérifie que le simplexe est bien une arête (dimension 1)
            if edge.dimension() != 1:
                raise ValueError("La chaîne doit être composée d'arêtes (1-simplexes).")
            # Comptage des sommets de l'arête
            for v in edge.vertices:
                boundary_counts[v] = boundary_counts.get(v, 0) + 1
                
        # Vérifier que chaque sommet apparaît un nombre pair de fois
        for count in boundary_counts.values():
            if count % 2 != 0:
                return False
        return True

    def sort_chain(self, chain):
        """
        Trie les arêtes de la chaîne pour que chaque arête ait un sommet en commun avec la précédente.
        """
        if not chain:
            return []

        sorted_chain = [chain[0]]
        remaining = set(chain[1:])

        while remaining:
            last_edge = self.simplices[sorted_chain[-1]]
            last_vertices = set(last_edge.vertices)

            found = False
            for edge_index in list(remaining):
                edge = self.simplices[edge_index]
                if any(v in last_vertices for v in edge.vertices):
                    sorted_chain.append(edge_index)
                    remaining.remove(edge_index)
                    found = True
                    break

            if not found:
                raise ValueError("Impossible d'ordonner la chaîne avec des arêtes connectées.")

        return sorted_chain
    def display(self, chain=None):
        """
        Permet de visualiser le complexe avec la librairie pyVista
        """
        import numpy as np
        import pyvista as pv

        plotter = pv.Plotter()

        # conversion des sommets
        vertices_np = np.array(self.vertices)

        # recuperation des faces
        triangles = self.surface()

        # Construction de la liste de connectivités attendue par PyVista :
        # Pour chaque triangle, on doit avoir : [3, i, j, k]
        faces = []
        for tri in triangles:
            faces.append(3)
            faces.extend(tri)
        faces_np = np.array(faces)

        # Création du maillage avec une transparence pour pouvoir voir les cycles
        mesh = pv.PolyData(vertices_np, faces_np)
        plotter.add_mesh(mesh, color="lightblue", opacity=0.7, show_edges=True)

        # affichage des sommets
        plotter.add_points(vertices_np, color="black", point_size=5, render_points_as_spheres=True)

        # Si une chaîne d'arêtes est fournie en param , on affiche alors ce cycle
        if chain is not None:
            for edge_index in chain:
                if edge_index < len(self.simplices):
                    simplex = self.simplices[edge_index]
                    if simplex.dimension() == 1:
                        pts = np.array([self.vertices[i] for i in simplex.vertices])
                        line = pv.Line(pts[0], pts[1])
                        # Affichage en bleu avec une épaisseur plus importante pour bien distinguer
                        plotter.add_mesh(line, color="blue", line_width=4)

        # Affichage interactif
        plotter.show()

    def polyscopeView(self, chain):
        """
        Permet de visualiser le complexe en utilisant la librairie Polyscope.
        Les faces sont affichées en transparence pour pouvoir voir la chaîne.
        """
        import polyscope as ps
        import numpy as np

        ps.init()

        # Conversion des sommets et récupération des triangles (faces) de la surface
        vertices = np.array(self.vertices)
        triangles = np.array(self.surface())  # Chaque face doit être une liste de 3 indices

        # Enregistrement du maillage de surface avec transparence
        ps_mesh = ps.register_surface_mesh("Complexe", vertices, triangles)
        ps_mesh.set_transparency(0.5)  

        # Construction de la chaîne d'arêtes passée en paramètre (on prend uniquement les arêtes)
        chain_edges = []
        for edge_index in chain:
            edge = self.simplices[edge_index]
            if edge.dimension() == 1:
                chain_edges.append(edge.vertices)
        chain_edges = np.array(chain_edges)

        # Enregistrement et affichage de la chaîne d'arêtes en rouge
        ps_chain = ps.register_curve_network("cycle", vertices, chain_edges)
        ps_chain.set_color((1, 0, 0))  # Couleur rouge

        ps.show()


class ChainOptimization:
    """Algorithms for computing a shorter homologuous cycle."""
    def __init__(self, complex: Complex):
        self.complex = complex
    
    def greedy(self, chain: list, cycle: bool, inside: bool):
        """
        Shorten the (co)cycle by adding a (co)boundary.
        - `chain` is a list of (?)
        - `cycle` is True if the chain is a cycle. Otherwise, it is a cocycle.
        - `inside` is True if the chain is in K. Otherwise, it is in L-K.
        """
        print(f"Improving chain of size {len(chain)} -> ", end='')
        short_chain = chain
        while True:
            chain2 = self.local_greedy(short_chain, cycle, inside)
            if chain2 is not None:
                short_chain = chain2
                # print(len(short_chain), ' -> ', end='')
            else:
                break
        print(len(short_chain))
        return short_chain
    
    def local_greedy(self, chain: list, cycle: bool, inside: bool):
        """
        Return a shorter homologuous chain.
        Let us say, that `chain` is a cycle of dimension 1, `cycle` is True and `inside` is True.
        Then, we look for a triangle adjacent to some edge in the chain that has at least 2 faces in the chain.
        If we find it, we compute the symmetric difference, which is a shorter cycle.
        If we do not find it, return None.
        
        This function actually can take cycles or cocycles, inside of outside K, of any dimension.
        """
        cur_chain = set(chain)
        for c in cur_chain:
            candidates = self.complex.simplices[c].coboundary if cycle else self.complex.simplices[c].boundary
            candidates = [s for s in candidates if self.complex.simplices[s].inside == inside]
            for t in candidates:
                other = self.complex.simplices[t].boundary if cycle else self.complex.simplices[t].coboundary
                other = [s for s in other if self.complex.simplices[s].inside == inside]
                other = set(other)
                if len(cur_chain & other) > len(other)/2:
                    cur_chain = cur_chain.symmetric_difference(other)
                    return list(cur_chain)
        return None
    
    def simulated_annealing(self, chain: list, cycle: bool, inside: bool, n_steps: int):
        """
        Optimize the size of a chain with simulated annealing.
        - `chain` is a list of simplices
        - `cycle` is True if the chain is a cycle (otherwise, it is a cocycle).
        - `inside` is True if the chain is in K (otherwise, it is in L-K).
        - `n_steps` is the number of steps in the simmulated annealing algorithm.
        """
        cur_sol = {'chain': copy.copy(chain), 'value': self.evaluate(chain, cycle)}
        best_sol = cur_sol
        print(f"Improving chain of size ({cur_sol['value']:.2f}, {len(chain)}) -> ", end='')
        for i in range(n_steps):
            next_sol = self.local_sa(cur_sol['chain'], cycle, inside)
            delta_size = next_sol['value'] - cur_sol['value']
            temperature = 1 - i/n_steps
            if delta_size < 0 or random.random() < math.exp(-delta_size / temperature):
                cur_sol = next_sol
                if cur_sol['value'] < best_sol['value']:
                    best_sol = copy.copy(cur_sol)
                    print(f"{cur_sol['value']:.2f} ({len(cur_sol['chain'])}, {i})", end=" -> ", flush=True)
        print('')
        return best_sol['chain']
    
    def evaluate(self, chain: list, cycle: bool):
        """Measure a chain of dimension 1 or 2."""
        if not cycle:
            return len(chain)
        return sum(self.volume(self.complex.simplices[c]) for c in chain)
    

    def local_sa(self, chain: list, cycle: bool, inside: bool):
        """
        Local perturbation of a chain for the simulated annealing algorithm.
        We choose a random triangle adjacent to the chain and add its boundary to the chain.
        """
        cur_chain = set(chain)
        count = 0
        while True:
            count += 1
            assert count < 1000
            c = random.choice(chain)
            candidates = self.complex.simplices[c].coboundary if cycle else self.complex.simplices[c].boundary
            candidates = [s for s in candidates if self.complex.simplices[s].inside == inside]
            if candidates:
                t = random.choice(candidates)
                other = self.complex.simplices[t].boundary if cycle else self.complex.simplices[t].coboundary
                other = [s for s in other if self.complex.simplices[s].inside == inside]
                other = set(other)
                next_chain = list(cur_chain.symmetric_difference(other))
                next_value = self.evaluate(next_chain, cycle)
                return {'chain': next_chain, 'value': next_value}
            
    def volume(self, cell: Simplex):
        """Compute the volume of a cell (length for 1-simplices, area for 2-simplices)"""
        if cell.dimension() == 1:
            p1 = self.complex.vertices[cell.vertices[0]]
            p2 = self.complex.vertices[cell.vertices[1]]
            return math.sqrt(sum((p2[i] - p1[i])**2 for i in range(3)))    
        if cell.dimension() == 2:
            p1 = self.complex.vertices[cell.vertices[0]]
            p2 = self.complex.vertices[cell.vertices[1]]
            p3 = self.complex.vertices[cell.vertices[2]]
            v1 = [p1[i] - p2[i] for i in range(3)]
            v2 = [p1[i] - p3[i] for i in range(3)]
            x = v1[1] * v2[2] - v1[2] * v2[1]
            y = v1[2] * v2[0] - v1[0] * v2[2]
            z = v1[0] * v2[1] - v1[1] * v2[0]
            return math.sqrt(x**2 + y**2 + z**2) / 2
        assert False # We should never reach this line
        return None
    
    def improve(self, pairs):
        """
        Improve the homology and cohomology basis of K and L-K. The steps are:
        1. Initialize candidates
        2. Improve outcycle and incocycles. This sets 
        3. Add their respective boundaries and coboundaries to the candidates
        4. Improve incycles and outcocycles
        5. Extract shortest inside homology basis and outside cohomology basis
        """
        f_in = [pair['inside'][1] for pair in pairs]
        g_in = [pair['inside'][0] for pair in pairs]
        f_out = [pair['outside'][1] for pair in pairs]
        g_out = [pair['outside'][0] for pair in pairs]
        n_steps = 100000
        incocycles = [self.simulated_annealing(y, False, True, n_steps) for y in f_in]
        outcycles = [self.simulated_annealing(x, True, False, n_steps) for x in g_out]
        candidate_incycles = [self.simulated_annealing(x, True, True, n_steps) for x in g_in]
        candidate_outcocycles = [self.simulated_annealing(y, False, False, n_steps) for y in f_out]
        for cycle in outcycles:
            x = self.complex.boundary(cycle)
            x = [c for c in x if self.complex.simplices[c].inside]
            candidate_incycles.append(self.simulated_annealing(x, True, True, n_steps))
        for cocycle in incocycles:
            y = self.complex.coboundary(cocycle)
            y = [c for c in y if not self.complex.simplices[c].inside]
            candidate_outcocycles.append(self.simulated_annealing(y, False, False, n_steps))
        incycles = self.shortest_basis(candidate_incycles, f_in)
        outcocycles = self.shortest_basis(candidate_outcocycles, g_out)
        return {'incycles': incycles, 'incocycles': incocycles, 'outcycles': outcycles, 'outcocycles': outcocycles}

    def shortest_basis(self, chains, map):
        """
        Extract the shortest basis from a set of cycles or cocycles.
        TODO This is not finished
        """
        # Sort chains by volume
        chains.sort(key=lambda x: self.evaluate(x, True))
        # Apply the map to the chains
        matrix = []
        for chain in chains:
            col = set()
            for i, row in enumerate(map):
                if len(set(chain) & set(row)) % 2 == 1:
                    col.add(i)
            matrix.append(col)
        # Obtain the first independent columns
        indices = []
        for j in range(len(matrix)): # for each column
            if matrix[j]:
                indices.append(j)
                i = next(iter(matrix[j]))
                for k in range(j, len(matrix)):
                    if i in matrix[k]:
                        matrix[k] = matrix[k].symmetric_difference(matrix[j])
        return [chains[j] for j in indices]
                    
    def optimize_cycle(self, chain: list, cycle: bool = True, inside: bool = True) -> list:
        """
        Notre premièere approche : optimisation locale
        """
        current_chain = set(chain)
        improvement_found = True

        while improvement_found:
            improvement_found = False
            best_chain = current_chain
            best_length = len(current_chain)
            
            # Parcours de chaque simplex du cycle actuel
            for s in current_chain:
                if cycle:
                    # Pour un cycle, on cherche dans la cobordure de s
                    candidates = [t for t in self.complex.simplices[s].coboundary 
                                  if self.complex.simplices[t].inside == inside]
                else:
                    # Pour un cocycle, on regarde dans le bord de s
                    candidates = [t for t in self.complex.simplices[s].boundary 
                                  if self.complex.simplices[t].inside == inside]

                # Pour chaque candidat, on calcule la modification proposée
                for t in candidates:
                    if cycle:
                        # On prend l'ensemble des faces de t qui sont dans le sous-complexe concerné
                        candidate_faces = set(u for u in self.complex.simplices[t].boundary 
                                              if self.complex.simplices[u].inside == inside)
                    else:
                        # Pour un cocycle, on prend la cobordure de t
                        candidate_faces = set(u for u in self.complex.simplices[t].coboundary 
                                              if self.complex.simplices[u].inside == inside)
                    
                    # La nouvelle chaîne est la différence symétrique entre la chaîne courante et candidate_faces
                    new_chain = current_chain.symmetric_difference(candidate_faces)
                    
                    # Si la nouvelle chaîne est strictement plus courte, on la retient comme meilleure amélioration
                    if len(new_chain) < best_length:
                        best_length = len(new_chain)
                        best_chain = new_chain
                        improvement_found = True
            
            # On met à jour la chaîne courante si une amélioration a été trouvée
            if improvement_found:
                current_chain = best_chain
        
        return list(current_chain)

    def optimize_cycle_global(self, chain: list, cycle: bool = True, inside: bool = True) -> list:
        """
        Optimise un cycle (ou cocycle) donné par une liste d'indices de simplexes, 
        en essayant de trouver le minimum global. 
        Si aucune amélioration locale n'est trouvée pendant un certain nombre d'itérations,
        une perturbation aléatoire est appliquée pour sortir du minimum local.
        """
        current_chain = set(chain)
        best_chain = current_chain.copy()
        best_length = len(current_chain)
        
        # Paramètres pour la recherche globale
        max_iter = 100000      # nombre maximum d'itérations globales
        max_no_improve = 5  # nombre d'itérations sans amélioration avant perturbation
        no_improve_count = 0

        for iteration in range(max_iter):
            improvement_found = False
            
            # Recherche locale : parcourir chaque simplex dans la chaîne actuelle
            for s in list(current_chain):
                # Sélectionner les candidats selon le type de chaîne
                if cycle:
                    candidates = [t for t in self.complex.simplices[s].coboundary 
                                if self.complex.simplices[t].inside == inside]
                else:
                    candidates = [t for t in self.complex.simplices[s].boundary 
                                if self.complex.simplices[t].inside == inside]
                
                # Pour chaque candidat, calculer la modification proposée
                for t in candidates:
                    if cycle:
                        candidate_faces = set(u for u in self.complex.simplices[t].boundary 
                                                if self.complex.simplices[u].inside == inside)
                    else:
                        candidate_faces = set(u for u in self.complex.simplices[t].coboundary 
                                                if self.complex.simplices[u].inside == inside)
                    
                    # Nouvelle chaîne obtenue par différence symétrique
                    new_chain = current_chain.symmetric_difference(candidate_faces)
                    
                    # Si la nouvelle chaîne est strictement plus courte, c'est une amélioration
                    if len(new_chain) < len(current_chain):
                        current_chain = new_chain
                        improvement_found = True
                        break  # on sort de la boucle sur les candidats
                
                if improvement_found:
                    break  # on recommence la recherche locale avec la nouvelle chaîne

            if improvement_found:
                no_improve_count = 0  # réinitialiser le compteur d'itérations sans amélioration
                if len(current_chain) < best_length:
                    best_chain = current_chain.copy()
                    best_length = len(current_chain)
            else:
                no_improve_count += 1
            
            # Si aucune amélioration n'est trouvée depuis longtemps, on effectue une perturbation
            if no_improve_count >= max_no_improve:
                # Choix aléatoire d'un simplex dans la chaîne actuelle
                s = random.choice(list(current_chain))
                if cycle:
                    candidates = [t for t in self.complex.simplices[s].coboundary 
                                if self.complex.simplices[t].inside == inside]
                else:
                    candidates = [t for t in self.complex.simplices[s].boundary 
                                if self.complex.simplices[t].inside == inside]
                if candidates:
                    t = random.choice(candidates)
                    if cycle:
                        candidate_faces = set(u for u in self.complex.simplices[t].boundary 
                                                if self.complex.simplices[u].inside == inside)
                    else:
                        candidate_faces = set(u for u in self.complex.simplices[t].coboundary 
                                                if self.complex.simplices[u].inside == inside)
                    # La perturbation est effectuée en prenant la différence symétrique
                    current_chain = current_chain.symmetric_difference(candidate_faces)
                no_improve_count = 0  # réinitialiser le compteur
            
        # Retourne la meilleure solution trouvée
        return list(best_chain)
    def _find_edge(self, v1: int, v2: int) -> int:
        """
        Retourne l'indice d'une arête (1-simplexe) dans le complexe reliant les sommets v1 et v2.
        Les sommets de l'arête sont comparés dans l'ordre croissant.
        """
        edge_vertices = sorted([v1, v2])
        for simplex in self.complex.simplices:
            if simplex.dimension() == 1 and sorted(simplex.vertices) == edge_vertices:
                return simplex.index
        return None

    def narrow(self, chain: list) -> list:
        """
        Améliore (resserre) un cycle donné par une liste d'indices d'arêtes,
        en cherchant à remplacer des portions du cycle par des arêtes directes
        dont le milieu est plus proche du centroïde du cycle, ce qui devrait
        rapprocher le cycle du trou.

        Retourne une nouvelle liste d'indices d'arêtes représentant le cycle resserré.
        """
        import math

        # 1. Reconstruction du cycle de sommets à partir de la chaîne d'arêtes.
        # On ne considère que les arêtes (1-simplexes).
        edges = [tuple(self.complex.simplices[c].vertices) 
                for c in chain if self.complex.simplices[c].dimension() == 1]

        # Construire un dictionnaire d'adjacence à partir des arêtes.
        adj = {}
        for (u, v) in edges:
            adj.setdefault(u, []).append(v)
            adj.setdefault(v, []).append(u)
        
        # reconstruction du cycle
        start = edges[0][0]
        cycle_vertices = [start]
        current = start
        prev = None
        while True:
            next_vertex = None
            for n in adj[current]:
                if n != prev:
                    next_vertex = n
                    break
            if next_vertex is None or next_vertex == start:
                break
            cycle_vertices.append(next_vertex)
            prev, current = current, next_vertex
        # Fermer le cycle
        if cycle_vertices[0] != cycle_vertices[-1]:
            cycle_vertices.append(cycle_vertices[0])
        
        # 2. Calcul du centroïde (approximatif centre du trou) en ignorant le dernier sommet redondant.
        coords = [self.complex.vertices[v] for v in cycle_vertices[:-1]]
        n = len(coords)
        centroid = [sum(c[i] for c in coords) / n for i in range(3)]
        
        # Fonction utilitaire pour la distance euclidienne entre deux points en 3D.
        def dist(p, q):
            return math.sqrt(sum((p[i] - q[i]) ** 2 for i in range(3)))
        
        # Fonction qui calcule le point milieu entre deux points.
        def midpoint(p, q):
            return [(p[i] + q[i]) / 2 for i in range(3)]
        
        # 3. Itération pour "narrow" le cycle.
        improved = True
        new_cycle = cycle_vertices[:-1]  # cycle sans la répétition finale
        while improved:
            improved = False
            L = len(new_cycle)
            # Pour chaque paire de sommets non consécutifs (mais pas la connexion finale)
            for i in range(L):
                for j in range(i + 2, L):
                    if i == 0 and j == L - 1:
                        continue  # ne pas rompre la fermeture du cycle
                    # Vérifier s'il existe une arête directe entre new_cycle[i] et new_cycle[j].
                    if self._find_edge(new_cycle[i], new_cycle[j]) is not None:
                        # Calcul du "milieu moyen" de la portion actuelle du cycle allant de i à j.
                        segment_midpoints = []
                        for k in range(i, j):
                            p1 = self.complex.vertices[new_cycle[k]]
                            p2 = self.complex.vertices[new_cycle[k + 1]]
                            segment_midpoints.append(midpoint(p1, p2))
                        avg_mid = [sum(m[l] for m in segment_midpoints) / len(segment_midpoints) for l in range(3)]
                        # Milieu de l'arête directe proposée.
                        p_direct1 = self.complex.vertices[new_cycle[i]]
                        p_direct2 = self.complex.vertices[new_cycle[j]]
                        direct_mid = midpoint(p_direct1, p_direct2)
                        # Si le milieu direct est plus proche du centroïde que le milieu moyen actuel,
                        # alors on remplace la portion du cycle par un raccourci direct.
                        if dist(direct_mid, centroid) < dist(avg_mid, centroid):
                            new_cycle = new_cycle[:i+1] + new_cycle[j:]
                            improved = True
                            break
                if improved:
                    break
        # fermeture pour obtenir un cycle
        if new_cycle[0] != new_cycle[-1]:
            new_cycle.append(new_cycle[0])
        
        # conversion du cycle de sommets en chaîne d'arêtes.
        edge_chain = []
        for i in range(len(new_cycle) - 1):
            edge_index = self._find_edge(new_cycle[i], new_cycle[i+1])
            edge_chain.append(edge_index)
        
        return edge_chain

    def find_shortest_path(self, A: int, C: int) -> list:
        """
        Construit un graphe à partir de toutes les arêtes (1-simplexes) du complexe,
        puis utilise l'algorithme de Dijkstra pour trouver le plus court chemin entre
        le sommet A et le sommet C, en pondérant chaque arête par la distance euclidienne.
        
        Retourne une liste d'indices de sommets qui représentent le chemin le plus court.
        Si aucun chemin n'est trouvé, retourne une liste vide.
        """
        import networkx as nx
        import numpy as np

        G = nx.Graph()
        # Parcourir tous les simplexes du complexe pour récupérer les arêtes
        for simplex in self.complex.simplices:
            if simplex.dimension() == 1:
                u, v = simplex.vertices
                p_u = np.array(self.complex.vertices[u])
                p_v = np.array(self.complex.vertices[v])
                weight = np.linalg.norm(p_u - p_v)
                G.add_edge(u, v, weight=weight)
        
        try:
            # Utiliser Dijkstra pour trouver le plus court chemin de A à C
            path = nx.dijkstra_path(G, source=A, target=C, weight='weight')
        except nx.NetworkXNoPath:
            path = []
        return path
    def find_shortest_path(self, A: int, C: int) -> list:
        """
        Construit un graphe à partir des arêtes (1-simplexes) du complexe,
        applique Dijkstra pour trouver le chemin le plus court entre A et C,
        et retourne une liste d'arêtes connectées.
        
        Chaque arête est représentée par son indice dans self.simplices.
        Si aucun chemin n'est trouvé, retourne une liste vide.
        """
        import networkx as nx
        import numpy as np

        G = nx.Graph()
        edge_map = {}  # Dictionnaire pour retrouver l'index de l'arête dans self.simplices

        # Construire le graphe avec les arêtes pondérées par la distance euclidienne
        for i, simplex in enumerate(self.complex.simplices):
            if simplex.dimension() == 1:  # Simplexe 1D = arête
                u, v = simplex.vertices
                p_u = np.array(self.complex.vertices[u])
                p_v = np.array(self.complex.vertices[v])
                weight = np.linalg.norm(p_u - p_v)
                G.add_edge(u, v, weight=weight)
                edge_map[(u, v)] = i
                edge_map[(v, u)] = i  # Pour permettre un accès dans les deux directions

        try:
            # Trouver le plus court chemin sous forme de liste de sommets
            vertex_path = nx.dijkstra_path(G, source=A, target=C, weight='weight')

            # Convertir la liste de sommets en liste d'arêtes connectées
            edge_path = [edge_map[(vertex_path[i], vertex_path[i+1])] for i in range(len(vertex_path) - 1)]
        
        except nx.NetworkXNoPath:
            edge_path = []

        return edge_path

    def process_edges(self, chain: list) -> list:
        """
        Trouve un cycle en reliant les arêtes de 'chain' partageant un sommet.
        Retourne uniquement les nouvelles arêtes ajoutées par l'optimisation.

        Retourne : list[int] -> indices des nouvelles arêtes formant un cycle.
        """
        # Extraction des arêtes sous forme (v1, v2) 
        edges = [tuple(self.complex.simplices[c].vertices)
                for c in chain if self.complex.simplices[c].dimension() == 1]

        new_edges_set = set()

        # Parcours des arêtes consécutives
        for i in range(len(edges) - 1):
            edge1 = edges[i]
            edge2 = edges[i + 1]

            # Trouver le sommet commun
            common = set(edge1).intersection(edge2)
            if common:
                # A et C sont les sommets distincts
                B = list(common)[0]
                A = edge1[0] if edge1[0] != B else edge1[1]
                C = edge2[0] if edge2[0] != B else edge2[1]

                # Calculer le plus court chemin entre A et C
                path = self.find_shortest_path(A, C)

                # Ajouter les nouvelles arêtes qui ne sont pas déjà dans la chaîne initiale
                for p in path:
                    if p not in new_edges_set:
                        new_edges_set.add(p)

        return list(new_edges_set)

if __name__ == "__main__":
    complex = Complex("data/socket.1.ele", "data/socket.1.node")
    complex.number_simplices()
    # complex.write("simp_complex.json")
    # complex.write_surface()


    chain = [7363, 16196, 16133, 25157, 7756, 7181, 13838, 7182, 9104, 8782, 11538, 21780, 18069, 8665, 7897, 11739, 18588, 24865, 3171, 23590, 23463, 9450, 21869, 2354, 5874, 13556, 6969]
    
    op = complex.sort_chain(chain)
    print("hum",op,len(op))
    
    optim = ChainOptimization(complex)
    #v = optim.greedy(chain, True, True)
    #recuit_sm = optim.simulated_annealing(chain, True, True, 1000000)
    k = optim.optimize_cycle(chain, True, True)
    c = optim.optimize_cycle_global(chain,True,True)
    l = optim.narrow(op)
    m = optim.process_edges(op)
    complex.display(m)
    #print("recuit simulé",recuit_sm,len(recuit_sm))
    #print("premiere approche",k,len(k))

    # en partant deux aretes incidentes 
    #trouver la liste des triangles qui les reéunits
    #ajouter les bords des triangles
    
    #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
    
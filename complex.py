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
    
    def findEdge(self, a: int, b: int) -> int:
        """
        Retourne l'indice de l'arête (1-simplexe) dans le complexe reliant les sommets a et b.
        """
        edge_vertices = sorted([a, b])
        for simplex in self.simplices:
            if simplex.dimension() == 1 and sorted(simplex.vertices) == edge_vertices:
                return simplex.index
        return None

    def distance(self, edge_list: list) -> float:
        """
        Calcule la distance euclidienne totale pour une liste d'arêtes.
        Chaque arête est représentée par son indice dans self.complex.simplices.
        """
        import math
        total_distance = 0.0
        for edge_index in edge_list:
            edge = self.simplices[edge_index]
            if edge.dimension() == 1:
                a, b = edge.vertices # indices des points a et b
                p_a = self.vertices[a]
                p_b = self.vertices[b]
                # Calcul de la distance euclidienne entre les points a et b
                dist_ab = math.sqrt(sum((p_a[i] - p_b[i])**2 for i in range(3)))
                total_distance += dist_ab
        return total_distance


    def isCycle(self, chain: list) -> bool:
        """
        Détermine si une chaîne de simplexes (supposés être des arêtes) forme un cycle.
        """
        boundary_counts = {}
        
        # Parcours de chaque arête dans la chaîne
        for c in chain:
            edge = self.simplices[c]
            # Comptage des sommets de l'arête
            for v in edge.vertices:
                boundary_counts[v] = boundary_counts.get(v, 0) + 1
                
        # checking que chaque sommet apparaît un nombre pair de fois
        for count in boundary_counts.values():
            if count % 2 != 0:
                return False
        return True

    def gravityCenter(self, cycle: list) -> list:
        """
        Calcule le centre de gravité d'un cycle d'arêtes.
        
        Paramètre:
        cycle : list[int]
            Liste d'indices d'arêtes (1-simplexes) formant un cycle.
        
        Retourne:
        list[float] : Coordonnées (x, y, z) du centre de gravité.
        
        La méthode procède en deux étapes :
        1. Reconstruire la chaîne ordonnée des sommets via get_sorted_vertices.
        2. Calculer la moyenne des coordonnées de ces sommets.
            Si le cycle est fermé (premier sommet répété en fin), le sommet final est ignoré.
        """
        # Récupération de la liste ordonnée des sommets
        vertices_cycle = self.get_sorted_vertices(cycle)
        
        # Si le cycle est fermé (le premier sommet est répété en fin), on l'enlève
        if len(vertices_cycle) > 1 and vertices_cycle[0] == vertices_cycle[-1]:
            vertices_cycle = vertices_cycle[:-1]
        
        # Calcul du centre de gravité en moyennant les coordonnées
        center = [0.0, 0.0, 0.0]
        for v in vertices_cycle:
            coords = self.vertices[v]
            center[0] += coords[0]
            center[1] += coords[1]
            center[2] += coords[2]
        
        n = len(vertices_cycle)
        if n == 0:
            return center
        center = [c / n for c in center]
        return center


    def get_sorted_vertices(self, chain):
        """
        Trie les arêtes de la chaîne pour former une séquence ordonnée de sommets connectés.
        
        La chaîne en entrée est une liste d'indices d'arêtes (1-simplexes). On commence par une arête
        quelconque et on la suit en enchaînant les sommets connectés, garantissant une orientation correcte.
        
        Retourne une liste ordonnée de sommets.
        
        Exemple :
        Si les arêtes sont (a,b), (b,c), (c,d), la fonction retourne [a, b, c, d].
        """
        if not chain:
            return []
        
        # Construire un dictionnaire sommet -> arêtes incidentes
        vertex_to_edges = {}
        for edge_idx in chain:
            edge = self.simplices[edge_idx]
            a, b = edge.vertices
            vertex_to_edges.setdefault(a, []).append(edge_idx)
            vertex_to_edges.setdefault(b, []).append(edge_idx)

        # Trouver un sommet de départ (un sommet d'une extrémité si cycle ouvert, sinon n'importe lequel)
        start_vertex = None
        for v, edges in vertex_to_edges.items():
            if len(edges) == 1:  # Un sommet d'extrémité
                start_vertex = v
                break
        if start_vertex is None:  # Si c'est un cycle fermé, on prend un sommet arbitraire
            start_vertex = next(iter(vertex_to_edges))
        
        # Reconstruction de la chaîne ordonnée des sommets
        sorted_vertices = [start_vertex]
        used_edges = set()
        
        while len(sorted_vertices) < len(chain) + 1:  # Nombre de sommets = nombre d'arêtes + 1
            current_vertex = sorted_vertices[-1]
            next_edge = None
            
            for edge_idx in vertex_to_edges.get(current_vertex, []):
                if edge_idx not in used_edges:
                    next_edge = edge_idx
                    break
            
            if next_edge is None:
                break  # Si on ne trouve plus d'arête à suivre, on arrête

            used_edges.add(next_edge)
            edge = self.simplices[next_edge]
            a, b = edge.vertices
            next_vertex = b if a == current_vertex else a
            sorted_vertices.append(next_vertex)

        return sorted_vertices

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
                    
    def _find_edge(self, a: int, b: int) -> int:
        """
        Retourne l'indice de l'arête (1-simplexe) dans le complexe reliant les sommets a et b.
        Les sommets de l'arête sont comparés dans l'ordre croissant.
        Si aucune arête n'est trouvée, retourne None.
        """
        edge_vertices = sorted([a, b])
        for simplex in self.complex.simplices:
            if simplex.dimension() == 1 and sorted(simplex.vertices) == edge_vertices:
                return simplex.index
        return None

    def find_shortest_path(self, A: int, C: int) -> list:
        """
        Construit un graphe à partir des arêtes (1-simplexes) du complexe,
        utilise l'algorithme de Dijkstra pour trouver le plus court chemin entre
        le sommet A et le sommet C en la distance euclidienne comme poids.
        
        Retourne une liste d'indices d'arêtes qui représentent le chemin le plus court.
        Si aucun chemin n'est trouvé, retourne une liste vide.
        """
        import networkx as nx
        import numpy as np

        G = nx.Graph()
        edge_map = {}  # Permet d'associer chaque paire de sommets à l'indice de l'arête correspondante.
        
        # Construction du graphe à partir de toutes les arêtes du complexe.
        for i, simplex in enumerate(self.complex.simplices):
            if simplex.dimension() == 1 and simplex.inside:  # On ne considère que les 1-simplexes (arêtes)
                a, b = simplex.vertices
                p_a = np.array(self.complex.vertices[a])
                p_b = np.array(self.complex.vertices[b])
                weight = np.linalg.norm(p_a - p_b) # j utilise la distance ecludienne comme poids
                G.add_edge(a, b, weight=weight)
                edge_map[(a, b)] = i
                edge_map[(b, a)] = i  # Pour un accès dans les deux sens

        try:
            # Utiliser Dijkstra pour trouver le chemin le plus court sous forme de liste de sommets.
            vertex_path = nx.dijkstra_path(G, source=A, target=C, weight='weight')
        except nx.NetworkXNoPath:
            return []

        # Convertir la liste de sommets en liste d'arêtes grâce à l'edge_map.
        edge_path = []
        for j in range(len(vertex_path) - 1):
            u = vertex_path[j]
            v = vertex_path[j+1]
            edge_index = edge_map.get((u, v))
            if edge_index is None:
                # Cette situation ne devrait normalement pas se produire.
                continue
            edge_path.append(edge_index)
        return edge_path

    def minPathOptmisationStep(self, chain: list, jump: int = 3) -> list:
        """
        À partir d'une chaîne ordonnée de sommets (liste d'indices de sommets),
        calcule le chemin le plus court (liste d'indices d'arêtes) entre chaque sommet et celui
        situé 'jump' positions plus loin dans l'ordre.
        
        Si le cycle n'est pas fermé (le premier sommet n'est pas répété en fin), il est fermé.
        Lorsqu'on atteint la fin, le segment relie le dernier sommet au premier.
        
        Paramètres:
        chain : list[int]
            Liste ordonnée des indices de sommets formant le cycle.
        jump : int, optionnel (par défaut 3)
            Nombre de positions à sauter pour définir le segment à optimiser.
            Par exemple, jump=3 signifie qu'on relie le sommet 0 au sommet 3, 
            le sommet 3 au sommet 6, etc.
        
        Retourne:
        list[int] : Liste des indices d'arêtes (1-simplexes) obtenue en concaténant les chemins
        les plus courts entre les paires de sommets ainsi déterminées.
        """
        # S'assurer que le cycle est fermé
        if chain[0] != chain[-1]:
            chain.append(chain[0])
        
        n = len(chain) - 1  # nombre de sommets uniques (sans la répétition finale)
        new_edges = []
        i = 0
        while i < n:
            target = i + jump
            if target >= n:
                target = 0  # fermer le cycle en reliant au premier sommet
            u = chain[i]
            v = chain[target]
            edge_path = self.find_shortest_path(u, v)
            new_edges.extend(edge_path)
            i += jump
        return new_edges



    
    def minPathOptmisation(self,chain):
        """
        """
        sorted_vertices = self.complex.get_sorted_vertices(chain) # liste des sommets ordonnées du cycle
        oldDist = self.complex.distance(chain) # distance ecludienne totale des aretes du cycle actuel
        cycle = self.minPathOptmisationStep(sorted_vertices) # on reduit le cycle en une étape

        while True: # tant qu on arrive à reduire le cycle alors on continue
            sorted_vertices = self.complex.get_sorted_vertices(cycle)
            cycle = self.minPathOptmisationStep(sorted_vertices,3)
            dist = self.complex.distance(cycle)
            if dist >= oldDist : # si ça ne s est pas amélioré après une étape alors on arête
                break
            oldDist = dist
        return cycle
    
    def approachCycleToCenter(self, cycle: list) -> list:
        """
        Raffine un cycle en approchant ses arêtes vers le centre de gravité.
        
        Pour chaque paire consécutive de sommets dans le cycle (obtenu via get_sorted_vertices),
        la fonction cherche un triangle (2-simplexe) contenant ces deux sommets dont le troisième
        sommet (x) est plus proche du centre de gravité du cycle que l'un des deux sommets (a ou b).
        
        Si un tel triangle est trouvé, le cycle est mis à jour en prenant la différence symétrique
        entre le cycle actuel et le bord de ce triangle, ce qui revient à "ajouter" le triangle.
        
        Paramètre:
        cycle : list[int]
            Liste d'indices d'arêtes (1-simplexes) formant le cycle.
        
        Retourne:
        list[int] : Liste d'indices d'arêtes formant le cycle raffiné.
        """
        import numpy as np

        # Calculer le centre de gravité du cycle (la méthode utilise get_sorted_vertices en interne)
        center = self.complex.gravityCenter(cycle)
        
        # Récupérer la liste ordonnée des sommets du cycle
        sorted_vertices = self.complex.get_sorted_vertices(cycle)
        
        # Représenter le cycle courant sous forme d'ensemble d'arêtes (pour faciliter la modification)
        current_cycle = set(cycle)
        
        improved = False
        # Parcourir chaque paire consécutive de sommets dans le cycle
        for i in range(len(sorted_vertices) - 1):
            a = sorted_vertices[i]
            b = sorted_vertices[i + 1]
            # Examiner tous les triangles du complexe
            for simplex in self.complex.simplices:
                if simplex.dimension() == 2:
                    verts = simplex.vertices
                    # Vérifier que le triangle contient bien a et b
                    if a in verts and b in verts:
                        # Déterminer le troisième sommet (x) du triangle
                        third = [v for v in verts if v != a and v != b]
                        if not third:
                            continue
                        x = third[0]
                        # Calculer la distance du centre aux points a, b et x
                        coords_center = np.array(center)
                        dist_a = np.linalg.norm(np.array(self.complex.vertices[a]) - coords_center)
                        dist_b = np.linalg.norm(np.array(self.complex.vertices[b]) - coords_center)
                        dist_x = np.linalg.norm(np.array(self.complex.vertices[x]) - coords_center)
                        # Si x est plus proche du centre que a ou b, on considère ce triangle
                        if dist_x < dist_a or dist_x < dist_b:
                            # Récupérer les indices des arêtes du triangle
                            tri_edges = set()
                            for pair in [(verts[0], verts[1]), (verts[1], verts[2]), (verts[0], verts[2])]:
                                edge_idx = self.complex.findEdge(pair[0], pair[1])
                                if edge_idx is not None:
                                    tri_edges.add(edge_idx)
                            # Mettre à jour le cycle : la différence symétrique ajoute le triangle
                            current_cycle = current_cycle.symmetric_difference(tri_edges)
                            improved = True
                            # On passe au segment suivant dès qu'un triangle est ajouté pour celui-ci
                            break
        if improved:
            return list(current_cycle)
        else:
            return cycle

if __name__ == "__main__":
    complex = Complex("data/socket.1.ele", "data/socket.1.node")
    complex.number_simplices()
    # complex.write("simp_complex.json")
    # complex.write_surface()

    # chaine initiale encerclant un trou (donné par le prof)
    chain = [7363, 16196, 16133, 25157, 7756, 7181, 13838, 7182, 9104, 8782, 11538, 21780, 18069, 8665, 7897, 11739, 18588, 24865, 3171, 23590, 23463, 9450, 21869, 2354, 5874, 13556, 6969]
    
    optim = ChainOptimization(complex)

    print("longueur initiale ",complex.distance(chain))
    cycle = optim.minPathOptmisation(chain)
    cycle = optim.approachCycleToCenter(cycle)
    cycle = optim.minPathOptmisation(cycle)
    cycle = optim.approachCycleToCenter(cycle)
    cycle = optim.minPathOptmisation(cycle)
    print("longeur après optimisation",len(cycle),complex.distance(cycle))

    #v = optim.greedy(chain, True, True)
    #recuit_sm = optim.simulated_annealing(chain, True, True, 1000000)

    #k = optim.optimize_cycle(chain, True, True)
    #c = optim.optimize_cycle_global(chain,True,True)
    complex.display(cycle)
    #print("recuit simulé",recuit_sm,len(recuit_sm))
    #print("premiere approche",k,len(k))

    # en partant deux aretes incidentes 
    #trouver la liste des triangles qui les reéunits
    #ajouter les bords des triangles
    
    #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
    
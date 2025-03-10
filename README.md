# Homological Quartets

Scripts for computing the homology quartets. These tools take a surface mesh, build a tetrahedralization and compute the homology quartets. You can view them [here](https://pageperso.lis-lab.fr/aldo.gonzalez-lorenzo/papers/alexander-duality/quartets.html).

The basis steps are :

1. Compute a tetrahedralization of the convex hull of a 3D surface mesh with *TetGen*
2. Compute the homology quartets using random HDVF with `hdvf.py`
3. View them with `website/index.html`

See the next sections for more details.

## TetGen

Download it from [here](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) and install it using `cmake`.

There are some examples files in `data/`. Move them to the folder where you built tetgen and run `./tetgen -pcz socket.poly`. Other options with refinements are:

- `./tetgen -pcz -q socket.poly` for a refinement using angles
- `./tetgen -pcz -q -a100.0 socket.poly` for a refinement using the volume of the tetrahedra

TetGen produces the files `socket.1.node` and `socket.1.ele` (among others), that you must move to `data`.

## HDVF

Rename the input file in `hdvf.py` and run it. It computes the homology quartets and produces the file `output.json`.

## TODO

- [x] I do not really compute the stitching pairs
- [x] Use thickness-breadth balls to obtain better cycles and cocycles
- [x] Improve the inside 1-cycles by adding boundaries
- [x] Improve chain with simulated annealing
# cycle_optimisation

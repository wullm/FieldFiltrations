Reads a 3D binary data cube of floats (e.g. an output box from 21cmFAST),
calculates a periodic Delaunay triangulation of the corresponding grid
using CGAL, and then computes a field filtration using the data cube values
with GUDHI.

A field filtration is a superlevel (or sublevel) filtration of the triangulation
based on the grid values. See section 2.3 of Elbers & van de Weygaert (2018).

Elbers, Willem, and Rien van de Weygaert. "Persistent topology of the reionisation bubble network. I: Formalism & Phenomenology.", Monthly Notices of the Royal Astronomical Society, in print (2018), arXiv:1812.00462.

The data cube format is INDEX = z + DIM * (y + DIM * x).

Dependencies:
-- CGAL (https://www.cgal.org)
-- GUDHI (http://gudhi.gforge.inria.fr)
Field Filtrations
====
Reads a 3D binary data cube of floats (e.g. an output box from 21cmFAST),
calculates a periodic Delaunay triangulation of the corresponding grid
using CGAL, and then computes a field filtration using the data cube values
with GUDHI.

A field filtration is a superlevel (or sublevel) filtration of the triangulation
based on the grid values. See section 2.3 of Elbers & van de Weygaert (2018).

Elbers, Willem, and Rien van de Weygaert. "Persistent topology of the reionisation bubble network. I: Formalism & Phenomenology.", Monthly Notices of the Royal Astronomical Society, in print (2018), arXiv:1812.00462.

The data cube format is INDEX = z + DIM * (y + DIM * x).

Radial Filtrations
====
We also provide Radial Filtrations, an alternative method based on a
course-graining approach. Starting from a 3D data cube, we first reduce
the input to a binary field (1=active, 0=inactive) using a user-specified
threshold. The method works by radially growing (alpha>0) or shrinking
(alpha<0) the active regions.

We first grow the active regions by turning on cells that are adjacent to
active cells. We repeat this procedure a number of times until the entire
box is active. Cells that were activated after i steps are assigned a filtration
parameter alpha=i/N, where N is the width of the box. Next, we return to the
original binary field and shrink the regions in a similar manner. The remaining
cells are now assigned a negative filtration parameter alpha=−i/N, where i is
the step in which the cell was deactivated. The field thus obtained can be fed
into the Field Filtration method described above.



Dependencies:
-- CGAL (https://www.cgal.org)
-- GUDHI (http://gudhi.gforge.inria.fr)
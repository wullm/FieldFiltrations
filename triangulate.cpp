/**
 * @trian.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 *
 * Import 3D binary data cube and calculate a field filtration
 * of a Delaunay triangulation of a regular periodic grid.
 * See Elbers & van de Weygaert (2018).
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Computes the persistent homology of a superlevel filtration of a 
 * three-dimensional grid, using specified grid values.
 */
 
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_3_io.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

#include <field_io.h>

using Kernel		= CGAL::Exact_predicates_inexact_constructions_kernel;

using Point_3 		= Kernel::Point_3;

//We assign each vertex a pair of values (index, filtration value) as addintional info
using Info			= std::pair<int, double>;
using Point_info_pair = std::pair<Point_3, Info>;

//Equip the CGAL 3D Periodic Delaunay Triangulation to handle points_with_info
using Gt			= CGAL::Periodic_3_Delaunay_triangulation_traits_3<Kernel>;
using VbDS			= CGAL::Periodic_3_triangulation_ds_vertex_base_3<>;
using Vb			= CGAL::Triangulation_vertex_base_3<Gt, VbDS>;
using CbDS			= CGAL::Periodic_3_triangulation_ds_cell_base_3<>;
using Cb			= CGAL::Triangulation_cell_base_3<Gt, CbDS> ;
using VbInfo		= CGAL::Triangulation_vertex_base_with_info_3<Info, Gt, Vb>;
using TDS			= CGAL::Triangulation_data_structure_3<VbInfo, Cb>;
using P3DT3			= CGAL::Periodic_3_Delaunay_triangulation_3<Gt, TDS>;

//Short form CGAL handles
using Point			= P3DT3::Point;
using Iso_cuboid	= P3DT3::Iso_cuboid;
using Vertex_handle = P3DT3::Vertex_handle;
using Cell_handle 	= P3DT3::Cell_handle;
using Locate_type	= P3DT3::Locate_type;

//Custom sort traits for CGAL spatial_sort, allowing for points_with_info
using Search_traits_3 = CGAL::Spatial_sort_traits_adapter_3<Kernel,CGAL::First_of_pair_property_map<Point_info_pair>> ;

//The GUDHI types
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

//The 21cmFAST Box/Array format is INDEX = z + HII_D * (y + HII_D*x).
int box_idx(int width, int x, int y, int z) {
	return z + width * (y + width * x);	
}

int main(int argc, char * const argv[]) {
	//Help option
	if(argc == 2 && strcmp(argv[1], "--help")==0) {
		std::cout << "Usage: ./triangulate FILE [PRIME] [OPTIONS]" << std::endl;
		std::cout << "PRIME is the field coefficient (defaults is 2)." << std::endl;
		std::cout << "" << std::endl;
		std::cout << "  -i      changes to top-down mode (default is bottom-up)" << std::endl;
		std::cout << "  -o      write results to separate files per dimension" << std::endl;
		std::cout << "  --help  displays this message" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "Imports 3D binary periodic data cube and calculates a" << std::endl;
		std::cout << "field filtration of a Delaunay triangulation of the grid." << std::endl;
		std::cout << "See Elbers & van de Weygaert (2018) for details." << std::endl;
		return 0;
    }
	
	//FILE specified?
	if (argc == 1) {
		std::cout << "Usage: ./triangulate FILE [PRIME] [OPTIONS]" << std::endl;
		std::cout << "PRIME is the field coefficient (defaults is 2)." << std::endl;
		std::cout << "Try --help for more information." << std::endl;
		return 0;
	}
	
	//Check for options
	bool top_down = false;
	bool output_files = false;	
	for (int argi = 2; argi < argc; argi++) {
		if (strcmp(argv[argi], "-i")==0) {
			top_down = true;
		}
		if (strcmp(argv[argi], "-o")==0) {
			output_files = true;
		}
	}
	
	//The field characteristic (default is 2)
	int coeff_field_characteristic = 2;	
	//Determine whether a field coefficient was specified
	if (argc >= 3) {
		//Is it an integer?
		bool is_int = true;
		for (int stri = 0; stri < strlen(argv[2]); stri++) {
			if (!isdigit(argv[2][stri])) {
				is_int = false;
			}
		}
		
		if (is_int) {
			coeff_field_characteristic = std::stoi(argv[2]);
		}
	}
	
	//File name entered as first argument
	std::string fname(argv[1]);	
	std::vector<float> floats = read_floats(fname);
	
	int number_of_floats = floats.size();
	int width = cbrt(number_of_floats); //cube root
	
	std::cout << "#width: " << width << std::endl;
	std::cout << "p d birth death " << std::endl;
	
	//The fundamental domain
	Iso_cuboid domain(0,0,0,1,1,1);	
	double spacing = 1./width;
	
	//Determine the maximum field value (only needed for top-down mode)
	double max_val;
	if (top_down) {
		for (int x=0; x<width; x++) {
			for (int y=0; y<width; y++) {		
				for (int z=0; z<width; z++) {	
					if ((x==0 && y==0 && z==0) || floats[box_idx(width,x,y,z)] > max_val) {
						max_val = floats[box_idx(width,x,y,z)];
					}
				}
			}
		}
	}
	
		
	int i = 0;
	std::vector<Point_info_pair> points;
	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {		
			for (int z=0; z<width; z++) {								
				double fil_val;
				if (top_down) {
					fil_val = max_val - floats[box_idx(width,x,y,z)];					
				} else {
					fil_val = floats[box_idx(width,x,y,z)];
				}
				
				Info info = {i, fil_val};
				
				points.push_back({Point(x*spacing, y*spacing, z*spacing), info});
				
				i++;
			}
		}
	}

	//Sorting the points for fast insertion & faster runtime of CGAL algorithms
	std::random_shuffle (points.begin(), points.end()); //algorithm assumes shuffled set
	CGAL::spatial_sort(points.begin(), points.end(), Search_traits_3 ());

	//The periodic triangulation
	P3DT3 T(domain);
   
	//Custom fast point insertion (modified from the usual, allowing for points_with_info)
	Cell_handle hint;
	for (std::vector<Point_info_pair>::const_iterator p = points.begin();p != points.end(); ++p) {
		Locate_type lt;
		Cell_handle c;
		int li, lj;
		c = T.locate (p->first, lt, li, lj, hint);
		Vertex_handle v = T.insert (p->first, lt, c, li, lj);
		if ( v==P3DT3::Vertex_handle() ) {
			hint=c;
		} else {
			v->info() = p->second;
			hint=v->cell();
		}
	}
	
	if (T.is_triangulation_in_1_sheet()) {
		T.convert_to_1_sheeted_covering();
	} else {
		std::cerr << "ERROR: we were not able to construct a triangulation within a single periodic domain." << std::endl;
		exit(-1);
	}
	
	//Export the triangulation (e.g. for plotting)
	//std::ofstream to_off("output_regular.off"); // as a .off file
	//CGAL::write_triangulation_to_off(to_off, T);
		
  
	//Now iterate through all the simplices and insert them at the right value
	Simplex_tree simplexTree;
		
	//Vertices	
	P3DT3::Vertex_iterator vit;
	for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit) {		
		Info info  = vit->info();		
		std::vector<int> simplexVector = {info.first};				
		simplexTree.insert_simplex(simplexVector, Filtration_value(info.second));
	}
	
	//Edges
	P3DT3::Finite_edges_iterator eit;
	for(eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit) {		
		Info info1 = eit->first->vertex( eit->second )->info();
		Info info2 = eit->first->vertex( eit->third )->info();
		
		//Find the maximum filtration value of the vertices
		double fil_val = fmax(info1.second, info2.second);
					
		//Construct the simplex
		std::vector<int> simplexVector = {info1.first, info2.first};
		
		//Ignore trivial edges
		if (info1.first != info2.first) {
			simplexTree.insert_simplex(simplexVector, Filtration_value(fil_val));	
		}		
	}
	
	//Facets
	P3DT3::Finite_facets_iterator fit;
	for(fit = T.finite_facets_begin(); fit != T.finite_facets_end(); ++fit) {
		//Facets are specified by a cell (four vertices) and an opposite vertex
		int opposite_vertex = fit->second;
		
		Info vs[3];
		vs[0] = fit->first->vertex((opposite_vertex + 1) % 4)->info();
		vs[1] = fit->first->vertex((opposite_vertex + 2) % 4)->info();
		vs[2] = fit->first->vertex((opposite_vertex + 3) % 4)->info();		
		
		//Find the maximum filtration value of the vertices
		double fil_val = fmax(fmax(vs[0].second, vs[1].second), vs[2].second);
		
		//Construct the simplex
		std::vector<int> simplexVector = {vs[0].first, vs[1].first, vs[2].first};
		
		//Ignore trivial facets (i.e. only insert if no vertices coincide)
		if (vs[0].first != vs[1].first && vs[0].first != vs[2].first && vs[1].first != vs[2].first) {					
			simplexTree.insert_simplex(simplexVector, Filtration_value(fil_val));	
		}
	}
	
	//Cells
	P3DT3::Finite_cells_iterator cit;
	for(cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit) {
		//Get all four vertices
		Info vs[4];
		vs[0] = cit->vertex(0)->info();
		vs[1] = cit->vertex(1)->info();
		vs[2] = cit->vertex(2)->info();
		vs[3] = cit->vertex(3)->info();
		
		//Find the maximum filtration value of the vertices
		double fil_val = fmax(fmax(fmax(vs[0].second, vs[1].second), vs[2].second), vs[3].second);
		
		//Ignore trivial cells (i.e. only insert if no vertices coincide)
		if (vs[0].first != vs[1].first && vs[0].first != vs[2].first && vs[0].first != vs[3].first && vs[1].first != vs[2].first && vs[1].first != vs[3].first && vs[2].first != vs[3].first) {		
		
			//Construct the simplex
			std::vector<int> simplexVector = {vs[0].first, vs[1].first, vs[2].first, vs[3].first};
		
			simplexTree.insert_simplex(simplexVector, Filtration_value(fil_val));	
		}
	}
		
	//Ignore features with persistence less than
	double min_persistence = 0;
	
	//Compute the persistence diagram of the complex
	Persistent_cohomology pcoh(simplexTree);
	//Initialises the coefficient field for homology
	pcoh.init_coefficients(coeff_field_characteristic);

	pcoh.compute_persistent_cohomology(min_persistence);

	//Output the diagrams for d=0,1,2
	for (int dim = 0; dim < 3; dim++) {
		
		//The resulting (birth,death) intervals (only features of this dimension)
		std::vector<std::pair<double, double>> persistence_intervals = pcoh.intervals_in_dimension(dim);
			
		std::ofstream out_file;
		if (output_files) {
			out_file.open("persistence_"+std::to_string(dim)+".csv");
			out_file << "birth death\n";
		}
		
		for (auto pair : persistence_intervals) {
			double birth = pair.first;
			double death = pair.second;
			
			//Ignore trivial features (zero persistence)
			if (birth != death) {
				
				if (top_down) {
					birth = max_val - birth;				
					if (isinf(death)) {
						death = -INFINITY;
					} else {
						death = max_val - death;
					}				
				}	
				
				if (output_files) {
					out_file << birth << " " << death << "\n";
				}
				std::cout << coeff_field_characteristic << " " << dim << " " << birth << " " << death << std::endl;
			}
		}
		
		if (output_files) {		
			out_file.close();
		}	
	}
	
	//Output the diagram in filediag
	//pcoh.output_diagram();
		
	return 0;
}

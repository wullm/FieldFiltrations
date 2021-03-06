/**
 * @radial.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Computes the persistent homology of a radial filtration of a 
 * binary field specified on a periodic 3D grid. See the README
 * for more details.
 */
 
#include <iostream>
#include <sstream>
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
		std::cout << "Usage: ./radial FILE THRESHOLD [PRIME] [OPTIONS]" << std::endl;
		std::cout << "THRESHOLD is the value above (below in -i mode) which a cell is activated." << std::endl;
		std::cout << "PRIME is the field coefficient (defaults is 2)." << std::endl;
		std::cout << "" << std::endl;
		std::cout << "  -i      changes to inverted mode (active if <= threshold)" << std::endl;
		std::cout << "  -o      write results to separate files per dimension" << std::endl;
		std::cout << "  --help  displays this message" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "Imports 3D binary periodic data cube and calculates a" << std::endl;
		std::cout << "spatial filtration of a Delaunay triangulation of the grid." << std::endl;
		std::cout << "See Elbers & van de Weygaert (2018) for details." << std::endl;
		return 0;
    }
	
	//FILE specified?
	bool exit_with_message = false;
	if (argc == 1) { //no
		exit_with_message = true;
	}
	
	//Check for options
	bool inverted = false;
	bool output_files = false;	
	for (int argi = 2; argi < argc; argi++) {
		if (strcmp(argv[argi], "-i")==0) {
			inverted = true;
		}
		if (strcmp(argv[argi], "-o")==0) {
			output_files = true;
		}
	}
	
	//Determine the activation threshold
	float threshold = 1.0;
	if (argc >= 3) {
		//Is it a number?
		std::istringstream is(argv[2]);		
		is >> threshold;
		
		if (!is.eof() || is.fail()) {
			exit_with_message = true;
		}
	} else {
		exit_with_message = true;
	}
	
	//The field characteristic (default is 2)
	int coeff_field_characteristic = 2;	
	//Determine whether a field coefficient was specified
	if (argc >= 4) {
		//Is it a positive integer?
		bool is_int = true;
		for (int stri = 0; stri < strlen(argv[3]); stri++) {
			if (!isdigit(argv[3][stri])) {
				is_int = false;
			}
		}
		
		if (is_int) {
			coeff_field_characteristic = std::stoi(argv[3]);
		}
	}
		
	
	if (exit_with_message) {
		std::cout << "Usage: ./radial FILE THRESHOLD [PRIME] [OPTIONS]" << std::endl;
		std::cout << "THRESHOLD is the value above (below in -i mode) which a cell is activated." << std::endl;
		std::cout << "PRIME is the field coefficient (defaults is 2)." << std::endl;
		std::cout << "Try --help for more information." << std::endl;
		return 0;
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
	
	
	
	
	//First determine the active cells (points that exceed the threshold)
	std::vector<bool> active_cells(width*width*width);
	bool any_active_cells = false;
	bool any_inactive_cells = false;
	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {		
			for (int z=0; z<width; z++) {
								
				if ((!inverted && floats[box_idx(width,x,y,z)] >= threshold)
					|| (inverted && floats[box_idx(width,x,y,z)] <= threshold)) {
						
					active_cells[box_idx(width,x,y,z)] = true;
					any_active_cells = true;
				} else {
					active_cells[box_idx(width,x,y,z)] = false;
					any_inactive_cells = true;
				}
			}
		}
	}
	
	//At least one cell active?
	if (!any_active_cells) {
		std::cerr << "Not a single cell is active at this threshold." << std::endl;
		return 0;		
	}
	//At least one cell inactive?
	if (!any_inactive_cells) {
		std::cerr << "Not a single cell is inactive at this threshold." << std::endl;
		return 0;		
	}
	
	//The filtration parameter is alpha, which is the radial growing/shrinking scale
	
	//We will determine the filtration value (alpha) of all cells
	std::vector<float> fil_vals(width*width*width);
	
	//First consider the points with alpha = 0 (the active cells)
	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {		
			for (int z=0; z<width; z++) {
								
				if (active_cells[box_idx(width,x,y,z)]) {
					fil_vals[box_idx(width,x,y,z)] = 0.0;
				}				
			}
		}
	}
	
	//Next radially shrink the active regions to determine the cells with negative alpha
	bool done = false;
	float alpha = 0.0;
	while (!done) {		
		int active_cells_left = 0;
		alpha -= spacing;		
				
		//Find active cells that are adjacent to inactive cells and deactivate them
		std::vector<bool> next_active_cells(width*width*width);
		for (int x=0; x<width; x++) {
			for (int y=0; y<width; y++) {		
				for (int z=0; z<width; z++) {
		
					if (!active_cells[box_idx(width,x,y,z)]) {
						//Inactive cells remain inactive
						next_active_cells[box_idx(width,x,y,z)] = false;
					} else {
						//Active cells become inactive if adjacent to inactive
						bool inactive_adjacent = false;
		
						for (int sx = -1; sx <= 1; sx++) {
							for (int sy = -1; sy <= 1; sy++) {
								for (int sz = -1; sz <= 1; sz++) {
									if (x + sx >= 0 && y + sy >= 0 && z + sz >= 0
										&& x + sx < width && y + sy < width && z + sz < width) {
											
										if (!active_cells[box_idx(width,x+sx,y+sy,z+sz)]) {
											inactive_adjacent = true;
										}
									}
								}
							}
						}
						
						
						if (inactive_adjacent) {							
							next_active_cells[box_idx(width,x,y,z)] = false;
							
							//Deactivated during this step => assign filtration value
							fil_vals[box_idx(width,x,y,z)] = alpha;							
						} else {
							next_active_cells[box_idx(width,x,y,z)] = true;
							active_cells_left++;
						}				
						
					}
				}
			}
		}
		
		active_cells = next_active_cells;
		
		if (active_cells_left==0) {
			done = true;			
		}
	}
	
	//Done with points with negative alpha. Proceed with points with positive alpha.
	
	//Again determine the active cells (points that exceed the threshold)
	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {		
			for (int z=0; z<width; z++) {
								
				if ((!inverted && floats[box_idx(width,x,y,z)] >= threshold)
					|| (inverted && floats[box_idx(width,x,y,z)] <= threshold)) {
						
					active_cells[box_idx(width,x,y,z)] = true;
				} else {
					active_cells[box_idx(width,x,y,z)] = false;
				}
			}
		}
	}
	
	//Next radially grow the active regions to determine the cells with positive alpha
	done = false;
	alpha = 0.0;
	while (!done) {		
		int inactive_cells_left = 0;
		alpha += spacing;		
				
		//Find inactive cells that are adjacent to active cells and activate them
		std::vector<bool> next_active_cells(width*width*width);
		for (int x=0; x<width; x++) {
			for (int y=0; y<width; y++) {		
				for (int z=0; z<width; z++) {
		
					if (active_cells[box_idx(width,x,y,z)]) {
						//Active cells remain active
						next_active_cells[box_idx(width,x,y,z)] = true;
					} else {
						//Inactive cells become active if adjacent to active
						bool active_adjacent = false;
		
						for (int sx = -1; sx <= 1; sx++) {
							for (int sy = -1; sy <= 1; sy++) {
								for (int sz = -1; sz <= 1; sz++) {
									if (x + sx >= 0 && y + sy >= 0 && z + sz >= 0
										&& x + sx < width && y + sy < width && z + sz < width) {
											
										if (active_cells[box_idx(width,x+sx,y+sy,z+sz)]) {
											active_adjacent = true;
										}
									}	
								}
							}
						}
						
						
						if (active_adjacent) {							
							next_active_cells[box_idx(width,x,y,z)] = true;
							
							//Activated during this step => assign filtration value
							fil_vals[box_idx(width,x,y,z)] = alpha;							
						} else {
							next_active_cells[box_idx(width,x,y,z)] = false;
							inactive_cells_left++;
						}				
						
					}
				}
			}
		}
		
		active_cells = next_active_cells;
		
		if (inactive_cells_left==0) {
			done = true;			
		}
	}
	
	//Done with creating the radial field
	//write_floats("fil_vals.box", fil_vals);
		
	//Now that all cells have a filtration value, build the Delaunay triangulation + filtration
	int i = 0;
	std::vector<Point_info_pair> points;
	
	//First, insert the points
	for (int x=0; x<width; x++) {
		for (int y=0; y<width; y++) {		
			for (int z=0; z<width; z++) {
				double fil_val = fil_vals[box_idx(width,x,y,z)];
				
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

/**
 * @field_io.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
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
 * Loads the grid values from file.
 */
 
#include <iostream>
#include <fstream>
#include <vector>

//Read a vector of floats from a binary file
std::vector<float> read_floats(std::string fname) {			
	//Open the file (in, binary, open at end)
	std::fstream file (fname, std::fstream::in | std::fstream::binary | std::fstream::ate);
					
	if (file.is_open()) {
		int size = file.tellg(); //current pos at end
		int float_size = sizeof(float); //bytes
		int number_of_floats = size / float_size;
		int width = cbrt(number_of_floats); //cube root
		
		if (width < 1 || width*width*width != number_of_floats) {
			std::cerr << "File is not a data cube of floats." << std::endl;
			exit(0);			
		}
				
		std::vector<float> the_floats(number_of_floats);

		//Go to the beginning
		file.seekg(0, std::ios::beg);
		float f; int i = 0;
		while(file.read(reinterpret_cast<char*>(&f), float_size)) {
			the_floats[i] = f;
			i++;
		}	
			
		file.close();	
		
		return the_floats;
	} else {
		std::cerr << "File not found." << std::endl;
		exit(0);		
	}
		
}
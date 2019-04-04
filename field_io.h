/**
 * @field_io.h
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
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

void write_floats(std::string fname, std::vector<float> floats) {
	//Open the file (out, binary)
	std::fstream file (fname, std::fstream::out | std::fstream::binary);
	
	if (file.is_open()) {
		int float_size = sizeof(float); //bytes
		
		for (auto f : floats) {
			file.write(reinterpret_cast<const char*>(&f), float_size);
		}
		
		file.close();
	} else {
		std::cerr << "Could not write file." << std::endl;
		exit(0);		
	}
}
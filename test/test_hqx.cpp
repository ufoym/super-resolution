#include "../image.h"
#include "../segment.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>

int main(int argc, char* argv[])
{
	std::string fn_input = argv[1];
	std::string fn_output = argv[2];
	image<rgb>* im = load(fn_input.c_str());
	if (im) {
		time_t start = clock();
		Segment seg(im);
		std::cout << std::left << std::setw(50) << fn_output
				  << std::setw(7)
				  << float(clock()-start)/CLOCKS_PER_SEC 
				  << std::setw(3) << "sec" << std::endl;
		image<rgb> *seg_img = seg.vis();
		save(fn_output.c_str(), seg_img);
		delete seg_img;
	}
	else {
		std::cout << "file not support" << std::endl;
	}	
	return 0;
}
#include "image.h"
#include "segment.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{
	if(argc != 2) {
		std::cout 	<< "drag and drop an image file on me" 
					<< std::endl;
	}
	else {
		std::string fn_input = argv[1];
		std::string fn_seg = fn_input + ".seg.ppm";
		image<rgb>* im = load(fn_input.c_str());
		if (im) {
			Segment seg(im);
			image<rgb> *seg_img = seg.vis();
			save(fn_seg.c_str(), seg_img);
			delete seg_img;
		}
		else {
			std::cout << "file not support" << std::endl;
		}		
	}
	system("pause");
	return 0;
}
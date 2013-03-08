#include "image.h"
#include "segment.h"
#include "curvefit.h"
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


			std::vector<float>				nodes;
			std::vector<std::vector<int>>	indices;
			std::vector<bool>				junction_map;

			seg.trace(nodes, indices, junction_map);
			CurveFitter fitter(nodes, indices, junction_map);
			std::string fn_fit = fn_input + ".fit.svg";
			fitter.saveToSVG(fn_fit);


			delete seg_img;
		}
		else {
			std::cout << "file not support" << std::endl;
		}		
	}
	system("pause");
	return 0;
}
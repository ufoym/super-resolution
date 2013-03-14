#include <iostream>
#include <fstream>
#include <string>

#include "image.h"
#include "segment.h"
#include "spline.h"
#include "render.h"


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
			BSplineFitter fitter(nodes, indices, junction_map);
			std::string fn_fit = fn_input + ".fit.svg";
			fitter.saveToSVG(fn_fit);

			
			std::vector<BSpline> splines;
			fitter.getSplines(splines);
			fitter.clear();
			Renderer renderer(nodes, indices, junction_map, 
				splines, im, seg_img, seg.get_color_map());
			delete seg_img;
			delete im;
		}
		else {
			std::cout << "file not support" << std::endl;
		}		
	}
	system("pause");
	return 0;
}
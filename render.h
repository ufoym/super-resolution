#pragma once

#include <map>
#include <vector>
#include <string>
#include "image.h"
#include "spline.h"


class Renderer 
{
private:
	#define UP_SCALE		3
	#define PATCH_RADIUS	UP_SCALE

public:
	Renderer(std::vector<float>&			nodes,
			 std::vector<std::vector<int>>&	indices,
			 std::vector<bool>&				junction_map,
			 std::vector<BSpline>&			splines,
			 image<rgb>*					im,
			 image<rgb>*					seg_img,
			 std::map<int, rgb>&			color_map);

private:
	rgb _linearComb(rgb& fg, rgb& bg, double alpha);

	int _rgbDist(rgb& c1, rgb& c2);

	double _pointDist(double x1, double y1, double x2, double y2);

	int _findNeighbor(float x, float y, std::vector<float>& points);

	bool _pointInPolygon(float x, float y, std::vector<float>& polygon);
};
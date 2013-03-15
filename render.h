#pragma once

#include <map>
#include <vector>
#include <string>
#include "image.h"
#include "spline.h"


class Renderer 
{
private:
	#define UP_SCALE			3
	#define PATCH_RADIUS		UP_SCALE
	#define HIRES_PATCH_RADIUS	UP_SCALE
	#define LORES_PATCH_RADIUS	1

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
	
	double _calcLikelihood( int xc, int yc, BSpline& spline, image<rgb>* im, rgb& fg_color, rgb& bg_color, std::vector<float>& shape_polygon );

	rgb _calcBgColor(int xc, int yc, image<rgb>* seg_img, rgb fg_color);

	void _mergeCurvePolygon( std::vector<float>& samples, std::vector<float>& shape_polygon, std::vector<float>& new_polygon );

	void _getPatch( int xc, int yc, image<rgb> *im, rgb* patch, int radius );

	double _clacPatchDiff( rgb *p, rgb *q, int radius );

	void _renderPatch( int xc, int yc, int radius, std::vector<float>& new_polygon, rgb& fg_color, rgb& bg_color, rgb *patch );

	void _bicubicInterp( rgb *hi_patch, int hi_side, rgb *lo_patch, int lo_side );
};
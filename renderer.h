#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "image.h"
#include "bsplinefit.h"


class Renderer 
{
private:
	#define UP_SCALE	5

public:
	Renderer(std::vector<float>&			nodes,
			 std::vector<std::vector<int>>&	indices,
			 std::vector<bool>&				junction_map,
			 std::vector<BSpline>&			splines,
			 image<rgb>*					im,
			 image<rgb>*					seg_img,
			 std::map<int, rgb>&			color_map)
	{
		// allocate memory for the result.
		int width  = UP_SCALE * im->width();
		int height = UP_SCALE * im->height();
		image<rgb> *result = new image<rgb>(width, height);

		// scale up all the coordinates.
		for (int i = 0; i < nodes.size(); i++)
			nodes[i] *= 3;
		for (int i = 0; i < splines.size(); i++) {
			for (int j = 0; j < splines[i].c.size(); j++)
				splines[i].c[j] *= UP_SCALE;
			for (int j = 0; j < splines[i].data.size(); j++)
				splines[i].data[j] *= UP_SCALE;
			for (int j = 0; j < splines[i].samples.size(); j++)
				splines[i].samples[j] *= UP_SCALE;
		}

		// render patches on shape boundaries.
		for (int id_spline = 0; id_spline < splines.size(); id_spline++) {

			// step 1: determine the central node of this spline.
			BSpline& spline = splines[id_spline];
			int num_samples = spline.samples.size() / 2;
			int id_central = num_samples / 2;
			int xc = (int)(0.5 + spline.samples[2 * id_central]);		
			int yc = (int)(0.5 + spline.samples[2 * id_central + 1]); 
			xc = min(xc, width - 1);
			yc = min(yc, height - 1);

			// step 2: determine background color of the current 5x5 patch using the seg_img.
			rgb fg_color = color_map[spline.id_shape];
			rgb bg_color = seg_img->access[yc / UP_SCALE][xc / UP_SCALE];
			int yy = yc / UP_SCALE;
			int xx = xc / UP_SCALE;
			for (int dx = -1; dx <= 0; dx++) {
				for (int dy = -1; dy <= 0; dy++) {
					if (xx + dx >= 0 && yy + dy >= 0) {
						if (_rgbDist(seg_img->access[yy + dy][xx + dx], fg_color) > 1)
							bg_color = seg_img->access[yy + dy][xx + dx];
					}
				}
			}

			//for (int x = max(0, xc - 2); x < min(xc + 2, width); x++) {
			//	for (int y = max(0, yc - 2); y < min(yc + 2, height); y++) {
			//		if (_rgbDist(seg_img->access[y / UP_SCALE][x / UP_SCALE], fg_color) 
			//			> _rgbDist(bg_color, fg_color)) {
			//			bg_color = seg_img->access[y / UP_SCALE][x / UP_SCALE];
			//		}
			//	}
			//}

			// step 3: merge the spline with the polygon of current shape.
			int id_shape = spline.id_shape;
			int shape_size = indices[id_shape].size();
			std::vector<float> polygon(2 * shape_size);
			for (int i = 0; i < shape_size; i++) {
				int idx = indices[id_shape][i];
				polygon[2 * i] = nodes[2 * idx];
				polygon[2 * i + 1] = nodes[2 * idx + 1];
			}

			std::vector<float>& samples = spline.samples;
			std::pair<float, float> lo_endpt(samples[0], samples[1]);
			std::pair<float, float> hi_endpt(samples[samples.size() - 2], samples[samples.size() - 1]);
			std::vector<float> new_polygon = samples;

			int idx_lo_ngbr = _findNeighbor(lo_endpt.first, lo_endpt.second, polygon);
			int idx_hi_ngbr = _findNeighbor(hi_endpt.first, hi_endpt.second, polygon);

			for (int i = idx_hi_ngbr; i % shape_size != idx_lo_ngbr; i++) {
				int idx = i % shape_size;
				new_polygon.push_back(polygon[2 * idx]);
				new_polygon.push_back(polygon[2 * idx + 1]);
			}
			new_polygon.push_back(polygon[2 * idx_lo_ngbr]);
			new_polygon.push_back(polygon[2 * idx_lo_ngbr + 1]);

			// step 4: assign color to each pixel in the 5x5 patch by pointInPolygon test.
			for (int x = max(0, xc - 2); x < min(xc + 2, width); x++) {
				for (int y = max(0, yc - 2); y < min(yc + 2, height); y++) {
					if (_pointInPolygon(x + 0.5, y + 0.5, new_polygon))
						result->access[y][x] = fg_color;
					else
						result->access[y][x] = bg_color;
				}
			}
		}

		// fill color of the interior of each shape.

		save("render.ppm", result);
		delete result;
	}


	inline int _rgbDist(rgb& c1, rgb& c2) 
	{
		return ((int)c1.r - c2.r) * ((int)c1.r - c2.r)
			 + ((int)c1.g - c2.g) * ((int)c1.g - c2.g)
			 + ((int)c1.b - c2.b) * ((int)c1.b - c2.b);
	}


	inline double _pointDist(double x1, double y1, double x2, double y2)
	{
		return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
	}


	int _findNeighbor(float x, float y, std::vector<float>& points)
	{
		int n = points.size() / 2;
		int idx = 0;
		for (int i = 1; i < n; i++) {
			if (_pointDist(x, y, points[2 * i], points[2 * i + 1]) <
				_pointDist(x, y, points[2 * idx], points[2 * idx + 1]))
				idx = i;
		}
		return idx;
	}


	bool _pointInPolygon(float x, float y, std::vector<float>& polygon)
	{
		int poly_sides = polygon.size() / 2;
		int i, j = poly_sides - 1;
		bool odd_nodes = false;

		for (i = 0; i < poly_sides; i++) {
			if ((polygon[2 * i + 1] < y && polygon[2 * j + 1] >= y
			 ||  polygon[2 * j + 1] < y && polygon[2 * i + 1] >= y)
			 && (polygon[2 * i] <= x    || polygon[2 * j] <= x)) {
				if (polygon[2 * i] + (y - polygon[2 * i + 1]) / (polygon[2 * j + 1] - polygon[2 * i + 1]) 
					* (polygon[2 * j] - polygon[2 * i]) < x) {
					 odd_nodes = !odd_nodes;
				}
			}
			j = i;
		}
		return odd_nodes;
	}
};
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "svg.h"
#include "render.h"

Renderer::Renderer( std::vector<float>&				nodes, 
					std::vector<std::vector<int>>&	indices, 
					std::vector<bool>&				junction_map, 
					std::vector<BSpline>&			splines, 
					image<rgb>*						im, 
					image<rgb>*						seg_img, 
					std::map<int, rgb>&				color_map )
{
	bool debug_mode = false;
	int DEBUG_ID_SPLINE = 200;
	SVGWriter writer("debug.svg");

	// allocate memory for the result.
	int width  = UP_SCALE * im->width();
	int height = UP_SCALE * im->height();
	image<rgb> *result = new image<rgb>(width, height);
	//image<char> *visited = new image<char>(width, height);

	typedef struct{
		int cnt, r, g, b;
		void operator += (rgb& c) {
			cnt++;
			r += c.r;
			g += c.g;
			b += c.b;
		}
		rgb simplify() {
			rgb c = {0, 0, 0};
			if (cnt != 0) {
				c.r = r / cnt;
				c.g = g / cnt;
				c.b = b / cnt;
			}
			return c;
		}
	} xrgb;
	image<xrgb> *im_tmp = new image<xrgb>(width, height);

	// scale up all the coordinates.
	for (int i = 0; i < nodes.size(); i++)
		nodes[i] *= UP_SCALE;
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
		int xc = (int)(0.5 + spline.data[2 * id_central]);		
		int yc = (int)(0.5 + spline.data[2 * id_central + 1]); 
		xc = std::min(xc, width - 1);
		yc = std::min(yc, height - 1);

		// step 2: determine background color of the current 5x5 patch using the seg_img.
		rgb fg_color = color_map[spline.id_shape];
		int xx = xc / UP_SCALE;	 xx = std::min(xx, width / UP_SCALE - 1);
		int yy = yc / UP_SCALE;  yy = std::min(yy, height / UP_SCALE - 1);

		rgb bg_color = seg_img->access[yy][xx];
		for (int dx = -1; dx <= 0; dx++) {
			for (int dy = -1; dy <= 0; dy++) {
				if (xx + dx >= 0 && yy + dy >= 0) {
					if (debug_mode && id_spline == DEBUG_ID_SPLINE) {
						rgb color = seg_img->access[yy + dy][xx + dx];
						std::cout << "nb_color: " << (int)color.r << "," << (int)color.g << "," << (int)color.b << "\n";
					}
					if (_rgbDist(seg_img->access[yy + dy][xx + dx], fg_color) > 1)
						bg_color = seg_img->access[yy + dy][xx + dx];
				}
			}
		}


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

		if (debug_mode && id_spline == DEBUG_ID_SPLINE) {
			writer.writePolygon(new_polygon, "#00FF00");
			writer.writeDots(new_polygon, "#0000FF");
		}

		// step 4: assign color to each pixel in the 5x5 patch by pointInPolygon test.	
		for (int x = std::max(0, xc - PATCH_RADIUS); x < std::min(xc + PATCH_RADIUS, width - 1); x++) {
			for (int y = std::max(0, yc - PATCH_RADIUS); y < std::min(yc + PATCH_RADIUS, height - 1); y++) {
				// subdivide a pixel into 5x5 grid with 25 sensor points.
				// perform pointInPolygon test to these 25 points and calculate
				// the pixel color as linear combination of fore- & back-ground colors.
				int num_inside = 0;
				for (float dx = 0.1; dx < 1; dx += 0.2) {
					for (float dy = 0.1; dy < 1; dy += 0.2) {
						if (debug_mode && id_spline == DEBUG_ID_SPLINE
							&& x == xc + 1 && y == yc + 1) {
								bool in = _pointInPolygon(x + dx, y + dy, new_polygon);
								std::string color;
								if (in) {color = "#FFA700";}
								else    {color = "#000000";}
								writer.writeDot(x + dx, y + dy, color);
						}
						if (_pointInPolygon(x + dx, y + dy, new_polygon))
							num_inside++;
					}
				}
				// a naive way to deal with patch overlap: just average the pixel value
				// if current patch has already been visited, assuming that the current
				// patch will only be overlapped with the previous reconstructed patch.
				//if (visited->access[y][x]) {
				//	result->access[y][x] = _linearComb(
				//		_linearComb(fg_color, bg_color, num_inside / 25.0),
				//		result->access[y][x], 0.5);
				//} 
				//else {
				//	result->access[y][x] = _linearComb(fg_color, bg_color, num_inside / 25.0);
				//	visited->access[y][x] = 1;
				//}
				// a more sophisticated way to deal mix overlapped regions:
				// just accumulate the color value for each pixel, and normalize at the last step.
				im_tmp->access[y][x] += _linearComb(fg_color, bg_color, num_inside / 25.0);
				// dotted the current patch in debug mode.
				if (debug_mode && id_spline == DEBUG_ID_SPLINE) {
					writer.writeDot(x + 0.5, y + 0.5, "#FF0000");
				}
			}
		}
	}

	// average pixel values in im_tmp to result.
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			result->access[y][x] = im_tmp->access[y][x].simplify();
		}
	}

	// fill color of the interior of each shape.
	// TODO: to be implemented.

	writer.close();
	save("render.ppm", result);
	//delete visited;
	delete im_tmp;
	delete result;
}

rgb Renderer::_linearComb( rgb& fg, rgb& bg, double alpha )
{
	rgb c;
	c.r = alpha * fg.r + (1 - alpha) * bg.r;
	c.g = alpha * fg.g + (1 - alpha) * bg.g;
	c.b = alpha * fg.b + (1 - alpha) * bg.b;
	return c;
}

bool Renderer::_pointInPolygon( float x, float y, std::vector<float>& polygon )
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

int Renderer::_findNeighbor( float x, float y, std::vector<float>& points )
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

double Renderer::_pointDist( double x1, double y1, double x2, double y2 )
{
	return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

int Renderer::_rgbDist( rgb& c1, rgb& c2 )
{
	return ((int)c1.r - c2.r) * ((int)c1.r - c2.r)
		+ ((int)c1.g - c2.g) * ((int)c1.g - c2.g)
		+ ((int)c1.b - c2.b) * ((int)c1.b - c2.b);
}



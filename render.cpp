#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "svg.h"
#include "render.h"

struct xrgb{
	int cnt, r, g, b;
	xrgb() {cnt = r = g = b = 0;}
	void operator += (rgb& c) {
		cnt++;
		r += c.r;
		g += c.g;
		b += c.b;
	}
	rgb simplify() {
		rgb c(0, 0, 0);
		if (cnt != 0) {
			c.r = r / cnt;
			c.g = g / cnt;
			c.b = b / cnt;
		}
		return c;
	}
	void clear() {cnt = r = g = b = 0;}
};

bool debug_mode = false;

Renderer::Renderer( std::vector<float>&				nodes, 
					std::vector<std::vector<int>>&	indices, 
					std::vector<bool>&				junction_map, 
					std::vector<BSpline>&			splines, 
					image<rgb>*						im, 
					image<rgb>*						seg_img, 
					std::map<int, rgb>&				color_map )
{
	
	
	int DEBUG_ID_SPLINE = 1500;
	SVGWriter writer("debug.svg");

	// allocate memory for the result.
	int width  = UP_SCALE * im->width();
	int height = UP_SCALE * im->height();
	image<rgb> *result = new image<rgb>(width, height);
	//image<char> *visited = new image<char>(width, height);

	
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


	std::vector<int> idx_winners;
	std::vector<float> shape_polygon;
	int last_shape = -1;
	for (int id_spline = 0; id_spline < splines.size(); id_spline += NUM_CURVE_TYPES) {
		// update shape_polygon only if the shape changes.
		// this saves a lot of computation by not recomputing the shape_polygon 
		// for each node in the zig-zag boundary.
		int cur_shape = splines[id_spline].id_shape;
		if (cur_shape != last_shape) {
			shape_polygon.clear();
			for (int i = 0; i < indices[cur_shape].size(); i++) {
				int idx = indices[cur_shape][i];
				shape_polygon.push_back(nodes[2 * idx]);
				shape_polygon.push_back(nodes[2 * idx + 1]);
			}
		}
		last_shape = cur_shape;
		if (id_spline == DEBUG_ID_SPLINE) {debug_mode = true;}

		// determine foreground and background color for the current node.
		int n = splines[id_spline].data.size() / 2;
		int id_central = n / 2;
		int xc = splines[id_spline].data[2 * id_central];
		int yc = splines[id_spline].data[2 * id_central + 1];
		xc = std::min(xc, width - 1);
		yc = std::min(yc, height - 1);
		rgb fg_color = color_map[cur_shape];
		rgb bg_color = _calcBgColor(xc / UP_SCALE, yc / UP_SCALE, seg_img, fg_color);

		

		double max_likelihood = -1;
		int idx_best = -1;
		
		for (int i = 0; i < NUM_CURVE_TYPES; i++) {
			double lh = _calcLikelihood(xc, yc, splines[id_spline + i], im, fg_color, bg_color, shape_polygon);
			if (debug_mode) std::cout << std::endl;
			// std::cout << 1e+7 - lh << std::endl;
			if (lh > max_likelihood) {
				idx_best = i + id_spline;
				max_likelihood = lh;
			}
		}
		if (debug_mode) {
			std::cout << "fg_color: "
				<< (int)fg_color.r << ","
				<< (int)fg_color.g << ","
				<< (int)fg_color.b << std::endl;
			std::cout << "bg_color: "
				<< (int)bg_color.r << ","
				<< (int)bg_color.g << ","
				<< (int)bg_color.b << std::endl;
		}
		if (id_spline == DEBUG_ID_SPLINE) {
			debug_mode = false;
		}
		//std::cout << idx_best - id_spline << ",    " << idx_best << std::endl;
		idx_winners.push_back(idx_best);
	}

















	debug_mode = true;
	// render patches on shape boundaries.
	//for (int id_spline = 0; id_spline < splines.size(); id_spline++) {
	for (int idx_set = 0; idx_set < idx_winners.size(); idx_set++) {
		int id_spline = idx_winners[idx_set];


		if (id_spline == 1500) id_spline = 1504;
#undef DEBUG_ID_SPLINE
#define DEBUG_ID_SPLINE 1504


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
			writer.writePolyline(new_polygon, "#00FF00");
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

double Renderer::_calcLikelihood( int xc, int yc, BSpline& spline, image<rgb>* im, rgb& fg_color, rgb& bg_color, std::vector<float>& shape_polygon )
{
	if (debug_mode){
		int a = 3;
	}
	// merge curve and polygon.
	std::vector<float> new_polygon;
	_mergeCurvePolygon(spline.samples, shape_polygon, new_polygon);

	

	rgb hi_patch[4 * HIRES_PATCH_RADIUS * HIRES_PATCH_RADIUS];
	rgb lo_patch[4 * LORES_PATCH_RADIUS * LORES_PATCH_RADIUS];
	rgb ds_patch[4 * LORES_PATCH_RADIUS * LORES_PATCH_RADIUS];

	if (debug_mode) {
		for (int y = 0; y < 6; y++) {
			for (int x = 0; x < 6; x++) {
				int idx = y * 6 + x;
				//if (((int)hi_patch[idx].r) != 255
				//	|| ((int)hi_patch[idx].g) != 255
				//	|| ((int)hi_patch[idx].b) != 255)
				//	std::cout << "here" << std::endl;
				/*std::cout << "(" << (int)hi_patch[idx].r 
					      << "," << (int)hi_patch[idx].g
						  << "," << (int)hi_patch[idx].b
						  << ")\t";*/
			}
			//std::cout << "\n";
		}
	}

	
	_renderPatch(xc, yc, HIRES_PATCH_RADIUS, new_polygon, fg_color, bg_color, hi_patch);
	xrgb tmp[4 * HIRES_PATCH_RADIUS * HIRES_PATCH_RADIUS];
	for (int i = 0; i < 4 * HIRES_PATCH_RADIUS * HIRES_PATCH_RADIUS; i++)
		tmp[i] += hi_patch[i];

	_bicubicInterp(hi_patch, 2 * HIRES_PATCH_RADIUS, 
				   ds_patch, 2 * LORES_PATCH_RADIUS);
	_getPatch(xc / 3, yc / 3, im, lo_patch, LORES_PATCH_RADIUS);

	double err = _clacPatchDiff(lo_patch, ds_patch, LORES_PATCH_RADIUS);
	return 1e+7 - err;
}

void Renderer::_mergeCurvePolygon( std::vector<float>& samples, std::vector<float>& shape_polygon, std::vector<float>& new_polygon )
{
	// step 3: merge the spline with the polygon of current shape.
	int shape_size = shape_polygon.size() / 2;
	std::pair<float, float> lo_endpt(samples[0], samples[1]);
	std::pair<float, float> hi_endpt(samples[samples.size() - 2], samples[samples.size() - 1]);
	new_polygon = samples;

	int idx_lo_ngbr = _findNeighbor(lo_endpt.first, lo_endpt.second, shape_polygon);
	int idx_hi_ngbr = _findNeighbor(hi_endpt.first, hi_endpt.second, shape_polygon);

	for (int i = idx_hi_ngbr; i % shape_size != idx_lo_ngbr; i++) {
		int idx = i % shape_size;
		new_polygon.push_back(shape_polygon[2 * idx]);
		new_polygon.push_back(shape_polygon[2 * idx + 1]);
	}
	new_polygon.push_back(shape_polygon[2 * idx_lo_ngbr]);
	new_polygon.push_back(shape_polygon[2 * idx_lo_ngbr + 1]);
}

void Renderer::_getPatch( int xc, int yc, image<rgb> *im, rgb* patch, int radius )
{
	for (int dx = -radius; dx < radius; dx++) {
		for (int dy = -radius; dy < radius; dy++) {
			if (xc + dx >= 0 && xc + dx < im->width()
				&& yc + dy >= 0 && yc + dy < im->height())
				patch[(dy + radius) * (2 * radius) + (dx + radius)]
					= im->access[yc + dy][xc + dx];
			else
				patch[(dy + radius) * (2 * radius) + (dx + radius)]
					= rgb(0, 0, 0);
		}
	}
}

rgb Renderer::_calcBgColor( int xc, int yc, image<rgb>* seg_img, rgb fg_color )
{
	rgb bg_color = seg_img->access[yc][xc];
	for (int dx = -1; dx <= 0; dx++) {
		for (int dy = -1; dy <= 0; dy++) {
			if (xc + dx >= 0 && xc + dx < seg_img->width() 
				&& yc + dy >= 0 && yc + dy < seg_img->height()) {
				if (seg_img->access[yc + dy][xc + dx] != fg_color)
					bg_color = seg_img->access[yc + dy][xc + dx];
			}
		}
	}
	return bg_color;
}

double Renderer::_clacPatchDiff( rgb *p, rgb *q, int radius )
{
	int dist = 0;
	for (int i = 0; i < 4 * radius * radius; i++) {
		dist += _rgbDist(p[i], q[i]);
	}
	return dist;
}

void Renderer::_renderPatch( int xc, int yc, int radius, std::vector<float>& new_polygon, rgb& fg_color, rgb& bg_color, rgb *patch )
{
	for (int x = xc - radius; x < xc + radius; x++) {
		for (int y = yc - radius; y < yc + radius; y++) {
			// subdivide a pixel into 5x5 grid with 25 sensor points.
			// perform pointInPolygon test to these 25 points and calculate
			// the pixel color as linear combination of fore- & back-ground colors.
			int num_inside = 0;
			for (float dx = 0.1; dx < 1; dx += 0.2) {
				for (float dy = 0.1; dy < 1; dy += 0.2) {
					if (_pointInPolygon(x + dx, y + dy, new_polygon))
						num_inside++;
				}
			}
			patch[(y + radius - yc) * (2 * radius) + (x + radius - xc)]
				= _linearComb(fg_color, bg_color, num_inside / 25.0);
		}
	}
}

void Renderer::_bicubicInterp( rgb *hi_patch, int hi_side, rgb *lo_patch, int lo_side )
{
	// FIXME: i did not implement a bicubic interpolation.
	// i just do ad-hoc for the case where UP_SCALE = 3.
	xrgb c;
	for (int x = 0; x < 3; x++) 
		for (int y = 0; y < 3; y++)
			c += hi_patch[y * 6 + x];
	lo_patch[0] = c.simplify();
	c.clear();

	for (int x = 3; x < 6; x++) 
		for (int y = 0; y < 3; y++)
			c += hi_patch[y * 6 + x];
	lo_patch[1] = c.simplify();
	c.clear();

	for (int x = 0; x < 3; x++) 
		for (int y = 3; y < 6; y++)
			c += hi_patch[y * 6 + x];
	lo_patch[2] = c.simplify();
	c.clear();

	for (int x = 3; x < 6; x++) 
		for (int y = 3; y < 6; y++)
			c += hi_patch[y * 6 + x];
	lo_patch[3] = c.simplify();
	c.clear();
}




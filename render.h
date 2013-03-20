#pragma once

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "image.h"
#include "spline.h"
#include "TableGen/table.h"

#include "ANN.h"


class Renderer 
{
private:
	#define UP_SCALE			3
	#define PATCH_RADIUS		UP_SCALE
	#define HIRES_PATCH_RADIUS	UP_SCALE
	#define LORES_PATCH_RADIUS	1
	#define NUM_NN				300

public:
	Renderer(std::vector<float>&			nodes,
		std::vector<std::vector<int>>&		indices,
		std::vector<bool>&					junction_map,
		image<rgb>*							im,
		image<rgb>*							seg_img,
		std::map<int, rgb>&					color_map)
	{
		//SVGWriter writer("T.svg", "T.png");
		//for (int ids = 0; ids < indices.size(); ids++) {
		//	std::vector<float> polygon;
		//	for (int idn = 0; idn < indices[ids].size(); idn++) {
		//		int idx = indices[ids][idn];
		//		polygon.push_back(nodes[2 * idx]);
		//		polygon.push_back(nodes[2 * idx + 1]);
		//	}
		//	writer.writePolygonWithText(polygon, "#0000FF");
		//}
		//writer.close();
		//return;



		SVGWriter writer("T_result.svg", "T_result.bmp");
		const int num_cols = 20;
		image<rgb> *bmp = new image<rgb>(10 + 8 * num_cols, 18 + 14 * (NUM_NN / num_cols));

		std::fstream err_file("nn_dist.txt", std::ios::out);
		// construct kd-tree.
		TableGenerator gen;
		//gen.loadPatchTableBinary("patch_table.dat");
		//gen.loadPolygonTableBinary("polygon_table.dat");
		gen.exec();
		gen.savePatchTableBinary("30_patch_table.dat");
		gen.savePolygonTableBinary("30_polygon_table.dat");
		std::vector<Entry>& table = gen.table;

		int num_points = table.size();
		ANNcoord *tree_data = new ANNcoord[4 * num_points];
		for (int i = 0; i < table.size(); i++) {
			tree_data[4 * i + 0] = table[i].lo_alphas[0];
			tree_data[4 * i + 1] = table[i].lo_alphas[1];
			tree_data[4 * i + 2] = table[i].lo_alphas[2];
			tree_data[4 * i + 3] = table[i].lo_alphas[3];
		}

		ANNpointArray pa = new ANNpoint[num_points];
		for (int i = 0; i < num_points; i++) {
			pa[i] = tree_data + 4 * i;
		}

		ANNkd_tree kd_tree(pa, num_points, 4);



		int nn_idx[NUM_NN];
		double nn_dist[NUM_NN];


		for (int ids = 0; ids < indices.size(); ids++) {
			for (int idn = 0; idn < indices[ids].size(); idn++) {

				if (ids != 1 || idn != 2) {
					continue;
				}

				// get the 2x2 original patch for this node.
				int idx = indices[ids][idn];
				int xc = 0.5 + nodes[2 * idx];
				int	yc = 0.5 + nodes[2 * idx + 1];

				rgb fg_color = color_map[ids];
				rgb bg_color = _calcBgColor(xc, yc, seg_img, fg_color);

				rgb origin_patch[4];
				// FIXME check index out of range.
				origin_patch[0] = im->get(xc - 1, yc - 1);
				origin_patch[1] = im->get(xc, yc - 1);
				origin_patch[2] = im->get(xc - 1, yc);
				origin_patch[3] = im->get(xc, yc);
				//for (int y = yc - 1; y <= yc; y++) {
				//	for (int x = xc - 1; x <= xc; x++) {
				//		if (0 <= x && x < im->width() && 0 <= y && y < im->height()) {
				//			origin_patch[y * 2 + x] = im->get(x, y);
				//		}
				//	}
				//}


				// calculate the alpha value of this patch.
				double origin_alpha[4];
				_recoverAlphaPatch(origin_alpha, origin_patch, 4, fg_color, bg_color);


				// search kd-tree for the nearest neighbors.
				kd_tree.annkSearch(origin_alpha, NUM_NN, nn_idx, nn_dist);


				// calculate the 6x6 patches by the nns.
				rgb hi_patch[L * L];
				rgb lo_patch[4];
				
				_drawPatch(bmp, 1, 9, origin_patch, 2, 2);
				for (int i = 0; i < NUM_NN; i++) {
					_blendPatch(hi_patch, table[nn_idx[i]].hi_alphas, fg_color, bg_color, L * L);
					_blendPatch(lo_patch, table[nn_idx[i]].lo_alphas, fg_color, bg_color, 4);
					
					int ulx = 5 +  8 * (i % num_cols);
					int uly = 1 + 14 * (i / num_cols);

					_drawPatch(bmp, ulx, uly, hi_patch, L, L);
					_drawPatch(bmp, ulx, uly + 9, lo_patch, 2, 2);
					
					Polyitem& item = gen.polygon_table[nn_idx[i] / 2];
					std::vector<float> polygon(2 * item.num);
					for (int j = 0; j < item.num; j++) {
						polygon[2 * j]     = ulx + item.polygon[2 * j];
						polygon[2 * j + 1] = uly + item.polygon[2 * j + 1];
					}
					
					//memcpy(&polygon[0], &(item.polygon[0]), 2 * item.num * sizeof(float));
				
					writer.writePolygon(polygon, "#0000FF");
					err_file << nn_dist[i] << "\n";
				}

			}
		}
		err_file.close();
		gen._saveBMP("T_result.bmp", bmp);
		writer.close();
		delete tree_data;
		delete bmp;
	}

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

	void _recoverAlphaPatch( double *alpha, rgb *pixel, int num_pixels, rgb fg_color, rgb bg_color )
	{
		for (int i = 0; i < num_pixels; i++) {
			int cnt = 0;
			alpha[i] = 0;
			if (fg_color.r != bg_color.r) {
				alpha[i] += ((double)pixel[i].r - bg_color.r) / ((double)fg_color.r - bg_color.r);
				cnt++;
			}
			if (fg_color.g != bg_color.g) {
				alpha[i] += ((double)pixel[i].g - bg_color.g) / ((double)fg_color.g - bg_color.g);
				cnt++;
			}
			if (fg_color.b != bg_color.b) {
				alpha[i] += ((double)pixel[i].b - bg_color.b) / ((double)fg_color.b - bg_color.b);
				cnt++;
			}
			if (cnt > 0) {
				alpha[i] /= cnt;
			}
		}
	}

	void _drawPatch(image<rgb> *im, int xc, int yc, rgb *origin_patch, int w, int h)
	{
		assert(0 <= xc && xc + w <= im->width() && 0 <= yc && yc + h <= im->height());
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				im->set(xc + x, yc + y, origin_patch[y * w + x]);
			}
		}	
	}

	void _blendPatch(rgb *pixel, float *alpha, rgb fg_color, rgb bg_color, int num_pixels)
	{
		for (int i = 0; i < num_pixels; i++) {
			pixel[i].r = (int)(alpha[i] * fg_color.r + (1 - alpha[i]) * bg_color.r);
			pixel[i].g = (int)(alpha[i] * fg_color.g + (1 - alpha[i]) * bg_color.g);
			pixel[i].b = (int)(alpha[i] * fg_color.b + (1 - alpha[i]) * bg_color.b);
		}
	}
};
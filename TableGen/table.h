#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <list>
#include "gpc.h"

typedef struct
{
	float lo_alphas[4];
	float hi_alphas[L * L];
} Entry;

class TableGenerator
{
private:
	#define L			6
	#define PI			3.1415926535897932
	#define PADDING		0.2

	std::list<Entry>	table;
	std::vector<std::pair<float, float>> rot_markers;

public:
	TableGenerator()
	{
		rot_markers.push_back(std::pair<float, float>(0, 0));
		rot_markers.push_back(std::pair<float, float>(0, L));
		rot_markers.push_back(std::pair<float, float>(L, L));
		rot_markers.push_back(std::pair<float, float>(L, 0));
	}

	void generate() 
	{
		for (float d1 = -28 * 0.3; d1 <= 28 * 0.3; d1 += 0.3) {
			for (float theta1 = 0.0; theta1 <= 175.0; theta1 += 5.0) {
				for (float d2 = -28 * 0.3; d2 <= 28 * 0.3; d2 += 0.3) {
					for (float theta2 = 0.0; theta2 <= 175.0; theta2 += 5.0) {

						// eliminate parallel lines.
						if (theta1 == theta2 && d1 != d2) {
							continue;
						}

						std::vector<std::pair<float, float>> 
							marker_set1,
							marker_set2, 
							tmp;
						std::pair<float, float> intersec(-1, -1);
						// case 1: the two lines coincide.
						if (d1 == d2 && theta1 == theta2) {
							if (d1 == 0 && theta1 > 90)
								continue;
							_calcMarkersFromLine(d1, theta1, tmp);
							if (tmp.size() != 2) 
								continue;
							marker_set1.push_back(tmp[0]);
							marker_set2.push_back(tmp[1]);
							
						}
						// case 2: the two lines intersect.
						else 
						{
							intersec = _calcIntersection(d1, theta1, d2, theta2);
							if (!(0 + PADDING < intersec.first && intersec.first < L - PADDING
								&& 0 + PADDING < intersec.second && intersec.second < L - PADDING)) {
									continue;
							}
							_calcMarkersFromLine(d1, theta1, marker_set1);
							_calcMarkersFromLine(d2, theta2, marker_set2);
						}

						

						
			
						
						for (int idx1 = 0; idx1 < marker_set1.size(); idx1++) {
							for (int idx2 = 0; idx2 < marker_set2.size(); idx2++) {
								std::pair<float, float> start = marker_set1[idx1];
								std::pair<float, float>   end = marker_set2[idx2];
								std::vector<float> polygon;
								Entry entry = {0};

								_getFgPolygon(start, end, intersec, polygon);	

								for (int x = 0; x < L; x++) {
									for (int y = 0; y < L; y++) {
										entry.hi_alphas[y * L + x] = _calcPixelAlpha(x, y, polygon);
										entry.lo_alphas[(y/3) * 2 + (x/3)] += entry.hi_alphas[y * L + x];
									}
								}
								for (int i = 0; i < 4; i++) {
									entry.lo_alphas[i] /= ((L/2) * (L/2));
								}

								table.push_back(entry);
								for (int i = 0; i < 4; i++) 
									entry.lo_alphas[i] = 1 - entry.lo_alphas[i];
								for (int i = 0; i < L * L; i++)
									entry.hi_alphas[i] = 1 - entry.hi_alphas[i];
								table.push_back(entry);
							}
						}
					}
				}
			}
		}
	}

	void _getFgPolygon( std::pair<float, float>&	start, 
						std::pair<float, float>&	end, 
						std::pair<float, float>&	intersec, 
						std::vector<float>&			polygon ) 
	{
		polygon.push_back(start.first);
		polygon.push_back(start.second);

		int rot_id_start = _getRotationID(start);
		int rot_id_end   = _getRotationID(end);
		for (int rot_id = rot_id_start; rot_id != rot_id_end; 
			rot_id = (rot_id + 1) % 4) {
			polygon.push_back(rot_markers[rot_id].first);
			polygon.push_back(rot_markers[rot_id].second);
		}

		polygon.push_back(end.first);
		polygon.push_back(end.second);

		if (intersec.first > 0 && intersec.second > 0) {
			polygon.push_back(intersec.first);
			polygon.push_back(intersec.second);
		}
	}

	void _calcMarkersFromLine( float d, float theta, std::vector<std::pair<float, float>>& markers ) 
	{
		if (theta == 90) {
			if (d < 0 || d > L)
				return;
			markers.push_back(std::pair<float, float>(d, 0));
			markers.push_back(std::pair<float, float>(d, L));
			return;
		}
		if (theta == 0) {
			if (d > 0 || d < -L)
				return;
			markers.push_back(std::pair<float, float>(0, -d));
			markers.push_back(std::pair<float, float>(L, -d));
		}
		if (d == 0 && theta > 90)
			return;

		double t = theta / PI,
			   k = std::tan(t),
			   x0 = d * std::sin(t),
			   y0 = -d * std::cos(t);
		double xu = (0 - y0) / k + x0,
			   xd = (L - y0) / k + x0,
			   yl = k * (0 - x0) + y0,
			   yr = k * (L - x0) + y0;

		if (0 <= xu && xu <= L)
			markers.push_back(std::pair<float, float>(xu, 0));
		if (0 <= xd && xd <= L)
			markers.push_back(std::pair<float, float>(xd, L));
		if (0 <  yr && yr <  L)
			markers.push_back(std::pair<float, float>(L, yr));
		if (0 <  yl && yl <  L)
			markers.push_back(std::pair<float, float>(0, yl));
	}

	std::pair<float, float> _calcIntersection( float d1, float theta1, float d2, float theta2 ) 
	{
		// this method assumes that the two lines from the input must intersect.
		double t1 = theta1 / PI,
			t2 = theta2 / PI,
			x1 = d1 * std::sin(t1),
			x2 = d2 * std::sin(t2),
			y1 = -d1 * std::cos(t1),
			y2 = -d2 * std::cos(t2),
			k1 = std::tan(t1),
			k2 = std::tan(t2),
			x, y;

		bool vertical1 = (theta1 == 90),
			 vertical2 = (theta2 == 90);

		if (vertical1 || vertical2) {
			if (vertical1) {
				std::swap(x1, x2);
				std::swap(y1, y2);
				std::swap(k1, k2);
			}
			x = x2;
			y = y1 + k1 * (x - x1);
		}
		else {
			x = (y2 - y1 + k1*x1 - k2*x2) / (k1 - k2);
			y = y1 + k1 * (x - x1);
		}

		return std::pair<float, float>(x, y);
	}

	int _getRotationID( std::pair<float, float>& marker ) 
	{
		float& x = marker.first;
		float& y = marker.second;
		if (y == 0 && 0 < x && x <= L)
			return 0;
		if (x == 0 && 0 <= y && y < L)
			return 1;
		if (y == L && 0 <= x && x < L)
			return 2;
		if (x == L && 0 < y && y <= L)
			return 3;
		return 0;  // control should never reach here.
	}

	double _calcPixelAlpha( int x, int y, std::vector<float>& polygon ) 
	{
		gpc_polygon pixel_polygon, corner_polygon, result;
		gpc_vertex_list pixel_contour, corner_contour;

		gpc_vertex pixel_vertices[4];
		pixel_vertices[0].x = x + 0;	pixel_vertices[0].y = y + 0;
		pixel_vertices[1].x = x + 1;	pixel_vertices[1].y = y + 0;
		pixel_vertices[2].x = x + 1;	pixel_vertices[2].y = y + 1;
		pixel_vertices[3].x = x + 0;	pixel_vertices[3].y = y + 1;

		gpc_vertex corner_vertices[7];  // there are at most 7 markers in the polygon.
		for (int i = 0; i < polygon.size() / 2; i++) {
			corner_vertices[i].x = polygon[2 * i];
			corner_vertices[i].y = polygon[2 * i + 1];
		}
		
		pixel_contour.num_vertices =4;
		pixel_contour.vertex = pixel_vertices;
		corner_contour.num_vertices = polygon.size() / 2;
		corner_contour.vertex = corner_vertices;

		pixel_polygon.contour = &pixel_contour;
		pixel_polygon.hole = NULL;
		pixel_polygon.num_contours = 1;
		corner_polygon.contour = &corner_contour;
		corner_polygon.hole = NULL;
		corner_polygon.num_contours = 1;

		gpc_polygon_clip(GPC_INT, &pixel_polygon, &corner_polygon, &result);
		double area = 0;
		if (result.num_contours > 0) {
			int n = result.contour[0].num_vertices;
			gpc_vertex *v = result.contour[0].vertex;
			for (int i = 0; i < n; i++) {
				area += v[i].x * v[(i + 1) % n].y - v[(i + 1) % n].x - v[i].y;
			}
			area /= 2;
		}
		gpc_free_polygon(&result);
		
		return area;
	}


};
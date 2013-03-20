#include "table.h"



double TableGenerator::_calcPixelAlpha( int x, int y, std::vector<float>& polygon )
{
	gpc_polygon		pixel_polygon, corner_polygon, result;
	gpc_vertex_list pixel_contour, corner_contour;

	gpc_vertex pixel_vertices[4];
	pixel_vertices[0].x = x + 0;	pixel_vertices[0].y = y + 0;
	pixel_vertices[1].x = x + 1;	pixel_vertices[1].y = y + 0;
	pixel_vertices[2].x = x + 1;	pixel_vertices[2].y = y + 1;
	pixel_vertices[3].x = x + 0;	pixel_vertices[3].y = y + 1;

	gpc_vertex corner_vertices[7];		// there are at most 7 markers in the polygon.
	assert(polygon.size() <= 2 * 7);	
	for (int i = 0; i < polygon.size() / 2; i++) {
		corner_vertices[i].x = polygon[2 * i];
		corner_vertices[i].y = polygon[2 * i + 1];
	}

	pixel_contour.num_vertices =4;
	pixel_contour.vertex = pixel_vertices;
	pixel_polygon.contour = &pixel_contour;
	pixel_polygon.hole = NULL;
	pixel_polygon.num_contours = 1;

	corner_contour.num_vertices = polygon.size() / 2;
	corner_contour.vertex = corner_vertices;
	corner_polygon.contour = &corner_contour;
	corner_polygon.hole = NULL;
	corner_polygon.num_contours = 1;

	gpc_polygon_clip(GPC_INT, &pixel_polygon, &corner_polygon, &result);
	double area = 0;
	if (result.num_contours > 0) {
		int n = result.contour[0].num_vertices;
		gpc_vertex *v = result.contour[0].vertex;
		for (int i = 0; i < n; i++) {
			area += v[i].x * v[(i + 1) % n].y - v[(i + 1) % n].x * v[i].y;
		}
		area /= 2;
	}
	gpc_free_polygon(&result);

	return std::abs(area);
}

int TableGenerator::_getRotationID( std::pair<float, float>& marker )
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
	return 0;			// control should never reach here.
}

void TableGenerator::_getFgPolygon( std::pair<float, float>& start, std::pair<float, float>& end, std::pair<float, float>& intersec, std::vector<float>& polygon )
{
	int rot_id_start = _getRotationID(start);
	int rot_id_end   = _getRotationID(end);

	polygon.push_back(start.first);
	polygon.push_back(start.second);

	for (int rot_id = rot_id_start; rot_id != rot_id_end; rot_id = (rot_id + 1) % 4) {
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

std::pair<float, float> TableGenerator::_calcIntersection( float d1, float theta1, float d2, float theta2 )
{
	// this method assumes that the two lines from the input must intersect.
	double  t1 = PI * (theta1 / 180.0),
		t2 = PI * (theta2 / 180.0),
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

void TableGenerator::_calcMarkersFromLine( float d, float theta, std::vector<std::pair<float, float>>& markers )
{
	if (theta == 90) {					// vertical case.
		if (d < 0 || d > L)
			return;
		markers.push_back(std::pair<float, float>(d, 0));
		markers.push_back(std::pair<float, float>(d, L));
		return;
	}

	if (theta == 0) {					// horizontal case.
		if (-d < 0 || -d > L)
			return;
		markers.push_back(std::pair<float, float>(0, -d));
		markers.push_back(std::pair<float, float>(L, -d));
		return;
	}

	if (d == 0 && theta > 90) {			// case that only intersect at the origin.
		return;
	}

	double t = PI * (theta / 180.0),	// case that have 2 intersections.
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

void TableGenerator::exec()
{
	const int num_samples = 30;
	double sample_gap = L * std::sqrt(2.0) / num_samples;

	//for (float d1 = -num_samples * sample_gap; d1 <= num_samples * sample_gap; d1 += sample_gap) {
	//	for (float theta1 = 0.0; theta1 <= 175.0; theta1 += 5.0) {

	//		for (float d2 = /*notice*/d1; d2 <= num_samples * sample_gap; d2 += sample_gap) {
	//			for (float theta2 = /*notice*/(d1 == d2 ? theta1 : 0.0); theta2 <= 175.0; theta2 += 5.0) {
	for (float theta1 = 0.0; theta1 <= 175.0; theta1 += 5.0) {
		for (float d1 = -num_samples * sample_gap; d1 <= num_samples * sample_gap; d1 += sample_gap) {
		
			for (float theta2 = /*notice*/theta1; theta2 <= 175.0; theta2 += 5.0) {
				for (float d2 = /*notice*/(theta1 == theta2 ? d1 : -num_samples * sample_gap); d2 <= num_samples * sample_gap; d2 += sample_gap) {
				

					// eliminate parallel lines.
					if (theta1 == theta2 && d1 != d2) {
						continue;
					}

					std::vector<std::pair<float, float>> marker_set1, marker_set2, tmp;
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
						if (marker_set1.size() != 2 || marker_set2.size() != 2)
							continue;
					}


					//std::cout << table.size() << std::endl;

					// in case 1: a line contributes 2 patches.
					// in case 2: a line pair contributes 8 patches.						
					for (int idx1 = 0; idx1 < marker_set1.size(); idx1++) {
						for (int idx2 = 0; idx2 < marker_set2.size(); idx2++) {
							std::pair<float, float> start = marker_set1[idx1];
							std::pair<float, float>   end = marker_set2[idx2];
							std::vector<float> polygon;
							Entry entry = {0};

							_getFgPolygon(start, end, intersec, polygon);	

							for (int y = 0; y < L; y++) {
								for (int x = 0; x < L; x++) {
									entry.hi_alphas[y * L + x] = _calcPixelAlpha(x, y, polygon);
									entry.lo_alphas[(y / UP_SCALE) * 2 + (x / UP_SCALE)] += entry.hi_alphas[y * L + x];
								}
							}
							for (int i = 0; i < 4; i++) {
								entry.lo_alphas[i] /= (UP_SCALE * UP_SCALE);
							}

							table.push_back(entry);
							for (int i = 0; i < 4; i++) 
								entry.lo_alphas[i] = 1 - entry.lo_alphas[i];
							for (int i = 0; i < L * L; i++)
								entry.hi_alphas[i] = 1 - entry.hi_alphas[i];
							table.push_back(entry);

							polygon_table.push_back(Polyitem(polygon.size() / 2, &polygon[0]));
						}
					}
				}
			}
		}
	}

	std::cout << "num_entries: " << table.size() << std::endl;
}

TableGenerator::TableGenerator()
{
	table.reserve(200000);
	polygon_table.reserve(200000);
	rot_markers.push_back(std::pair<float, float>(0, 0));
	rot_markers.push_back(std::pair<float, float>(0, L));
	rot_markers.push_back(std::pair<float, float>(L, L));
	rot_markers.push_back(std::pair<float, float>(L, 0));
}


#pragma once
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include "image.h"
#include "helper.h"
#include <Eigen/Dense>



struct ConicSection
{
	double A, B, C, D, E, F;
	float minx, miny, maxx, maxy;
};

class CurveFitter
{
#define NUM_NEIGHBORS 15
private:
	std::vector<ConicSection> conics;

public:
	CurveFitter(std::vector<float>&				nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map)
	{
		// fit curves.
		for (int ids = 0; ids < indices.size(); ids++) {
			for (int idn = 0; idn < indices[ids].size(); idn += 10) {
				std::vector<float> points(2 * (1 + 2 * NUM_NEIGHBORS));
				int shape_size = indices[ids].size();
				float minx = 99999, miny = minx;
				float maxx = -1, maxy = maxx;

				for (int i = idn - NUM_NEIGHBORS; i <= i + NUM_NEIGHBORS; i++) {
					// in the case shape_size is much smaller than NUM_NEIGHBORS, 
					// i + shape_size may still be negative, 
					// so we use shape_size * NUM_NEIGHBORS instead.
					int idx = indices[ids][(i + shape_size * NUM_NEIGHBORS) % shape_size];
					points.push_back(nodes[2 * idx]);
					points.push_back(nodes[2 * idx + 1]);

					minx = std::min(minx, points[2 * idx]);
					maxx = std::max(maxx, points[2 * idx]);
					miny = std::min(miny, points[2 * idx + 1]);
					maxy = std::max(maxy, points[2 * idx + 1]);
				}

				ConicSection tmp;
				tmp.minx = minx;
				tmp.miny = miny;
				tmp.maxx = maxx;
				tmp.maxy = maxy;
				_fitConic(points, tmp);
				conics.push_back(tmp);
			}
		}
	}


	void saveToSVG(std::string& filename)
	{
		std::ofstream file(filename, std::ios::out | std::ios::binary);
		for (std::vector<ConicSection>::const_iterator it = conics.begin();
			 it != conics.end(); it++) {
				 // densely sample each conic section
				 // and ouput the sampled points to file.
		}
	}


	void _fitConic(std::vector<float>& points, ConicSection& conic)
	{
		const int n = 4;//1 + 2 * NUM_NEIGHBORS;  // must equal to points.size() / 2.
		Eigen::MatrixXd    b(n, 1);  // bi = -(yi)^2
		Eigen::MatrixXd    A(n, 5);  // A  = [a1, a2, ..., an]^T
									 // ai = [xi^2-yi^2, 2xiyi, 2xi, 2yi, 1]^T
		Eigen::MatrixXd   AT(5, n);
		Eigen::MatrixXd  ATA(5, 5);
		Eigen::MatrixXd  ATb(5, 1);
		Eigen::MatrixXd    p(5, 1);  // p = [A, B, D, E, F]^T
									 // p = (A^T * A)^{-1} * A^T * b
		double xi, yi;

		for (int i = 0; i < n; i++) {
			yi = points[2 * i + 1];
			b(i, 0) = -(yi * yi);
		}

		for (int i = 0; i < n; i++) {
			xi = points[2 * i];
			yi = points[2 * i + 1];
			A(i, 0) = xi * xi - yi * yi;	
			A(i, 1) = 2 * xi * yi;
			A(i, 2) = 2 * xi;
			A(i, 3) = 2 * yi;
			A(i, 4) = 1;
		}

		AT  = A.transpose();
		ATA = AT * A;
		ATb = AT * b;
		p   = ATA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		conic.A = p(0,0);
		conic.B = p(1,0);
		conic.C = 1 - conic.A;
		conic.D = p(2,0);
		conic.E = p(3,0);
		conic.F = p(4,0);
	}

};
#pragma once
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include "image.h"
#include "helper.h"
#include <Eigen/Dense>
#include <iostream>
#include <cmath>



struct ConicSection
{
	double A, B, C, D, E, F;
	float minx, miny, maxx, maxy;
	float lx, ly, ux, uy;
};

class ConicFitter
{
#define NUM_NEIGHBORS 15
private:
	std::vector<ConicSection> conics;

public:
	ConicFitter(std::vector<float>&				nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map)
	{
		for (int ids = 0; ids < indices.size(); ids++) {
			for (int idn = 0; idn < indices[ids].size(); idn += 10) {

				std::vector<float> points(2 * (1 + 2 * NUM_NEIGHBORS));
				int shape_size = indices[ids].size();
				float minx = 99999, miny = minx;
				float maxx = -1, maxy = maxx;
				float lx, ly, ux, uy;

				for (int i = idn - NUM_NEIGHBORS; i <= idn + NUM_NEIGHBORS; i++) {
					// in the case shape_size is much smaller than NUM_NEIGHBORS, 
					// i + shape_size may still be negative, 
					// so we use shape_size * NUM_NEIGHBORS instead.
					int idx = indices[ids][(i + shape_size * NUM_NEIGHBORS) % shape_size];
					int px = 2 * (i + NUM_NEIGHBORS - idn);
					int py = px + 1;

					points[px] = nodes[2 * idx];
					points[py] = nodes[2 * idx + 1];

					if (points[px] < minx) {
						minx = points[px];
						lx = points[px];
						ly = points[py];
					}
					if (points[py] < miny) {
						miny = points[py];
						ux = points[px];
						uy = points[py];
					}
					// minx = std::min(minx, points[px]);
					maxx = std::max(maxx, points[px]);
					// miny = std::min(miny, points[py]);
					maxy = std::max(maxy, points[py]);
				}

				ConicSection tmp;
				tmp.minx = minx;
				tmp.miny = miny;
				tmp.maxx = maxx;
				tmp.maxy = maxy;
				tmp.lx = lx;
				tmp.ly = ly; 
				tmp.ux = ux;
				tmp.uy = uy;
				_fitConic(points, tmp);
				conics.push_back(tmp);
			}
		}
		_saveConics();
	}

	void _saveConics()
	{
		std::ofstream file("conics.txt", std::ios::out);
		for (int i = 0; i < conics.size(); i++) {
			file << conics[i].A << "\t"
				 << conics[i].B << "\t"
				 << conics[i].C << "\t"
				 << conics[i].D << "\t"
				 << conics[i].E << "\t"
				 << conics[i].F << std::endl;
			file << "minx: " << conics[i].minx << std::endl
				 << "maxx: " << conics[i].maxx << std::endl
				 << "miny: " << conics[i].miny << std::endl
				 << "maxy: " << conics[i].maxy << std::endl;
		}
		for (int i = 0; i < conics.size() ; i++) {
			file << "ezplot('("
				 << conics[i].A << ")*x*x + ("
				 << conics[i].B << ")*x*y + ("
				 << conics[i].C << ")*y*y + ("
				 << conics[i].D << ")*x + ("
				 << conics[i].E << ")*y + ("
				 << conics[i].F << ") = 0')"
				 << std::endl;
		}

		file.close();
	}

	void saveToSVG(std::string& filename)
	{
		std::ofstream file(filename, std::ios::out | std::ios::binary);
		const std::string svg_header = "<?xml version=\"1.0\" standalone=\"no\"?>\n"
			"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
			"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
			"<svg width=\"224\" height=\"225\" version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns=\"http://www.w3.org/2000/svg\">\n";
		file << svg_header << std::endl;

		std::ofstream file2("sample_size.txt", std::ios::out);

		for (int i = 0; i < conics.size(); i++) {
			// densely sample each conic section
			// and ouput the sampled points to file.
			std::vector<float> samples;
			_sampleConicSection(conics[i], samples);
			file2 << samples.size() << std::endl;

			file << "<polyline points=\"";
			for (int j = 0; j < samples.size(); j += 2) {
				file << samples[j] << ","
					 << samples[j + 1] << " ";
			}
			file << "\" style=\"fill-opacity:0; "
				 << "stroke-opacity:0.8; "
				 << "stroke:#00ff00; stroke-width:0.2\"/>"
				 << std::endl;
		}
		file << "</svg>" << std::endl;
		file.close();
		file2.close();
	}


	void _sampleConicSection(ConicSection& conic, std::vector<float>& samples)
	{
		enum {FOREACH_X, FOREACH_Y};
		int mode;
		if ((conic.maxx - conic.minx) > (conic.maxy - conic.miny)) {
			mode = FOREACH_X;
			samples.push_back(conic.lx);
			samples.push_back(conic.ly);
		}
		else {
			mode = FOREACH_Y;
			samples.push_back(conic.ux);
			samples.push_back(conic.uy);
		}

		const double EPSILON = 0.1;
		double a, b, c;
		double A = conic.A, 
			B = conic.B,
			C = conic.C,
			D = conic.D,
			E = conic.E,
			F = conic.F;
			
		if (mode == FOREACH_X) {
			for (float x = conic.minx; x <= conic.maxx; x += 0.1f) {
				a = C;
				b = B*x + E;
				c = A*x*x + D*x + F;
				std::vector<float> y = _solveQuadratic(a, b, c);
				if (y.size() == 1) {
					samples.push_back(x);
					samples.push_back(y[0]);
				}
				else if (y.size() == 2) {
					float yprev = samples[samples.size() - 1];
					samples.push_back(x);
					if (std::abs(y[0] - yprev) < std::abs(y[1] - yprev))
						samples.push_back(y[0]);
					else
						samples.push_back(y[1]);
				}
			}
		}
		else {
			for (float y = conic.miny; y <= conic.maxy; y += 0.1f) {
				a = A;
				b = B*y + D;
				c = C*y*y + E*y + F;
				std::vector<float> x = _solveQuadratic(a, b, c);
				if (x.size() == 1) {
					samples.push_back(x[0]);
					samples.push_back(y);
				}
				else if (x.size() == 2) {
					float xprev = samples[samples.size() - 2];
					if (std::abs(x[0] - xprev) < std::abs(x[1] - xprev))
						samples.push_back(x[0]);
					else
						samples.push_back(x[1]);
					samples.push_back(y);
				}
			}
		}
		samples.erase(samples.begin());
		samples.erase(samples.begin());
	}


	std::vector<float> _solveQuadratic(double a, double b, double c)
	{
		const double eps = 1.0e-6;  // i wonder if there are some conventional epsilon.
		std::vector<float> roots;

		if (std::abs(a) < eps) {  // a == 0.
			if (std::abs(b) > eps) {
				roots.push_back(-c / b);
			}
			return roots;
		}

		double delta = b*b - 4*a*c;
		if (delta < -eps) {  // delta < 0.
			//std::cout << "zeros!!!" << std::endl;
			return roots;
		}

		if (delta < eps) {  // treate it as 0.
			roots.push_back(-b / (2*a));
			return roots;
		}

		double sqrt_delta = std::sqrt(delta);
		roots.push_back((-b - sqrt_delta) / (2*a));
		roots.push_back((-b + sqrt_delta) / (2*a));
		return roots;
	}


	void _fitConic(std::vector<float>& points, ConicSection& conic)
	{
		const int n = 1 + 2 * NUM_NEIGHBORS;  // must equal to points.size() / 2.
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
		 p   = ATA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ATb);
		//p = ATA.inverse() * ATb;
#ifdef MY_DEBUG
		std::cout << ATA << std::endl;
		std::cout << ATb << std::endl;
		std::cout << p	 << std::endl;
#endif
		conic.A = p(0, 0);
		conic.B = p(1, 0) * 2;
		conic.C = 1 - conic.A;
		conic.D = p(2, 0) * 2;
		conic.E = p(3, 0) * 2;
		conic.F = p(4, 0);
	}

};
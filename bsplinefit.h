#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "b-spline\f2c.h"


extern "C" {
/* Subroutine */ int curev_(integer *idim, real *t, integer *n, real *c__, 
	integer *nc, integer *k, real *u, integer *m, real *x, integer *mx, 
	integer *ier);
/* Subroutine */ int parcur_(integer *iopt, integer *ipar, integer *idim, 
	integer *m, real *u, integer *mx, real *x, real *w, real *ub, real *
	ue, integer *k, real *s, integer *nest, integer *n, real *t, integer *
	nc, real *c__, real *fp, real *wrk, integer *lwrk, integer *iwrk, 
	integer *ier);
}


struct BSpline
{
	int k;
	std::vector<float> t;
	std::vector<float> c;
	std::vector<float> u;
};

class BSplineFitter
{
#define NUM_NEIGHBORS 15
private:
	std::vector<BSpline> splines;

public:
	BSplineFitter(std::vector<float>&			nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map)
	{
		for (int ids = 0; ids < indices.size(); ids++) {
			for (int idn = 0; idn < indices[ids].size(); idn += 10) {

				std::vector<float> points(2 * (1 + 2 * NUM_NEIGHBORS));
				int shape_size = indices[ids].size();

				for (int i = idn - NUM_NEIGHBORS; i <= idn + NUM_NEIGHBORS; i++) {
					// in the case shape_size is much smaller than NUM_NEIGHBORS, 
					// i + shape_size may still be negative, 
					// so we use shape_size * NUM_NEIGHBORS instead.
					int idx = indices[ids][(i + shape_size * NUM_NEIGHBORS) % shape_size];
					int px = 2 * (i + NUM_NEIGHBORS - idn);
					int py = px + 1;

					points[px] = nodes[2 * idx];
					points[py] = nodes[2 * idx + 1];
				}
				
				splines.push_back(BSpline());
				_fitBSpline(points, splines[splines.size() - 1]);
			}
		}
		_saveBSpline();
	}


	void _fitBSpline(std::vector<float>& points, BSpline& spline)
	{
		long  iopt = 0, ipar = 0, idim = 2, k = 3,
			  m = points.size() / 2,
			  mx = idim * m,
			  nest = m + k + 1,
			  nc = nest * idim,
			  lwrk = m*(k+1)+nest*(6+idim+3*k),
			  n, ier;
		long  *iwrk = new long[nest];
		float ub, ue, fp, s = 5,
			  *x = &points[0];
		float *u = new float[m],
			  *w = new float[m],
			  *t = new float[nest],
			  *c = new float[nc],
			  *wrk = new float[lwrk];

		for (int i = 0; i < m; i++)
			w[i] = 1;

		  
		parcur_(&iopt, &ipar, &idim,
			&m, u, &mx, x, w, &ub, &ue,
			&k, &s, &nest, &n, t, &nc,
			c, &fp, wrk, &lwrk, iwrk,
			&ier);

		spline.k = k;
		spline.t.resize(n);
		memcpy(&spline.t[0], t, n * sizeof(float));
		spline.c.resize(nc);
		memcpy(&spline.c[0], c, nc * sizeof(float));
		spline.u.resize(m);
		memcpy(&spline.u[0], u, m * sizeof(float));

		delete[] iwrk, u, w, t, c, wrk;
	}


	void _sampleBSpline(BSpline& spline, std::vector<float>& samples)
	{
		long idim = 2, 
			 n = spline.t.size(),
			 nc = spline.c.size(),
			 k = spline.k,
		     m = spline.u.size(),
			 mx = idim * m,
			 ier;
		samples.resize(2 * m);
		float *t = &spline.t[0],
			  *c = &spline.c[0],
			  *u = &spline.u[0],
			  *x = &samples[0];
			 
		curev_(&idim, t, &n, c,
		&nc, &k, u, &m, x, &mx,
		&ier);

		// by here the results are already stored in samples.
	}


	void _saveBSpline()
	{
		std::ofstream file("bsplines.txt", std::ios::out);
		for (int i = 0; i < splines.size(); i++) {
			int n = splines[i].t.size(),
				k = splines[i].k;
			file << "Degree: " << k << "\n"
				 << "Knots: "  << n << "\n";
			for (int j = 0; j < splines[i].t.size(); j++)
				file << splines[i].t[j] << ", ";
			file << "\nControl points: " << n - k - 1 << "\n";
			for (int j = 0; j < n - k - 1; j++)
				file << "(" << splines[i].c[n * 0 + j]
			         << "," << splines[i].c[n * 1 + j]
					 << "), ";
			file << std::endl;
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
		std::cout <<  file.is_open() << std::endl;

		for (int i = 0; i < splines.size(); i++) {
			// densely sample each conic section
			// and ouput the sampled points to file.
			std::vector<float> samples;
			_sampleBSpline(splines[i], samples);

			file << "<polyline points=\"";
			for (int j = 0; j + 1 < samples.size(); j += 2) {
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
	}
};
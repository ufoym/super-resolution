#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "b-spline\f2c.h"
#include "writesvg.h"


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
#define NUM_NEIGHBORS		10
#define SAMPLE_INTERVAL		5
private:
	std::vector<BSpline> splines;

public:
	BSplineFitter(std::vector<float>&			nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map)
	{
		for (int ids = 0; ids < indices.size(); ids++) {
			for (int idn = 0; idn < indices[ids].size(); idn += SAMPLE_INTERVAL) {

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
		long  iopt = 0, ipar = 0, idim = 2, k = 5,
			  m = points.size() / 2,
			  mx = idim * m,
			  nest = m + k + 1,
			  nc = nest * idim,
			  lwrk = m*(k+1)+nest*(6+idim+3*k),
			  n, ier;
		long  *iwrk = new long[nest];
		float ub, ue, fp, s = 3.0,
			  *x = &points[0];
		float *u = new float[m],
			  *w = new float[m],
			  *t = new float[nest],
			  *c = new float[nc],
			  *wrk = new float[lwrk];

		// use gaussian weights that have pdf=1.0 at the central,
		// and have pdf = e^(-1.5)=0.2231 at the two ends.
		double variance = NUM_NEIGHBORS * NUM_NEIGHBORS / 3;
		int idx_central = m / 2;
		for (int i = 0; i < m; i++) {
			w[i] = 1.0;
			//w[i] = std::exp(-((i - idx_central) * (i - idx_central)) / (2 * variance));
			//w[i] = 1.0 / w[i];
			/*if (i == 0)
				std::cout << w[i] << std::endl;*/
		}
		  
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
		SVGWriter writer(filename);
		for (int i = 0; i < splines.size(); i++) {
			std::vector<float> samples;
			_sampleBSpline(splines[i], samples);
			writer.writePolygon(samples, "#00FF00");
		}
		writer.close();
	}
};
#include <sstream>
#include "spline.h"
#include "svg.h"
#include "b-spline\f2c.h"  // this header should be lastly added.

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

BSplineFitter::BSplineFitter( std::vector<float>& nodes, 
							  std::vector<std::vector<int>>& indices, 
							  std::vector<bool>& junction_map )
{
	std::pair<int, float> param_table[NUM_CURVE_TYPES];
	param_table[0] = std::pair<int, float>(1, 0.0);
	param_table[1] = std::pair<int, float>(1, 0.5);
	param_table[2] = std::pair<int, float>(3, 0.5);
	param_table[3] = std::pair<int, float>(3, 2.0);
	param_table[4] = std::pair<int, float>(5, 2.0);

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

			for (int i = 0; i < NUM_CURVE_TYPES; i++) {
				int k = param_table[i].first;
				float s = param_table[i].second;

				splines.push_back(BSpline());
				splines[splines.size() - 1].id_shape = ids;
				splines[splines.size() - 1].id_node = idn;
				splines[splines.size() - 1].data = points;
				_fitBSpline(k, s, points, splines[splines.size() - 1]);
			}
		}
	}
	_saveBSpline();
}


void BSplineFitter::clear()
{
	splines.clear();
}

void BSplineFitter::getSplines( std::vector<BSpline>& splines )
{
	splines = this->splines;
}

void BSplineFitter::saveToSVG( std::string& filename )
{
	SVGWriter writer(filename);
	for (int i = 0; i < splines.size(); i++) {
		//std::vector<float> samples;
		//_sampleBSpline(splines[i], samples);
		//writer.writePolygon(samples, "#00FF00");
		std::string color;
		std::stringstream stream;
		stream << i;
		std::string text = stream.str();

		int n = splines[i].samples.size() / 2;
		float x = splines[i].samples[2 * (n / 2)];
		float y = splines[i].samples[2 * (n / 2) + 1];
		if (i % 2 == 0) {color = "#00FF00";}
		else			{color = "#0000FF";}

		writer.writePolygon(splines[i].samples, color);
		writer.writeText(x, y, text, color);
	}
	writer.close();
}

void BSplineFitter::_saveBSpline()
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

void BSplineFitter::_sampleBSpline( BSpline& spline, std::vector<float>& samples )
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

void BSplineFitter::_fitBSpline( long k, float s, std::vector<float>& points, BSpline& spline )
{
	long  iopt = 0, ipar = 0, idim = 2,
		m = points.size() / 2,
		mx = idim * m,
		nest = m + k + 1,
		nc = nest * idim,
		lwrk = m*(k+1)+nest*(6+idim+3*k),
		n, ier;
	long  *iwrk = new long[nest];
	float ub, ue, fp,
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
	spline.samples.resize(2 * m);
	_sampleBSpline(spline, spline.samples);

	delete[] iwrk, u, w, t, c, wrk;
}

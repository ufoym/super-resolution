#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include "f2c.h"



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


class SVGWriter
{
private:
	std::ofstream file;

public:
	SVGWriter(std::string filename)
	{
		file.open(filename, std::ios::out | std::ios::binary);
		file << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
			"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
			"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
			"<svg width=\"30\" height=\"30\" version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	}

	void writeDots(std::vector<float> points, std::string color)
	{
		for (int i = 0 ; i < points.size(); i += 2) {
			file << "<circle cx=\"" << points[i] 
			     << "\" cy=\"" << points[i + 1] 
				 << "\" r=\"0.1\" fill=\"" << color << "\" />"
				 << std::endl;
		}
	}

	void writePolygon(std::vector<float> points, std::string color)
	{
		file << "<polyline points=\"";
		for (int i = 0; i < points.size(); i += 2) {
			file << points[i] << ","
				 << points[i + 1] << " ";
		}
		file << "\" style=\"fill-opacity:0; "
			 << "stroke-opacity:0.8; "
			 << "stroke:" << color << "; stroke-width:0.1\"/>"
			 << std::endl;
	}

	void close()
	{
		file << "</svg>" << std::endl;
		file.close();
	}
};

int main()
{
	

	const double PI = 3.14159265;
	std::vector<float> points;

	for (float y = 0; y <= 30; y += 0.3) {
		float x = 5 * std::sin(y);
		float noise =  + (double)std::rand() / RAND_MAX;
		x += (noise > 0.5 ? 1 : -1) * noise*2;
		points.push_back(x);
		points.push_back(y);
	}
	
	points.clear();
	for (float i = 10.0; i >= 0; i -= 0.2) {
		points.push_back(0);
		points.push_back(i);
	}
	for (float i = 0.2; i <= 10.0; i += 0.2) {
		points.push_back(i);
		points.push_back(0);
	}
	

	long iopt = 0, ipar = 0, idim = 2, k = 3,
		 m = points.size() / 2,
		 mx = idim * m,
		 nest = m + k + 1,
		 nc = nest * idim,
		 lwrk = m*(k+1)+nest*(6+idim+3*k),
		 n, ier;
	long *iwrk = new long[nest];
	float ub, ue, fp, s = 0.3;
	float *u = new float[m],
		  *w = new float[m],
		  *x = new float[mx],
		  *t = new float[nest],
		  *c = new float[nc],
		  *wrk = new float[lwrk];

	for (int i = 0; i < m; i++)
		w[i] = 1;
	for (int i = 0; i < mx; i++)
		x[i] = points[i];

		  
	parcur_(&iopt, &ipar, &idim,
		&m, u, &mx, x, w, &ub, &ue,
		&k, &s, &nest, &n, t, &nc,
		c, &fp, wrk, &lwrk, iwrk,
		&ier);

	std::cout << ier << std::endl;
	std::cout << fp << std::endl;
	for (int i = 1; i < m; i++)
		std::cout << u[i] - u[i - 1] << std::endl;


	/* Subroutine */ int curev_(integer *idim, real *t, integer *n, real *c__, 
	integer *nc, integer *k, real *u, integer *m, real *x, integer *mx, 
	integer *ier);

	//u = new float[1001];
	//x = new float[2002];
	//u[0] = 0;
	//for (int i = 1; i < 1001; i++)
	//	u[i] = u[i - 1] + (1.0/1001);
	//m = 1001;
	//mx = 2002;

	//std::cout << "**" << u[1000] << std::endl;

	curev_(&idim, t, &n, c,
		&nc, &k, u, &m, x, &mx,
		&ier);

	std::cout << std::endl;
	std::cout << ier << std::endl;
	

	std::vector<float> result;
	for (int i = 0; i < 2 * m; i++)
		result.push_back(x[i]);

	SVGWriter writer("test.svg");
	writer.writeDots(points, "#FF0000");
	//writer.writePolygon(points, "#0000FF");
	writer.writePolygon(result, "#0000FF");
	writer.close();

	delete[] iwrk, u, w, x, t, c, wrk;
}
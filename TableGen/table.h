#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <list>
#include <cassert>

#include "gpc.h"
#include "../image.h"
#include "../svg.h"


#define PI			3.1415926535897932
#define PADDING		0.6
#define UP_SCALE	3
#define L			6	// (2 * UP_SCALE)


typedef struct
{
	float lo_alphas[4];
	float hi_alphas[L * L];
} Entry;

struct Polyitem
{
	Polyitem() {}
	Polyitem(int n, float *ptr) {num = n; memcpy(polygon, ptr, 2 * num * sizeof(float));}
	int num;
	float polygon[14];
} ;

class TableGenerator
{
private:
	std::vector<Entry>						table;
	std::vector<std::pair<float, float>>	rot_markers;
	std::vector<Polyitem>					polygons;

public:
	TableGenerator();

	void exec();

	void _calcMarkersFromLine( float d, float theta, std::vector<std::pair<float, float>>& markers );

	std::pair<float, float> _calcIntersection( float d1, float theta1, float d2, float theta2 );

	void _getFgPolygon( std::pair<float, float>&	start, 
						std::pair<float, float>&	end, 
						std::pair<float, float>&	intersec, 
						std::vector<float>&			polygon );

	int _getRotationID( std::pair<float, float>& marker );

	double _calcPixelAlpha( int x, int y, std::vector<float>& polygon );

	void visualize() 
	{
		const int num_visualize = 1000;
		const int num_cols= 30;
		int num_rows = num_visualize / num_cols;
		if (num_visualize % num_cols > 0)
			num_rows++;

		int num_patches = table.size();
		int interval = num_patches / num_visualize;

		SVGWriter writer("cpp_table.svg");
		int width = 3 * num_cols;
		int height = 3 * num_rows;
		image<rgb> *im = new image<rgb>(width, height);
		for (int id = 0; id < num_patches; id += interval) {

			int i = id / interval;
			int x = 3 * (i % num_cols);
			int y = 3 * (i / num_cols);
			rgb tmp(0, 0, 0);

			Polyitem& item = polygons[id / 2];
			std::vector<float> polygon;
			polygon.resize(2 * item.num);
			for (int j = 0; j < item.num; j++) {
				polygon[2 * j] = x + item.polygon[2 * j] / UP_SCALE;
				polygon[2 * j + 1] = y + item.polygon[2 * j + 1] / UP_SCALE;
			}
			writer.writePolygon(polygon, "#FF0000");

			tmp.g = (int)(255.0 * table[id].lo_alphas[0]);
			im->set(x, y, tmp);
			tmp.g = (int)(255.0 * table[id].lo_alphas[1]);
			im->set(x + 1, y, tmp);
			tmp.g = (int)(255.0 * table[id].lo_alphas[2]);
			im->set(x, y + 1, tmp);
			tmp.g = (int)(255.0 * table[id].lo_alphas[3]);
			im->set(x + 1, y + 1, tmp);
		}

		save("cpp_table.ppm", im);
		_saveBMP("cpp_table.bmp", im);
		writer.close();
		delete im;
	}

	void saveTable(const std::string filename)
	{
		std::ofstream file(filename, std::ios::out | std::ios::binary);
		file << table.size() << "\n";
		for (int i = 0; i < table.size(); i++) {
			for (int j = 0; j < 4; j++) 
				file << table[i].lo_alphas[j] << " ";
			file << "\n";
			for (int j = 0; j < L * L; j++)
				file << table[i].hi_alphas[j] << " ";
			file << "\n";
		}
		file.close();
	}

	void loadTable(const std::string filename)
	{
		std::ifstream file(filename, std::ios::in | std::ios::binary);
		int num_entries;
		file >> num_entries;
		table.resize(num_entries);
		for (int i = 0; i < num_entries; i++) {
			for (int j = 0; j < 4; j++)
				file >> table[i].lo_alphas[j];
			for (int j = 0; j < L * L; j++)
				file >> table[i].hi_alphas[j];
		}
		file.close();
	}



	void _saveBMP(const std::string filename, image<rgb> *im)
	{
		typedef struct {
			char id[2];
			long filesize;
			short reserved[2];
			long headersize;
			long infoSize;
			long width;
			long height;
			short biPlanes;
			short bits;
			long biCompression;
			long biSizeImage;
			long biXPelsPerMeter;
			long biYPelsPerMeter;
			long biClrUsed;
			long biClrImportant;
		} BMPHEAD;

		std::cout << sizeof(BMPHEAD) << std::endl;

		BMPHEAD bh;
		memset ((char *)&bh, 0, sizeof(BMPHEAD)); /* sets everything to 0 */
		memcpy (bh.id, "BM", 2);
			// bh.filesize  =   calculated size of your file (see below)
			// bh.reserved  = two zero bytes
			bh.headersize  = 54L;				// (for 24 bit images)
			bh.infoSize  =  40L;				// (for 24 bit images)
			bh.width     = (long)im->width();	// width in pixels of your image
			bh.height    = (long)im->height();	// depth in pixels of your image
			bh.biPlanes  =  1;					// (for 24 bit images)
			bh.bits      = 24;					// (for 24 bit images)
			bh.biCompression = 0L;				// (no compression)

		int bytesPerLine = bh.width * 3;  /* (for 24 bit images) */
		if (bytesPerLine & 0x0003) {
			bytesPerLine |= 0x0003;
			++bytesPerLine;
		}
		bh.filesize = bh.headersize + (long)bytesPerLine * bh.height;


		unsigned char *linebuf = new unsigned char[bytesPerLine];
		memset(linebuf, 0, bytesPerLine * sizeof(unsigned char));
		FILE *f = fopen(filename.c_str(), "wb");
		assert(f != NULL);

		fwrite(bh.id, 1, 2, f);
		fwrite(&bh.filesize, 1, 52, f);
		for (int y = bh.height - 1; y >= 0; y--) {
			for (int x = 0; x < bh.width; x++) {
				rgb& c = im->get(x, y);
				linebuf[3 * x + 0] = c.b;
				linebuf[3 * x + 1] = c.g;
				linebuf[3 * x + 2] = c.r;
			}
			fwrite(linebuf, 1, bytesPerLine, f);
		}
		fclose(f);
		delete linebuf;
	}
};
#pragma once

#include <string>
#include <fstream>
#include <vector>

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
			 << "stroke:" << color << "; stroke-width:0.2\"/>"
			 << std::endl;
	}

	void close()
	{
		file << "</svg>" << std::endl;
		file.close();
	}
};
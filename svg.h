#pragma once

#include <string>
#include <fstream>
#include <vector>

class SVGWriter
{
private:
	std::ofstream file;

public:
	SVGWriter(std::string filename);

	void writeDot(float x, float y, std::string color);

	void writeDots(std::vector<float> points, std::string color);

	void writePolygon(std::vector<float> points, std::string color);

	void writeText(float x, float y, std::string text, std::string color);

	void close();
};